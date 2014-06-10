/*
 * File:   snapshot.cc
 * Author: ts337
 *
 * Created on 17 June 2011, 12:43
 */

#include "coela_core/src/image_utils.h"
#include "coela_core/src/pixel_array_routines.h"
#include "coela_core/src/DS9Region.h"

#include "coela_luckypipe/src/CCD_dataset_info.h"
#include "coela_luckypipe/src/camera_config_info.h"

#include "coela_luckypipe/src/gain_utils.h"
#include "coela_luckypipe/src/single_CCD_filters.h"
#include "coela_luckypipe/src/image_cleanup.h"
#include "coela_luckypipe/src/frame_info.h"
#include "coela_analysis/src/psf_generation.h"

#include "getoptpp/getopt_pp.h"

#include <tbb/tick_count.h>
#include <tbb/task_scheduler_init.h>

#include <boost/filesystem.hpp>
using namespace std;

struct SnapshotOptions {
    int n_threads;
    int n_frames;
    string camconf_path;
    string dataset_info_path;
};

SnapshotOptions default_options()
{
    SnapshotOptions s;
    s.n_threads = 6;
    s.n_frames = 100;
    return s;
}

void print_usage_and_exit()
{
    SnapshotOptions s = default_options();
    cout<<"Usage: lucky_snapshot CCD_datasetInfo.txt \n\n"
        <<"[-n --nthreads (number)]      --- Use n processor threads (default: "
        << s.n_threads << ")"<<endl
        <<"[-m --nframes (number)]      ---Create snapshot from m frames; default:"
        << s.n_frames<< ")"<<endl;
    exit(0);
}

SnapshotOptions load_options_from_command_line(int argc, char ** argv)
{
    if (argc<2) { print_usage_and_exit(); }
    SnapshotOptions opts = default_options();
    opts.dataset_info_path =argv[1];

    using namespace GetOpt;
    GetOpt_pp input(argc, argv);
    input >> Option('n', "nthreads", opts.n_threads, opts.n_threads);
    input >> Option('m', "nframes", opts.n_frames, opts.n_frames);

    return opts;
}

using namespace std;
using namespace coela;
using namespace single_CCD_filters;
using tbb::tick_count;

int main(int argc, char** argv)
{

    SnapshotOptions opts = load_options_from_command_line(argc, argv);
    CCD_DatasetInfo dataset_inf(opts.dataset_info_path);
    CameraConfigInfo camconf(dataset_inf.default_camera_config_file);

    vector<FrameInfo> frame_vec(
        FrameInfo::get_frames_vec(dataset_inf.CCD_inputdir,
                                  dataset_inf.in_filestem,
                                  dataset_inf.extension,
                                  dataset_inf.ccd_id,
                                  dataset_inf.dataset_output_base_dir
                                 ));

    CCD_calibration_info ccd_calibration =
        camconf.get_calibration_info_for_CCD_id(dataset_inf.ccd_id);

    PixelRange crop_box = ccd_calibration.cropped_PixelRange;

    PixelArray2d<float> bias_frm(crop_box.x_dim(), crop_box.y_dim(), 0.0);

    Sequential_File_Buffer_Filter<> input_filter(frame_vec, opts.n_frames);
    Decompress_Filter decompressor;
    Frame_Count_Display_Filter count_filter(true);
    FrameCropFilter crop_filter(crop_box); //Use these two if we have a good background region

    Col_Histogram_Gather col_hist_filter(crop_box,
                                         50, crop_box.y_dim());
    col_hist_filter.set_hist_min_value(-300);

    Frame_Summation sum_filter(bias_frm.range(),
                               ccd_calibration.crop_region);

    Serial_Decommission_Filter end_filter;

    cout<<"Initializing scheduler with " << opts.n_threads <<" threads."<<endl;

    tbb::task_scheduler_init scheduler(opts.n_threads);
    tbb::pipeline pipeline;

    pipeline.add_filter(input_filter);
    pipeline.add_filter(decompressor);
    pipeline.add_filter(count_filter);
    pipeline.add_filter(crop_filter);
    pipeline.add_filter(col_hist_filter);
    pipeline.add_filter(sum_filter);
    pipeline.add_filter(end_filter);

    cout<<"Gathering histograms for column bias frame from first "<< opts.n_frames
        <<" frames..."<<endl;


    tick_count t0;
    t0 = tick_count::now();

    pipeline.run(coela::single_CCD_filters::lucky_n_tokens);

    tick_count t1 = tick_count::now();
    double time =(t1-t0).seconds();

    const vector<HistogramContainer14bit>& col_hists = col_hist_filter.get_hists();
    cout<<"\n";
    cout<<"Collected histograms of all columns in " <<opts.n_frames
        << " in "<<time<<" seconds or "<< time/ opts.n_frames
        <<" secs per frame"<<endl;

    for (size_t x=0; x!=col_hists.size(); ++x) {
        cout<<"\r Fitting col " <<x+1 <<"\t\t"<<flush;
//        col_hists[x].write_to_file( col_hists_folder+"/col"+string_utils::itoa(x)+".txt");
        double dark_estimate =
            image_cleanup::determine_histogram_bias_pedestal_via_thresholded_centroid(col_hists[x],
                    0.8);
//        double dark_estimate = gain_utils::fit_readout_gaussian(col_hists[x].convert_to_value_count_map(), false).front();
        pixel_array_routines::set_col_to_value(bias_frm, x+1, (float)dark_estimate);
    }
    cout<<"\n";


    CCDImage<double> raw_sum = sum_filter.Avg();
    string output_dir = dataset_inf.dataset_output_base_dir+"snapshot/";
    boost::filesystem::create_directories(output_dir);

    string filestem = dataset_inf.out_filestem +
                      "_CCD"+string_utils::itoa(dataset_inf.ccd_id) +
                      "_"+string_utils::itoa(opts.n_frames)+"_frm_snapshot";

//    raw_sum.write_to_file(output_dir+ filestem + "raw.fits");
    CCDImage<double> coldb = raw_sum;
    coldb.pix -= PixelArray2d<double>(bias_frm);
    coldb.write_to_file(output_dir+ filestem + ".fits");

    return 0;
}

