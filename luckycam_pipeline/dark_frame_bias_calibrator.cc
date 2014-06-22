
#include "coela_utility/src/misc_math.h"

#include "coela_core/src/ccd_image.h"
#include "coela_core/src/image_utils.h"
#include "coela_core/src/ds9_region.h"

#include "coela_luckypipe/src/ccd_dataset_info.h"
#include "coela_luckypipe/src/camera_config_info.h"

#include "coela_luckypipe/src/gain_utils.h"
#include "coela_luckypipe/src/single_ccd_filters.h"
#include "coela_luckypipe/src/image_cleanup.h"
#include "coela_luckypipe/src/frame_info.h"
#include "coela_luckypipe/src/clean_and_register_subroutines.h"
#include "coela_analysis/src/psf_generation.h"

//#include "../threading/threading_tools.h"
#include <tbb/tick_count.h>
#include <tbb/task_scheduler_init.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "getoptpp/getopt_pp.h"


using namespace coela;
using namespace std;
using tbb::tick_count;

struct ProgramOptionSet {
    string dataset_info_filename;
    CcdDatasetInfo dataset;
};

void print_usage_and_exit()
{
    exit(0);
}

ProgramOptionSet load_options_from_command_line(int argc, char** argv)
{
    if (argc != 2) { print_usage_and_exit(); }
    ProgramOptionSet opts;
    opts.dataset_info_filename=argv[1];
    opts.dataset = CcdDatasetInfo(opts.dataset_info_filename);
    return opts;
}

struct RowHists10bit {
    RowHists10bit(size_t width, int min_val) {
        pixel_hists.resize(width);
        for (size_t i=0; i!=pixel_hists.size(); ++i) {
            pixel_hists[i].set_min_value(min_val);
        }
    }

    void clear() {
        for (size_t i=0; i!=pixel_hists.size(); ++i) {
            pixel_hists[i].clear();
        }
    }
    vector<HistogramContainer10bit> pixel_hists;
};

void process_n_rows(vector<FrameInfo>& frame_vec,
                    PixelArray2d<float>& bias_frm,
                    vector<RowHists10bit>& row_hist_vec,
                    const PixelRange& crop_box,
                    const PixelRange& temporal_debias_box,
                    size_t first_row,
                    size_t n_rows_to_process,
                    size_t n_frames_to_process
                   )
{
    if ((first_row + n_rows_to_process -1) > bias_frm.range().y_dim()) {
        throw runtime_error("!!!");
    }

    for (size_t i=0; i!=row_hist_vec.size(); ++i) {
        row_hist_vec[i].clear();
    }
    cout << n_rows_to_process << "rows starting at " << first_row << endl;

    HistogramContainer14bit single_frame_hist;

    for (size_t i=0; i!=n_frames_to_process; ++i) {
        cout<<"\rFrame "<<i<<" of "<< n_frames_to_process<<flush;
        FrameInfo& frm= frame_vec[i];
        CcdImage<float> img;

        img.pix = PixelArray2d<float>::load_from_file(frm.file_path,
                   frm.header_byte_offset);
        img = CcdImage<float>::sub_image(img, crop_box);

        single_frame_hist.clear();
        if (frm.bias_pedestal== 0.0) {
            frm.bias_pedestal =
                image_cleanup::determine_bias_pedestal_from_box_in_raw_image(
                    img.pix,
                    temporal_debias_box,
                    single_frame_hist);
        }


        img.pix -= frm.bias_pedestal;

        for (size_t offset = 0; offset!= n_rows_to_process; ++offset) {
            size_t row_num = first_row + offset;
            RowHists10bit& rowhist = row_hist_vec[offset];
            for (size_t x=1; x<=img.pix.range().x_dim(); ++x) {
                rowhist.pixel_hists[x-1].count_rounded_value(img.pix(x, row_num));
            }
        }
    }
    cout<<endl;

    for (size_t offset = 0; offset!= n_rows_to_process; ++offset) {
        size_t row_num = first_row + offset;

        for (size_t x=1; x<=bias_frm.range().x_dim(); ++x) {
            bias_frm(x, row_num) =
                image_cleanup::determine_histogram_bias_pedestal_via_thresholded_centroid(
                    row_hist_vec[offset].pixel_hists[x-1]);
        }
    }
}


int main(int argc, char** argv)
{

    ProgramOptionSet opts = load_options_from_command_line(argc, argv) ;
    const CcdDatasetInfo& cdi(opts.dataset);
    cout<<"CCD output dir: "<<cdi.CCD_outputdir<<endl;

    CameraConfigInfo camconf(cdi.default_camera_config_file);

    vector<FrameInfo> frame_vec(
        FrameInfo::get_frames_vec(cdi.CCD_inputdir,
                                  cdi.in_filestem,
                                  cdi.extension,
                                  cdi.ccd_id,
                                  cdi.dataset_output_base_dir));

    FrameInfo::set_bias_pedestal_estimates(frame_vec, 0.0);

    string output_folder=cdi.CCD_outputdir+"/calibration_data/";
    boost::filesystem::create_directories(output_folder);

    const size_t n_frames_to_process_limit = 3000;
    const int n_frames_to_process = min(n_frames_to_process_limit, frame_vec.size());

    size_t n_rows_to_process_at_once = 500;

    CcdImage<float> first;
    FrameInfo& frm= frame_vec.front();




    first.pix = PixelArray2d<float>::load_from_file(frm.file_path, frm.header_byte_offset);
    first.initialize_CCD_grid_for_raw_data();
    PixelRange crop_box = camconf.get_calibration_info_for_CCD_id(
                              frm.ccd_id).cropped_PixelRange;
    first = CcdImage<float>::sub_image(first, crop_box);

    CcdBoxRegion temporal_debias_rgn =
        camconf.get_calibration_info_for_CCD_id(
            frm.ccd_id).default_temporal_debiasing_histogram_region;

    PixelRange temporal_debias_box =
        first.CCD_grid.corresponding_pixel_region(temporal_debias_rgn).bounded_pixels();


    //============================================================================================
    //Per pixel bias estimation:
    cerr<<"Allocating pixel hists...";
    vector<RowHists10bit> row_hist_vec;


    for (size_t i=0; i!=n_rows_to_process_at_once; ++i) {
        row_hist_vec.push_back(RowHists10bit(first.pix.range().x_dim(), -300));
    }

    cerr<<"Done"<<endl;
    cout<<endl;

    CcdImage<float> per_pixel_bias_estimates(first);
    per_pixel_bias_estimates.pix.assign(0);
    //----------------------------------------------------------------------------------
    size_t first_row=1;
    while (n_rows_to_process_at_once!=0) {
        if (first_row + n_rows_to_process_at_once -1 <= first.pix.range().y_dim()) {
            process_n_rows(frame_vec,
                           per_pixel_bias_estimates.pix,
                           row_hist_vec,
                           crop_box,
                           temporal_debias_box,
                           first_row,
                           n_rows_to_process_at_once,
                           n_frames_to_process);

            per_pixel_bias_estimates.write_to_file(output_folder+"per_pixel_bias_frame.fits");
            first_row+=n_rows_to_process_at_once;

        } else {
            n_rows_to_process_at_once--;
        }

    }

    //============================================================================================
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------

    FrameInfo::write_list_to_file(frame_vec, output_folder+"frame_list.txt");

    double uniform_bias_est =
        FrameInfo::get_median_bias_pedestal_estimate(frame_vec);
    cout<<"Uniform bias est: "<< uniform_bias_est<<endl;

//
    CcdImage<float> standard_col_bias_frm =
        clean_and_register_subroutines::
        estimate_column_bias_pattern(
            frame_vec, output_folder,
            crop_box,
            false,
            temporal_debias_box,
            n_frames_to_process,
            7);
//
    standard_col_bias_frm.write_to_file(output_folder +"std_col_bias_frame.fits");

    CcdImage<float> row_bias_frm =
        clean_and_register_subroutines::
        estimate_row_bias_pattern(
            frame_vec, output_folder,
            crop_box,
            standard_col_bias_frm,
            n_frames_to_process,
            7);
//
    row_bias_frm.write_to_file(output_folder +"row_bias_frame.fits");

//----------------------------------------------------------------------------------
    CcdImage<float> blank_frame(per_pixel_bias_estimates);
    blank_frame.pix.assign(0.0);

    CcdImage<double> raw_avg =
        clean_and_register_subroutines::create_debiased_average(
            frame_vec,
            crop_box,
            first.CCD_grid.image_outline_,
            blank_frame,
            n_frames_to_process,
            7);

    CcdImage<double>pixel_debiased_avg(raw_avg);
    pixel_debiased_avg.pix -= PixelArray2d<double>(per_pixel_bias_estimates.pix);
    pixel_debiased_avg.write_to_file(output_folder+"pixel_debiased_avg.fits");

    CcdImage<float> combined_bias_frame(standard_col_bias_frm);
    combined_bias_frame.pix += row_bias_frm.pix;

    CcdImage<double> algorithm_debiased_avg(raw_avg);
//        clean_and_register_subroutines::create_debiased_average(
//            frame_vec,
//            crop_box,
//            first.CCD_grid.image_outline_,
//            combined_bias_frame,
//            n_frames_to_process,
//            7);

    algorithm_debiased_avg.pix -= PixelArray2d<double>(combined_bias_frame.pix);
    algorithm_debiased_avg.write_to_file(output_folder+"algorithm_debiased_avg.fits");

//----------------------------------------------------------------------------------
    CcdImage<float> diff_after_col_debias(per_pixel_bias_estimates);
    diff_after_col_debias.pix -= standard_col_bias_frm.pix;
    diff_after_col_debias.write_to_file(output_folder+"col_db_bias_diff.fits");

    CcdImage<float> diff_after_col_row_db(diff_after_col_debias);
    diff_after_col_row_db.pix -= row_bias_frm.pix;
    diff_after_col_row_db.write_to_file(output_folder + "col_and_row_db_bias_diff.fits");

    return (0);
}
