/*
 * File:   main.cc
 * Author: ts337
 *
 * Created on 15 March 2011, 15:27
 */

#include "coela_utility/src/misc_math.h"

#include "coela_core/src/image_utils.h"
#include "coela_core/src/ds9_region.h"
#include "coela_core/src/drizzle.h"

#include "coela_luckypipe/src/ccd_dataset_info.h"
#include "coela_luckypipe/src/camera_config_info.h"

#include "coela_luckypipe/src/gain_utils.h"
#include "coela_luckypipe/src/image_cleanup.h"
#include "coela_luckypipe/src/frame_info.h"

#include "coela_luckypipe/src/single_ccd_filters.h"
#include "coela_luckypipe/src/multi_ccd_filters.h"

#include "coela_luckypipe/src/drizzle_subroutines.h"

//#include "../threading/threading_tools.h"
#include <tbb/tick_count.h>
#include <tbb/task_scheduler_init.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "getoptpp/getopt_pp.h"
#include <map>

using namespace std;
using namespace coela;
using namespace drizzle;
using namespace pipeline;
using namespace multi_CCD_filters;
using namespace GetOpt;
using namespace coela::string_utils;
using tbb::tick_count;

///This struct is used with the intention of eventually enabling reduction across multiple observations of the same target.
struct ProcessedDatasetInfo {
    CCD_DatasetInfo info;
    string info_file_path;
    string frames_list_path;
    vector<FrameInfo> frames;
};

struct ProgramOptionSet {
    int n_threads;

    CameraConfigInfo camconf;
    DrizzleSettings drizzle_settings;
    vector<ProcessedDatasetInfo> guide_CCD_datasets;
    vector<ProcessedDatasetInfo> additional_CCD_datasets;


    bool separate_percentile_output;
    vector< pair<double,double> > output_percentile_pairs;

//     vector<string> slaved_CCD_frame_list_filenames;
    bool output_diagnostic_files;
    bool output_frame_by_frame_mosaics;  //for analysis of wide field tip-tilt anisoplanatism etc.
//     string cosmic_ray_output_folder;
//     bool quick_zeroing;


    bool autoguide_sim;
    bool pixel_average_mode;


};



ProgramOptionSet default_drizzle_options()
{

    ProgramOptionSet default_options;

    default_options.n_threads=7;
//    default_options.output_diagnostic_files=true;

    default_options.output_frame_by_frame_mosaics = false;
//    default_options.cosmic_ray_output_folder="drizzle_diagnostic_frames";

    default_options.separate_percentile_output=false;

    default_options.autoguide_sim=false;
    default_options.pixel_average_mode=false;

    return default_options;
}

void print_usage_and_exit()
{
    cerr<<"Typical Usage: \nlucky_drizzle -i /path/to/DrizzleSettings.txt -g dataset_info.txt \n"
        <<"[-n --gm -a --am --autoguider --bands -D -O]\n\n"
        <<"-g  /path/to/CCD_dataset_info.txt \n"
        <<"     (Guide star information files, may have multiple entries)\n"
        <<"--gm /path/to/CCD_dataset_info.txt /path/to/corresponding_frame_list.txt\n"
//            <<"     (Guide star information files paired with manually chosen frame lists)\n"
        <<"-a  path/to/additional_ccd_info_1.txt \n"
        << "        (drizzle additional CCDs (shifts same as main CCD) )\n"
//            <<"--am  path/to/additional_ccd_info_1.txt path/to/corresponding_frame_list.txt \n"
//            << "        (drizzle additional CCDs with non-default frame_lists_\n"
        <<"--pixelaverage   (Create drizzled average without any pixel shifts - used for removing cosmic rays when creating dark current frames, etc)"<<"\n"
        <<"--autoguider (autoguider sim, simple output only.) \n"
        <<"--bands  L1 U1 L2 U2 ....  (output percentile bands as separate, rather than cumulative, images.) \n"
        <<"-n --nthreads 4  (run the pipeline with m threads (default 6)) \n"
        <<"-D (true/false)     (Output files marked for further inspection (default: true))\n"
        <<"-O     (Do not perform drizzle, instead output mosaics of the individual frames (after debiasing))\n"
        <<endl;
    exit(0);
}

ProgramOptionSet get_options_from_command_line(int argc, char** argv)
{
    if (argc<3) { print_usage_and_exit(); }
    ProgramOptionSet default_opts = default_drizzle_options();
    ProgramOptionSet opts = default_opts;


    GetOpt_pp cmd_line_args(argc, argv);
    cmd_line_args.exceptions_all(); //turn on exceptions

    string drizzle_settings_path;
    cmd_line_args >> Option('i', drizzle_settings_path);
    opts.drizzle_settings = DrizzleSettings(drizzle_settings_path);

    if (opts.drizzle_settings.normalisation_on==false) {
        opts.drizzle_settings.thresholding_on=false;
    }

    cmd_line_args >> Option('n', "nthreads", opts.n_threads, default_opts.n_threads);

    //----------------------------------------------------------------------
    //Get guide star datasets:
    if (cmd_line_args>>OptionPresent('g')) {
        vector<string> guide_star_dataset_info_paths;
        cmd_line_args >> Option('g', guide_star_dataset_info_paths);
        for (size_t i=0; i!=guide_star_dataset_info_paths.size(); i++) {
            ProcessedDatasetInfo pdi;
            pdi.info_file_path = guide_star_dataset_info_paths[i];
            pdi.info = CCD_DatasetInfo(pdi.info_file_path);
            pdi.frames_list_path = pdi.info.most_recently_output_frame_list;
            cout<<"Attempting to load guide star frame list from \""<<pdi.frames_list_path<<"\""<<endl;
            pdi.frames = FrameInfo::load_frames_vec(pdi.frames_list_path);
            opts.guide_CCD_datasets.push_back(pdi);
        }
    } else if (cmd_line_args>>OptionPresent("gm")) {
        vector<string> gs_info_list_pairs;
        cmd_line_args >> Option("gm", gs_info_list_pairs);
        for (size_t i=0; i<gs_info_list_pairs.size()-1; i+=2) {
            ProcessedDatasetInfo pdi;
            pdi.info_file_path = gs_info_list_pairs[i];
            pdi.frames_list_path = gs_info_list_pairs[i+1];

            pdi.info = CCD_DatasetInfo(pdi.info_file_path);
            pdi.frames = FrameInfo::load_frames_vec(pdi.frames_list_path);
            opts.guide_CCD_datasets.push_back(pdi);
        }
    } else {
        throw runtime_error("No guide star / CCD specified");
    }



    assert(opts.guide_CCD_datasets.size()==
           1); //Fixme - currently only works with 1 dataset at a time
    //----------------------------------------------------------------------

    if (cmd_line_args >> OptionPresent('a'))  {
        vector<string> additional_dataset_info_paths;
        cmd_line_args >> Option('a', additional_dataset_info_paths);
        for (size_t i=0; i!=additional_dataset_info_paths.size(); i++) {
            if (additional_dataset_info_paths[i]!= opts.guide_CCD_datasets.front().info_file_path) {
                //Will need to make this more sophisticated to deal with multiple obs.
                ProcessedDatasetInfo pdi;
                pdi.info = CCD_DatasetInfo(additional_dataset_info_paths[i]);
                pdi.frames_list_path = pdi.info.most_recently_output_frame_list;
                pdi.frames = FrameInfo::load_frames_vec(pdi.frames_list_path);
                opts.additional_CCD_datasets.push_back(pdi);
            }
        }
    }

    if (cmd_line_args >> OptionPresent('O'))  {
        opts.output_frame_by_frame_mosaics=true;
    }

//    string output_diagnostic_files_str;
//    ops >> Option('D', "output_diagnostics", output_diagnostic_files_str,
//        string_utils::bool_to_string(default_opts.output_diagnostic_files));
//    output_opts.output_diagnostic_files=
//            string_utils::string_to_bool(output_diagnostic_files_str);

    if (cmd_line_args >> OptionPresent("autoguider")) { opts.autoguide_sim=true; }

    else if (cmd_line_args >> OptionPresent("pixelaverage")) {
        opts.pixel_average_mode=true;
        if (!opts.additional_CCD_datasets.empty()) {
            throw runtime_error("PixelAverage mode is only intended to be used with a single dataset at a time");
        }
    }
    //If doing autoguider sim or simple average, just use all the frames:
    //Otherwise, can specify quality banded output,
    //Or the most common option, which is to produce cumulative images at the levels specified in drizzle settings.

    if (opts.autoguide_sim || opts.pixel_average_mode) {
        opts.output_percentile_pairs.clear();
        opts.output_percentile_pairs.push_back(pair<double,double>(0.0,1.0));
    } else if (cmd_line_args >> OptionPresent("bands")) {
        opts.separate_percentile_output=true;

        vector<string> percentile_pair_strings;
        cmd_line_args >> Option("bands", percentile_pair_strings);

        if (percentile_pair_strings.size() % 2 !=0) { throw runtime_error("Check percentile pairings"); }
        for (size_t i=0; i<percentile_pair_strings.size(); i+=2) {
            double lower = atof(percentile_pair_strings[i]);
            double upper = atof(percentile_pair_strings[i+1]);
            if (lower>=upper) { throw runtime_error("Check percentile pairings"); }
            opts.output_percentile_pairs.push_back(pair<double,double>(lower,upper));
        }
    } else {
        double upper_limit=1.00;
        for (size_t i=0; i!=opts.drizzle_settings.output_percentiles.size(); ++i) {
            double lower_limit = 1.0 - opts.drizzle_settings.output_percentiles[i]/100.0;
            opts.output_percentile_pairs.push_back(pair<double,double>(lower_limit, upper_limit));
            upper_limit=lower_limit;
        }
    }

    string camconf_path = opts.guide_CCD_datasets.front().info.default_camera_config_file;
    opts.camconf=CameraConfigInfo(camconf_path);
    cout<<"Frame lists for guide star CCD:  "<< endl;
    for (size_t i=0; i!=opts.guide_CCD_datasets.size(); ++i) {
        cout << opts.guide_CCD_datasets[i].info.most_recently_output_frame_list << endl;
        if (opts.guide_CCD_datasets[i].info.default_camera_config_file!=camconf_path) {
            throw runtime_error("Mismatching camera config path for "+
                                opts.guide_CCD_datasets[i].info.CCD_outputdir);
        }
    }
    cout<<"Additional CCD frame lists: "<<endl;
    for (size_t i=0; i!=opts.additional_CCD_datasets.size(); ++i) {
        cout << opts.additional_CCD_datasets[i].info.most_recently_output_frame_list << endl;
        if (opts.additional_CCD_datasets[i].info.default_camera_config_file!=camconf_path) {
            throw runtime_error("Mismatching camera config path for "+
                                opts.additional_CCD_datasets[i].info.CCD_outputdir);
        }
    }

//    cmd_line_args.end_of_options();

    return opts;
}

std::vector<MultiFrame> get_multi_frame_vec(const std::vector<FrameInfo>& gs_frm_vec,
        const std::vector< ProcessedDatasetInfo >& additional_datasets,
        bool perform_timestamp_correlation_and_correction)
{
    vector< vector<FrameInfo> > additional_frames;
    for (size_t i=0; i!=additional_datasets.size(); ++i) {
        additional_frames.push_back(additional_datasets[i].frames);
    }
    return MultiFrame::get_multi_frame_vec(gs_frm_vec,
                                           additional_frames,
                                           perform_timestamp_correlation_and_correction);
}

int main(int argc, char** argv)
{
    get_options_from_command_line(argc, argv);

    cout <<"\n***** Single GS Drizzle MT ***** " <<endl;

    ProgramOptionSet opts=get_options_from_command_line(argc, argv);

//==========================================================================================================
//*********** General Prep *********************
//==========================================================================================================

    //==================== Load the additional frame lists, unify lists =================================
    assert(opts.guide_CCD_datasets.size()==
           1);    //FIXME - currently only works for a single observation run.
    vector<FrameInfo>& gs_CCD_frms = opts.guide_CCD_datasets.front().frames;

    vector<MultiFrame> multi_frm_vec =
        get_multi_frame_vec(gs_CCD_frms,
                            opts.additional_CCD_datasets,
                            opts.camconf.simulated_data);


    //Set to the collated list:
    gs_CCD_frms = MultiFrame::retrieve_frames_for_CCD_id(
                      gs_CCD_frms.front().ccd_id,
                      multi_frm_vec);

    if (opts.pixel_average_mode) {
        FrameInfo::set_blank_guide_star_positions(gs_CCD_frms);
    }

    //==================== //Add percentile ranking =================================

    FrameInfo::add_percentile_rankings_based_on_first_star(gs_CCD_frms);

    const DrizzleSettings& ds(opts.drizzle_settings);

    //=================== Create directory structure, output copies of unified frame lists =============================
    string drizzle_output_dir = ds.create_output_dirs(
                                    opts.guide_CCD_datasets.front().frames_list_path, gs_CCD_frms, multi_frm_vec);

    //================ Check a few settings, dump some sample frames: =============================
    cout <<"First file path " << gs_CCD_frms.front().file_path<<endl;
    cout<<"==========================================================="<<endl;
    cout<<"Normalisation on: " << bool_to_string(opts.drizzle_settings.normalisation_on)
        <<endl;
    cout<<"Thresholding on: " << bool_to_string(opts.drizzle_settings.thresholding_on)<<endl;



    MosaicImage<float>first_gs_frm(
        CCDImage<float>(gs_CCD_frms.front().file_path,
                        gs_CCD_frms.front().header_byte_offset));
    first_gs_frm.initialize_CCD_grid_for_raw_data();

    const CCD_calibration_info gs_ccd_properties
        =opts.camconf.get_calibration_info_for_CCD_id(gs_CCD_frms.front().ccd_id);

    if (!opts.camconf.simulated_data) {
        first_gs_frm = CCDImage<float>::sub_image(first_gs_frm,
                       gs_ccd_properties.cropped_PixelRange);
    }
//    p.first_gs_frm.add_comment("Composed from run "+(string_utils::pull_filestem(p.first_gs_frm.associated_filename())));
    first_gs_frm.write_to_file(drizzle_output_dir+"first_GS_frame_raw_cropped.fits");

    //=================  Prep CCD specific drizzle info structs: =================================================================

    vector<CCD_DatasetInfo> all_datasets_info;

    for (size_t i=0; i!=opts.guide_CCD_datasets.size(); ++i) {
        all_datasets_info.push_back(opts.guide_CCD_datasets[i].info);
    }
    for (size_t i=0; i!=opts.additional_CCD_datasets.size(); ++i) {
        all_datasets_info.push_back(opts.additional_CCD_datasets[i].info);
    }
    vector<CCD_specific_drizzle_inf> drizzle_inf_vec =
        CCD_specific_drizzle_inf::prep_drizzle_data_vec(ds,
                all_datasets_info,
                opts.camconf, multi_frm_vec.front(), drizzle_output_dir);

    if (opts.pixel_average_mode) {
        assert(drizzle_inf_vec.size()==1);
        CCD_specific_drizzle_inf& ccd_inf = drizzle_inf_vec.front();
        ccd_inf.output_mosaic_region = ccd_inf.input_mosaic_region;
        ccd_inf.drizzled_output_PixelRange = first_gs_frm.pix.range();
        ccd_inf.output_CCD_low_corner = ccd_inf.input_CCD_region.low;

    }

    //=================  Get target posn for guide star =============================
    CCD_Position gs_target_CCD_posn =FrameInfo::mean_guide_star_estimate(
                                         gs_CCD_frms).Position;
    cout<<"GS target CCD posn: " <<gs_target_CCD_posn<<endl;

//   ================== Figure out the output region mosaic coords ==================================
    CCD_specific_drizzle_inf gs_drizzle_inf =
        CCD_specific_drizzle_inf::find_drizzle_info_for_ccd_id(gs_CCD_frms.front().ccd_id,
                drizzle_inf_vec);
    first_gs_frm.initialize_mosaic_grid_to_specific_region(
        gs_drizzle_inf.input_mosaic_region);

//==========================================================================================================
//*********** Create Pipeline *********************
//==========================================================================================================
//
    cout <<"Creating pipeline"<<endl;
//
    tbb::task_scheduler_init init_manual(opts.n_threads);
////    tbb::task_scheduler_init init_auto;
////    tbb::task_scheduler_init init_single(1);
    tbb::pipeline pipeline;
//
//
//    if (ds.guide_regions.empty()) ds.guide_regions.push_back(multi_ccd_region());
//
    Quality_Selected_Load_to_Buffer_Filter input_buffer(
        gs_CCD_frms,
        multi_frm_vec,
        drizzle_inf_vec,
        gs_target_CCD_posn
    );

    single_CCD_filters::Decompress_Filter decompressor;

    single_CCD_filters::Frame_Count_Display_Filter count_filter(true);

    pipeline.add_filter(input_buffer);
    pipeline.add_filter(decompressor);
    pipeline.add_filter(count_filter);



    Drizzle_Token_Crop_and_Debias_Filter db_filter;

    Cosmic_Ray_Downweighting_Filter ray_weight_filter(20);

//    single_CCD_filters::Float_To_Double_Filter sim_data_convert_filter;

    if (!opts.camconf.simulated_data) {
        pipeline.add_filter(db_filter);
        pipeline.add_filter(ray_weight_filter);
    } else {
//        pipeline.add_filter(sim_data_convert_filter);
    }

    normalization_filter norm_filter;

    const double threshold_in_photo_electrons=0.33;
    Photon_Thresholding_Filter thresh_filter(threshold_in_photo_electrons);
//    //amplifier glow subtraction
    Dark_Current_Subtraction_Filter DC_filter(ds.thresholding_on);
    if (ds.normalisation_on) { pipeline.add_filter(norm_filter); }
    if (ds.normalisation_on && ds.thresholding_on) { pipeline.add_filter(thresh_filter); }
    if (ds.normalisation_on) { pipeline.add_filter(DC_filter); }
//
//    //Parallel drizzle filters
    drizzle_filter drizzler(ds.drizzle_scale_factor,
                            ds.drizzle_pixel_fraction,
                            ds.thresholding_on);


//    raw_data_mosaic_filter raw_mosaic_output_filter(drizzle_inf_vec,
//            ds.output_base_folder+"raw_mosaics/",
//             full_output_region);

    if (!opts.output_frame_by_frame_mosaics) { pipeline.add_filter(drizzler); }
//    else if (opts.output_frame_by_frame_mosaics) pipeline.add_filter(raw_mosaic_output_filter);

//    //Rapid summation filter, creates all appropriate output frames and chooses the right one based on CCD_id
    output_summation_filter sum_filter(drizzle_inf_vec, ds.thresholding_on);
    if (!opts.output_frame_by_frame_mosaics) { pipeline.add_filter(sum_filter); }
//
    single_CCD_filters::Serial_Decommission_Filter end_filter;
    pipeline.add_filter(end_filter);
//

    tick_count t0 = tick_count::now();
    cout<<"Drizzling from "<<gs_CCD_frms.size()<<" exposures over "<<
        multi_frm_vec.front().synchronized_CCD_frames.size()<<" CCDs ("
        << multi_frm_vec.front().synchronized_CCD_frames.size()*gs_CCD_frms.size()
        <<" frames)"<<endl;
    cout<<"Running with "<< opts.n_threads<< " theads"<<endl;

    if (opts.output_frame_by_frame_mosaics) {
        opts.output_percentile_pairs.clear();
        opts.output_percentile_pairs.push_back(pair<double,double>(0.0,1.0));
    }

    for (size_t percentiles_index=0; percentiles_index!=opts.output_percentile_pairs.size();
            ++percentiles_index) {
        if (opts.separate_percentile_output) { sum_filter.reset(); }
        input_buffer.set_selection_limits(opts.output_percentile_pairs[percentiles_index].first,
                                          opts.output_percentile_pairs[percentiles_index].second);
        cout <<"Running pipeline with limits "<<input_buffer.lower_limit()
             <<", "<<input_buffer.upper_limit()<<endl;

        //******************************************      PIPELINE READY...   *********************************************


        //*********************************************** RUN PIPELINE **************************
        pipeline.run(single_CCD_filters::lucky_n_tokens);


        //*********************************************** PROCESS RESULTS...     **************************
        tick_count t1 = tick_count::now();
        float time =(t1-t0).seconds();

        int counter= count_filter.n_counted();
        cout<<"Drizzled " << counter <<" frames in " << time << " seconds, a rate of " <<
            counter/time <<" fps or " << time/counter<<" secs/frame"<<endl;
        cout<<"Compressed Data throughput rate is approx"<<input_buffer.bytes_loaded()/1e6/time
            <<"MB /sec"<<endl;

        vector< MosaicImage<double> >
        CCD_drizzle_sums; //a convenient copy for the "Append" mosaic creation
        vector< MosaicImage<double> > CCD_drizzle_weights;
        vector< MosaicImage<double> > threshed_CCD_drizzle_sums;
        //    string filename_percentile_text = ftoa(100.0-input_buffer.upper_limit())+"p_"+ftoa(100.0-input_buffer.lower_limit())+"p"; //separate
        string filename_percentile_text = ftoa(input_buffer.lower_limit(),
                                               3) +"_1.0_pC";  //cumulative
        if (opts.separate_percentile_output) {
            filename_percentile_text =ftoa(input_buffer.lower_limit(),
                                           3) +"_"+ftoa(input_buffer.upper_limit(),3)+"pS";
        }


        //******************* Unweight, col-debias and save to file the individual CCD drizzle images ****************************
        for (size_t i=0; i!=drizzle_inf_vec.size(); ++i) {
            int current_CCD_id=drizzle_inf_vec[i].CCD_id;


            CCD_drizzle_sums.push_back(sum_filter.drizzled_vals_for_CCD(current_CCD_id));
            CCD_drizzle_weights.push_back(sum_filter.drizzled_weights_for_CCD(current_CCD_id));

            FitsHeader additional_hdr_info_per_CCD;
            additional_hdr_info_per_CCD.add_keyword("CCD_ID", itoa(current_CCD_id), "CCD ID number");

            MosaicImage<double>  analog_CCD_image;
            analog_CCD_image.pix = drizzle::unweight_drizzle_results(CCD_drizzle_sums.back().pix,
                                   CCD_drizzle_weights.back().pix);
            analog_CCD_image.initialize_CCD_grid_to_specific_offset_and_scale(
                drizzle_inf_vec[i].output_CCD_low_corner, ds.drizzle_scale_factor);
            analog_CCD_image.initialize_mosaic_grid_to_specific_region(
                drizzle_inf_vec[i].output_mosaic_region);
            //            if (analog_CCD_image.CCD_id_from_keyval()==-1) analog_CCD_image.set_keyword("CAMERAID", itoa(current_CCD_id));
            analog_CCD_image.write_to_file(
                drizzle_output_dir+"/CCD"+itoa(current_CCD_id)+"/"+filename_percentile_text+"_AM.fits" ,
                additional_hdr_info_per_CCD);

            if (ds.post_drizzle_col_debias_on && !opts.camconf.simulated_data) {
                //FIXME
//                cerr<<"Will not do post drizzle col debias - FIXME!"<<endl;
                image_cleanup::col_10p_debias_preserving_bg_level(&analog_CCD_image.pix,100);
                analog_CCD_image.write_to_file(drizzle_output_dir+"/CCD"+itoa(current_CCD_id)+"/"
                                               +filename_percentile_text+"_AM_col-db.fits");
                image_cleanup::row_10p_debias_preserving_bg_level(&analog_CCD_image.pix);
                analog_CCD_image.write_to_file(drizzle_output_dir+"/CCD"+itoa(current_CCD_id)+"/"
                                               +filename_percentile_text+"_AM_col-row-db.fits");

                analog_CCD_image.pix *= CCD_drizzle_weights.back().pix;
                CCD_drizzle_sums.back() =
                    analog_CCD_image; //replace the weighted CCDImage with its column debiased counterpart
            }
            CCD_drizzle_weights.back().write_to_file(drizzle_output_dir+"/CCD"+itoa(
                        current_CCD_id)+"/"+filename_percentile_text+"_weights.fits");

            if (ds.thresholding_on) {
                threshed_CCD_drizzle_sums.push_back(sum_filter.drizzled_threshed_vals_for_CCD(
                                                        current_CCD_id));
                CCDImage<double>  threshed_CCD_image;
                threshed_CCD_image.pix = drizzle::unweight_drizzle_results(
                                             threshed_CCD_drizzle_sums.back().pix, CCD_drizzle_weights.back().pix);
                threshed_CCD_image.write_to_file(drizzle_output_dir+"/CCD"+itoa(
                                                     current_CCD_id)+"/"+filename_percentile_text+"_PC.fits");

                if (ds.post_drizzle_col_debias_on) {
                    //FIXME
                    image_cleanup::col_10p_debias_preserving_bg_level(&threshed_CCD_image.pix,100);
                    threshed_CCD_image.write_to_file(drizzle_output_dir+"/CCD"+itoa(
                                                         current_CCD_id)+"/"+filename_percentile_text+"_PC_col-db.fits");
                    threshed_CCD_drizzle_sums.back() = threshed_CCD_image;
                    threshed_CCD_drizzle_sums.back().pix *=
                        CCD_drizzle_weights.back().pix; //replace the weighted CCDImage with its column debiased counterpart
                }

                //To do:
                //                psf_models::reference_psf combined_image = combine_regular_and_threshed_bitmaps(
                //                        analog_CCD_image, threshed_CCD_image,
                //                        1.0/current_ccd_gain_inf.calculate_threshold_pass_fraction(),
                //                        0.20,
                //                        ds.drizzle_zoom_factor*1.5);
                ////                if (ds.post_drizzle_col_debias_on) image_cleanup::col_10p_debias(combined_image);
                //                combined_image.write_to_file(drizzle_output_dir+"/CCD"+itoa(current_CCD_id)+"/"+filename_percentile_text+"_combined.fits" );
                //                combined_image.write_mask_to_file(drizzle_output_dir+"/CCD"+itoa(current_CCD_id)+"/"+filename_percentile_text+"_combined_PC_mask.fits" );
            }

        }


        if (ds.create_drizzle_mosaic && !opts.additional_CCD_datasets.empty()) {
            MosaicImage<double> mosaic_sum = init_blank_mosaic_frame(
                                                 CCD_specific_drizzle_inf::get_mosaic_region(drizzle_inf_vec),
                                                 opts.drizzle_settings.drizzle_scale_factor);
            ImageGrid<coordinate_types::mosaic> mosaic_ref_frame(
                CCD_specific_drizzle_inf::get_mosaic_region(drizzle_inf_vec),
                opts.drizzle_settings.drizzle_scale_factor);
            MosaicImage<double> mosaic_weights(mosaic_sum);

            create_quick_drizzle_mosaic(CCD_drizzle_sums,
                                        CCD_drizzle_weights,
                                        mosaic_sum,
                                        mosaic_weights);

            MosaicImage<double> unweighted_mosaic(mosaic_sum);
            unweighted_mosaic.pix =
                drizzle::unweight_drizzle_results(mosaic_sum.pix, mosaic_weights.pix);

            unweighted_mosaic.write_to_file(
                drizzle_output_dir+filename_percentile_text+"_AM.fits");

            //================== Save a target marker to the output dir ===============================    //
            MosaicPosition gs_mos_posn = first_gs_frm.mosaic_grid.corresponding_grid_Position(
                                             first_gs_frm.CCD_grid.corresponding_pixel_Position(gs_target_CCD_posn));

            vector<ds9::DS9Region> gs_target_file;
            gs_target_file.push_back(ds9::DS9Region::point(gs_mos_posn));
            ds9::DS9Region::save_regions_to_file(
                drizzle_output_dir+"drizzle_targets.reg", gs_target_file);
            //============================================================================
            //
            //
            //            CCDImage<double>  unweighted_full_image(unweight_drizzle_results(full_output_sum, full_output_weights));
            ////            if (ds.normalisation_on) unweighted_full_image.add_normalisation_flag();
            //            unweighted_full_image.write_to_file(drizzle_output_dir+filename_percentile_text+"_AM.fits");
            //            full_output_weights.write_to_file(drizzle_output_dir+filename_percentile_text+"_weights.fits");
            //
            //            if (ds.thresholding_on){
            //                create_drizzle_mosaic(threshed_CCD_drizzle_sums, CCD_drizzle_weights, full_threshed_output, full_output_weights);
            //                CCDImage<double>  threshed_unweighted_image(unweight_drizzle_results(full_threshed_output, full_output_weights));
            //    //            if (ds.post_drizzle_col_debias_on) image_cleanup::col_10p_debias(threshed_unweighted_image);
            //                threshed_unweighted_image.write_to_file(drizzle_output_dir+filename_percentile_text+"_PC.fits");
            //
            ////                reference_psf combined_image = combine_regular_and_threshed_bitmaps(
            ////                            unweighted_full_image, threshed_unweighted_image,
            ////                            1.0/
            ////                            0.2,
            ////                            ds.drizzle_zoom_factor*2.0);
            ////                combined_image.write_to_file(drizzle_dir+filename_percentile_text+"_combined.fits" );
            ////                combined_image.write_mask_to_file(drizzle_dir+filename_percentile_text+"_combined_thresh_mask.fits" );
            //            }
        }
    }

//    tick_count t2 = tick_count::now();
//    float total_time =(t2-t0).seconds();
//    float megabytes_loaded=input_buffer.bytes_loaded()/1024./1024.;
//    cout <<"Average disk access speed of " << megabytes_loaded /total_time <<" MB/sec"<<endl;
//    cout<< (opts.output_diagnostic_files ?
//        "Marked frames written to folder"+diagnostics_output_folder_full_path :
//        "Rerun with diagnosis output on to save frames")<<endl;
//

    return 0;
}

