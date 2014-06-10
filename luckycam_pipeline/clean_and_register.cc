/*
 * File:   lucky_clean_and_register.cc
 * Author: ts337
 * stuff
 * Created on 18 August 2009, 11:42
 */

//#include <stdlib.h>

#include "coela_utility/src/lucky_math_funcs.h"

#include "coela_core/src/image_utils.h"
#include "coela_core/src/DS9Region.h"

#include "coela_luckypipe/src/CCD_dataset_info.h"
#include "coela_luckypipe/src/camera_config_info.h"

#include "coela_luckypipe/src/gain_utils.h"
#include "coela_luckypipe/src/single_CCD_filters.h"
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

/*
 * Program flow and functionality:
 * Load file list, either from index or generated.
 * Determine file format and adjust regions accordingly
 * Clean each file. Add to averages. Optionally, write decompressed to file.
 * Gather gain information for the CCD for this run. Write gain_info to CCD base dir.
 * If registering: Register all guide stars on this CCD, write results to frame_list in output base dir.
 *
 */

using namespace coela;
//using namespace coela::image_cleanup;
using namespace std;
//using psf_models::reference_psf;
using tbb::tick_count;
//using namespace  single_CCD_threading_tools;

struct xcorr_templates {
    static string airy1x_parabolic_fit;
    static string airy4x;
    static string airy4x_parabolic_fit;
    static string normed_airy1x_parabolic_fit;
    static string normed_airy4x;
    static string normed_airy4x_parabolic_fit;
    static vector<string> get_defined_template_types();

};

string xcorr_templates::airy1x_parabolic_fit="airy_core_1x_pb";
string xcorr_templates::airy4x="airy_core_4x";
string xcorr_templates::airy4x_parabolic_fit="airy_core_4x_pb";
string xcorr_templates::normed_airy1x_parabolic_fit="normalised_core_1x_pb";
string xcorr_templates::normed_airy4x="normalised_core_4x";
string xcorr_templates::normed_airy4x_parabolic_fit="normalised_core_4x_pb";

vector<string> xcorr_templates::get_defined_template_types()
{
    vector<string> v;
    v.push_back(xcorr_templates::airy1x_parabolic_fit);
    v.push_back(xcorr_templates::airy4x);
    v.push_back(xcorr_templates::airy4x_parabolic_fit);
    v.push_back(xcorr_templates::normed_airy1x_parabolic_fit);
    v.push_back(xcorr_templates::normed_airy4x);
    v.push_back(xcorr_templates::normed_airy4x_parabolic_fit);
    return v;
}


struct ProgramOptionSet {
    string dataset_info_filename;
    CCD_DatasetInfo dataset;

    CameraConfigInfo camconf;

    int n_threads;
    int n_frames_to_process;

    bool output_cleaned_files;
    bool output_diagnostic_files;
    string bias_pedestal_diagnostics_output_subfolder;

    bool output_windowed_images;
    string window_region_filename;
    string cleaned_output_folder;

    double cosmic_ray_candidate_threshold;
    string cosmic_ray_diagnostics_output_subfolder;

    bool imprint_file_time_stamp;

    bool generate_column_bias_frame;
    string precalc_col_bias_frame_path;
    size_t n_frames_to_use_for_column_bias_gen;

    bool frame_by_frame_temporal_debiasing_on;
    CCD_BoxRegion temporal_debiasing_hist_region;
    PixelRange temporal_debiasing_hist_box;

    bool build_faint_histogram;
    bool build_raw_histogram_only; //Grab histogram only, do not use for standard debiasing.
    string faint_hist_region_filename;
    CCD_BoxRegion
    faint_hist_region; //Usually describes a region of low flux suitable for EM gain determination.
    string faint_histogram_output_filename; //However, the user may choose to create a histogram of a bright region to check for saturation.

    bool dynamic_row_debiasing_on;

    bool gs_registration_on;
    string guide_star_region_file;
    vector<CCD_BoxRegion> guide_star_regions;
    double analysis_resample_factor;
    bool perform_centroiding;
//    bool using_log_kernel;
    bool use_parabola_fit;
    string xcorr_template_type;
    string custom_psf_model_path;
    double convolution_thresh;

    bool output_avg_images;

};

ProgramOptionSet get_default_options()
{
    ProgramOptionSet default_opts;
    //Set hard-coded defaults:
    default_opts.n_threads=6;
    default_opts.n_frames_to_process=0; //if zero, simply processes all frames in list

    default_opts.output_cleaned_files=false;
    default_opts.output_diagnostic_files=true;
    default_opts.bias_pedestal_diagnostics_output_subfolder = "bias_pedestal_diagnostics";

    default_opts.output_windowed_images=false;

    default_opts.imprint_file_time_stamp=false;

    default_opts.cleaned_output_folder = "clean_and_reg_output_frames";

    default_opts.cosmic_ray_diagnostics_output_subfolder = "cosmic_ray_candidates/";
    default_opts.cosmic_ray_candidate_threshold = 8000.0;

    default_opts.precalc_col_bias_frame_path="";
    default_opts.generate_column_bias_frame=true;
    default_opts.n_frames_to_use_for_column_bias_gen=3000;

    default_opts.frame_by_frame_temporal_debiasing_on = true;

    default_opts.build_faint_histogram=false;
    default_opts.build_raw_histogram_only=false;
    default_opts.faint_histogram_output_filename="bg_histogram.txt";


    default_opts.dynamic_row_debiasing_on=false;

    default_opts.gs_registration_on=false;
    default_opts.analysis_resample_factor=4.0;
    default_opts.perform_centroiding=false;

    default_opts.use_parabola_fit= false;
    default_opts.xcorr_template_type=xcorr_templates::airy4x;
    default_opts.convolution_thresh=0.25;

    default_opts.output_avg_images=true;

    return default_opts;
}

void print_usage_and_exit(const ProgramOptionSet& defaultopts)
{
    cerr<<"Usage: Cleanup CCD_DatasetInfo.txt  \n\n"
        "[-g --gs_registration [gs_rgn.reg] ]     ---Perform guide star registration. if no region file specified, use default from CCD_dataset_info.txt\n"
        "[-b --build_histogram [background_rgn.reg] ]     ---Gather histogram. If no region file specified, use default from CCD_dataset_info.txt\n"
        "[-n --nthreads (number)]           ---run the pipeline with n threads (default 6) \n\n"
        "[-m --nframes (number) ]           ---do a quick test run on the first ... frames\n\n"
//            "[-f]               ---Use frame mode zeroing (faster than hist fitting, but less accurate) \n\n"
        "[-d  true/false <n_frames> ]  ---determine bias frame from first <n_frames> frames\n"
        "(defaults):"<< string_utils::bool_to_string(defaultopts.generate_column_bias_frame)
        <<" , "
        <<defaultopts.n_frames_to_use_for_column_bias_gen <<" frames). \n\n"
        "[-c  path\\to\\col_bias_frame.fits]     ---use the specified column bias frame rather than generate a fresh one \n\n"

        "[-H           ---Just gather a histogram from the given region, \n"
        "                 don't perform debiasing or output avg files (handy for taking a histogram of a possibly saturated image region)\n\n"
        "[-h (histogram_output_filename)(.txt will be appended)]     ---Output gathered histogram to a non-standard location (default: gain_histogram.txt)\n\n"
        "[-p]                               ---Use parabolic fitting of the xcorr peak (default false)\n\n"
        "[--xcorr <template_type>]      ---Specify xcorr template type(default \"airy_core\")\n\n"
        "[-a --centroid]               ---Take a centroid of the guide star region for use with autoguider sims\n\n"
        "[-r threshval]   ---Cosmic ray threshold value: default 6000 \n"
        "[-O 0/1/2]     ---Output files for visual inspection (0=none, 1=marked_for_diagnosis_only, 2=all)(default:1)\n\n"

        "[-T  / --ftimestamp ]      --- Insert an approximate timestamp into the header of every file. \n"
        "[-w  window_Region.reg]    --Output a cleaned, windowed region from each file (Not yet implemented)"
        <<endl; exit(0);
}



ProgramOptionSet load_options(int argc, char** argv)
{
    ProgramOptionSet default_opts = get_default_options();
    if (argc<2) { print_usage_and_exit(default_opts); }

    ProgramOptionSet opts=default_opts;
    opts.dataset_info_filename = argv[1];
    opts.dataset = CCD_DatasetInfo(opts.dataset_info_filename);

    opts.camconf = CameraConfigInfo(opts.dataset.default_camera_config_file);

    using namespace GetOpt;
    GetOpt::GetOpt_pp args(argc, argv);
    args.exceptions(std::ios::failbit);

    try {
        if (args >> OptionPresent('g', "gs_registration")) {
            opts.gs_registration_on=true;
            args >> Option('g', "gs_registration", opts.guide_star_region_file);
            if (opts.guide_star_region_file.empty()) {
                opts.guide_star_region_file= opts.dataset.default_guide_star_region_file;
            }
            opts.guide_star_regions = ds9::load_vector_from_file<CCD_BoxRegion>
                                      (opts.guide_star_region_file);
            cerr<<"********    "<<opts.guide_star_regions.size()<< "         **********"<<endl;
            if (opts.guide_star_regions.empty()) {
                throw runtime_error("Error generation option set - no guide regions found in" +
                                    opts.guide_star_region_file);
            }
            cout<<"Will track "<< opts.guide_star_regions.size()<< " GS regions loaded from file:"<<
                opts.guide_star_region_file<<endl;
        }

        if (!opts.camconf.simulated_data) {
            opts.frame_by_frame_temporal_debiasing_on=
                default_opts.frame_by_frame_temporal_debiasing_on;
            if (opts.frame_by_frame_temporal_debiasing_on) {
                opts.temporal_debiasing_hist_region =
                    opts.camconf.get_calibration_info_for_CCD_id(
                        opts.dataset.ccd_id).default_temporal_debiasing_histogram_region;
            }
        } else {opts.frame_by_frame_temporal_debiasing_on=false;}

        if (args >> OptionPresent('b', "build_histogram")) {
            opts.build_faint_histogram=true;
            args >> Option('b', "build_histogram", opts.faint_hist_region_filename);
            if (opts.faint_hist_region_filename.empty()) {
                opts.faint_hist_region_filename= opts.dataset.default_faint_histogram_region_file;
            }
            opts.faint_hist_region = ds9::load_vector_from_file<CCD_BoxRegion>
                                     (opts.faint_hist_region_filename).front();
            cout<<"Will build histogram for region loaded from file:"<<
                opts.faint_hist_region_filename<<endl;
        }

        args >> Option('n', "nthreads", opts.n_threads, default_opts.n_threads);
        args >> Option('m', "max_frames", opts.n_frames_to_process,
                       default_opts.n_frames_to_process);

        int output_cleaned_files_option;
        if (args >> OptionPresent('O')) {
            args >> Option('O', "output_cleaned_files", output_cleaned_files_option);
            if (output_cleaned_files_option>1) { opts.output_cleaned_files=true; }
            if (output_cleaned_files_option<1) { opts.output_diagnostic_files=false; }
        }

        if (args >> OptionPresent('w')) {
            args >> Option('w', opts.window_region_filename);
        }



        if (args >> OptionPresent('d')) {
            vector<string> option_strings;
            args >> Option('d', option_strings);
            using namespace string_utils;

            if (option_strings.size()==2 && string_to_bool(option_strings.front())==true) {
                opts.generate_column_bias_frame=true;
                opts.n_frames_to_use_for_column_bias_gen = atoi(option_strings[1]);
            } else if (option_strings.size()==1 && string_to_bool(option_strings.front())==true) {
                opts.generate_column_bias_frame=true;
                opts.n_frames_to_use_for_column_bias_gen=default_opts.n_frames_to_use_for_column_bias_gen;
            } else {
                opts.generate_column_bias_frame=false;
            }

            if (opts.generate_column_bias_frame) {
                cout<<"Will generate column bias frame from first "<<
                    opts.n_frames_to_use_for_column_bias_gen<< "frames"<<endl;

            }
        }

        if (args >> OptionPresent('c')) {
            args >> Option('c', opts.precalc_col_bias_frame_path);
            opts.generate_column_bias_frame=false;
            cout<<"Will use user supplied col-bias frame at "<< opts.precalc_col_bias_frame_path
                <<endl;
            if (args >> OptionPresent('d')) { cout<<"WARNING: This overrides column bias frame generation options;"<<endl; }
        }


        args>> Option('r', opts.cosmic_ray_candidate_threshold,
                      default_opts.cosmic_ray_candidate_threshold);
        cout<<"Cosmic ray candidate threshold counts: "<<opts.cosmic_ray_candidate_threshold<<endl;


        if (args >> OptionPresent('p')) {
            opts.use_parabola_fit=true;
            cout<<"Will use parabola fitting of Xcorr Position."<<endl;
        }

        if (args >> OptionPresent("xcorr")) {
            args >> Option("xcorr", opts.xcorr_template_type);
            cout<<"Will use xcorr template type:"<< opts.xcorr_template_type <<endl;
            vector<string> v = xcorr_templates::get_defined_template_types();
            if (find(v.begin(), v.end(), opts.xcorr_template_type) == v.end()) {
                throw runtime_error(
                    "Invalid cross correlation template type specified, got:" +
                    opts.xcorr_template_type);
            }
        }



        if (args >> OptionPresent('a', "centroid")) {
            opts.perform_centroiding=true;
            cout<<"Will record GS centroids for Autoguider sim."<<endl;
        }


        if (args>>OptionPresent('H')) {
            opts.build_raw_histogram_only=true;
            opts.output_avg_images = false;
            opts.gs_registration_on=false;
            opts.build_faint_histogram=true;
            cout << "Will only gather histogram from " << opts.faint_hist_region_filename << endl
                 << " (No temporal debias, avgs or registration)" << endl;
        }

        args >> Option('h', "histogram_output_filename", opts.faint_histogram_output_filename,
                       default_opts.faint_histogram_output_filename);

        if (args >> OptionPresent('T', "ftimestamp")) {
            opts.imprint_file_time_stamp=true;
        }
    } catch (exception& e) {
        cerr<<"Error processing option set; please check command line args "
            "and configuration files. Error message reads:\n"<<e.what()<<endl;
        exit(1);
    }


    PixelRange primary_guide_star_box(1,1,1,1);



    return opts;
}



int main(int argc, char** argv)
{

    ProgramOptionSet opts=load_options(argc, argv);
    cout <<"\n*****Lucky Cleanup ***** " <<endl;

    const CCD_DatasetInfo& cdi(opts.dataset);
    cout<<"CCD output dir: "<<cdi.CCD_outputdir<<endl;




    if (opts.output_diagnostic_files) { boost::filesystem::create_directories(cdi.CCD_outputdir+opts.bias_pedestal_diagnostics_output_subfolder); }
    if (opts.output_cleaned_files) { boost::filesystem::create_directories(opts.cleaned_output_folder); }

//    boost::filesystem::create_directories(opts.cosmic_ray_diagnostics_output_folder);
//    boost::filesystem::create_directories(opts.cosmic_ray_diagnostics_output_folder+"/confirmed");
//    boost::filesystem::create_directories(opts.cosmic_ray_diagnostics_output_folder+"/rejected");

    const CameraConfigInfo& camconf(opts.camconf);

    //******   Load cleanup settings, get file list =================
    cout <<"Cleanup settings:\n" <<cdi<<endl;

    vector<FrameInfo> frame_vec(
        FrameInfo::get_frames_vec(cdi.CCD_inputdir,
                                  cdi.in_filestem,
                                  cdi.extension,
                                  cdi.ccd_id,
                                  cdi.dataset_output_base_dir));

    //Initialize the stable bias pedestal estimate to something sensible
    //If bias drift tracking is enable these should get overwritten before ever getting used.
    FrameInfo::set_bias_pedestal_estimates(frame_vec, cdi.uniform_bias_pedestal_estimate);

    FrameInfo::write_list_to_file(frame_vec, cdi.CCD_outputdir+"full_file_list.txt");

    //Effectively taken care of by loading filter... only need this if we want to write out a "files_to_process" list.
//    if (opts.n_frames_to_process!=0){ //{n_frames_to_process=frame_list.size(); }//process all by default//
//        if ( (size_t)out.n_frames_to_process < frame_list.size()){
//            vector<FrameInfo> temp(frame_list.begin(), frame_list.end());
//            temp = vector<FrameInfo>(temp.begin(), temp.begin()+n_frames_to_process);
//            frame_list = list<FrameInfo>(temp.begin(), temp.end());
//        }
//    }


//================================= Get initial frame for setting pixel boxes...  =================================================

    assert(!frame_vec.empty());
    HistogramContainer14bit temp_hist;
    CCDImage<float> first(frame_vec.front().file_path,frame_vec.front().header_byte_offset);
    first.initialize_CCD_grid_for_raw_data();

    FitsHeader first_header(frame_vec.front().file_path,frame_vec.front().header_byte_offset);

    first.write_to_file(cdi.CCD_outputdir+"first_input_image_raw.fits");
    PixelRange cropped_PixelRange;
    if (!opts.camconf.simulated_data) {
        cropped_PixelRange=camconf.get_calibration_info_for_CCD_id(cdi.ccd_id).cropped_PixelRange;
        first = CCDImage<float>::sub_image(first, cropped_PixelRange) ;
        first.write_to_file(cdi.CCD_outputdir+"first_input_image_cropped.fits");
    }

    if (opts.frame_by_frame_temporal_debiasing_on) {
        opts.temporal_debiasing_hist_box=
            first.CCD_grid.corresponding_pixel_region(
                opts.temporal_debiasing_hist_region).bounded_pixels();
    }

    PixelRange window_box;
    if (opts.output_windowed_images) {
        window_box =
            first.CCD_grid.corresponding_pixel_region(
                ds9::load_vector_from_file<CCD_BoxRegion>(opts.window_region_filename).front()
            ).expand_to_pixel_boundaries().bounded_pixels();
    }

    //=========================  Load any relevant background regions =========================================
    map<int,long> gain_histogram;
    CCD_BoxRegion gain_CCD_rgn;
    PixelRange gain_box;

    if (opts.build_faint_histogram) {
        CCD_BoxRegion gain_CCD_rgn = opts.faint_hist_region;
        gain_CCD_rgn.expand_to_pixel_boundaries() ;
        gain_box =  first.CCD_grid.corresponding_pixel_region(gain_CCD_rgn).bounded_pixels() ;
        cout <<"Building histogram from CCD_rgn " <<gain_CCD_rgn<<endl;
        cout<<"Will save this histogram to file: " <<opts.faint_histogram_output_filename<<endl;
    } else {
        gain_box=PixelRange(1,1,1,
                            1); //need this to be valid so we can still init. relevant filters, even if not in use.
    }

//---------------------------------------------------------------------------------------------------
//create / load column bias frame

    CCDImage<float> combined_bias_frame;
    if (!camconf.simulated_data) {
        if (opts.generate_column_bias_frame) {
            using namespace clean_and_register_subroutines;
            combined_bias_frame = estimate_column_bias_pattern(frame_vec,
                                  cdi.CCD_outputdir,
                                  cropped_PixelRange,
                                  opts.frame_by_frame_temporal_debiasing_on,
                                  opts.temporal_debiasing_hist_box,
                                  opts.n_frames_to_use_for_column_bias_gen, opts.n_threads);

            opts.dataset.uniform_bias_pedestal_estimate = combined_bias_frame.pix.min_val();
            combined_bias_frame.write_to_file(cdi.CCD_outputdir
                                              +"cleanup_generated_column_bias_frame.fits");
        } else if (!opts.precalc_col_bias_frame_path.empty()) {
            cout<<"\nLoading precal col bias frame from:  "<< opts.precalc_col_bias_frame_path
                <<"\n"<<endl;
            combined_bias_frame=CCDImage<float>(opts.precalc_col_bias_frame_path);
            combined_bias_frame.write_to_file(cdi.CCD_outputdir+"loaded_column_bias_frame.fits");
        } else {
            cout<<"\n**** WARNING: No column bias frame supplied, available or generated - column bias will not be corrected *****\n"<<endl;
            combined_bias_frame.pix = PixelArray2d<float>(first.pix.range().x_dim(),
                                      first.pix.range().y_dim(), 0.0);
        }

        combined_bias_frame.write_to_file(cdi.CCD_outputdir+"column_bias_frame_last_used.fits");
        opts.dataset.most_recently_used_column_bias_frame = cdi.CCD_outputdir
                +"column_bias_frame_last_used.fits";
        opts.dataset.write_to_file(opts.dataset_info_filename);

        if (opts.camconf.get_calibration_info_for_CCD_id(
                    cdi.ccd_id).precal_row_bias_frame_available) {
            string path = opts.camconf.get_calibration_info_for_CCD_id(
                              cdi.ccd_id).precal_row_bias_frame_path;
            cout<<"Loading row bias levels from file \""<<path<<"\""<<endl;
            CCDImage<float> precal_row_bias_frame(path);
            combined_bias_frame.pix += precal_row_bias_frame.pix;
        }

        CCDImage<float> first_debiased(first);
        cerr<<"Debiasing sample frame...";
        cout<<first_debiased.pix.range()<<" ; "<<combined_bias_frame.pix.range()<<endl;
        first_debiased.pix-=combined_bias_frame.pix;
        first_debiased.write_to_file(cdi.CCD_outputdir+"first_input_image_debiased.fits");
        cerr<<"Done"<<endl;

    } else {
        cout<<"\n**** Simulated data --- proceeding with blank bias frame *****\n"<<endl;
        combined_bias_frame.pix = PixelArray2d<float>(first.pix.range().x_dim(),
                                  first.pix.range().y_dim(), 0.0);
    }



//=========================  Load any relevant guide regions, generate our kernel =========================================

    psf_models::reference_psf kernel;

    PixelRange centroiding_box;

    if (!opts.guide_star_regions.empty()) {
        if (opts.gs_registration_on) {
            //*******Generate kernel if required
            //Deal with all template types requiring no resampling
            if (opts.xcorr_template_type == xcorr_templates::airy1x_parabolic_fit ||
                    opts.xcorr_template_type == xcorr_templates::normed_airy1x_parabolic_fit) {
                opts.analysis_resample_factor=1.0;
            } else {
                opts.analysis_resample_factor=4.0;
            }


            if (opts.xcorr_template_type == xcorr_templates::airy1x_parabolic_fit ||
                    opts.xcorr_template_type == xcorr_templates::airy4x ||
                    opts.xcorr_template_type == xcorr_templates::airy4x_parabolic_fit) {
                kernel=
                    clean_and_register_subroutines::generate_airy_core_template(
                        camconf,
                        cdi.ccd_id,
                        1.0 / opts.analysis_resample_factor
                    );
            } else if (opts.xcorr_template_type == xcorr_templates::normed_airy1x_parabolic_fit ||
                       opts.xcorr_template_type == xcorr_templates::normed_airy4x ||
                       opts.xcorr_template_type == xcorr_templates::normed_airy4x_parabolic_fit) {
                kernel = clean_and_register_subroutines::
                         generate_normalised_core_template(
                             camconf,
                             cdi.ccd_id,
                             1.0/opts.analysis_resample_factor
                         );

            } else {
                throw runtime_error("\n*** Error in xcorr template creation *** \n");
            }

            //Deal with all template types requiring no parabolic fitting
            if (opts.xcorr_template_type == xcorr_templates::airy4x ||
                    opts.xcorr_template_type == xcorr_templates::normed_airy4x) {
                opts.use_parabola_fit = false;
            } else {
                opts.use_parabola_fit = true;
            }

            kernel.psf_image.write_to_file(cdi.CCD_outputdir+"kernel.fits");
            kernel.mask.write_to_file(cdi.CCD_outputdir+"kernel_mask.fits");

            cout <<"Registering guide stars at ";
            for (size_t i = 0; i != opts.guide_star_regions.size(); ++i) {
                opts.guide_star_regions[i].shrink_to_pixel_boundaries();
                cout << "\t[" << opts.guide_star_regions[i] << "]\t";
            }
        }//end if (cs.gs_registration_on)
        if (opts.perform_centroiding) {
            centroiding_box = (first.CCD_grid.corresponding_pixel_region(
                                   opts.guide_star_regions.front())
                              ).shrink_to_pixel_boundaries().bounded_pixels();
        }
    }//end if (!gs_regions.empty())



//=========================  Set up machinery using parameters from first file in list =========================================


//    hist_min = first.get_min_value();
//    cout <<"First frame has min value of " << hist_min << " once cleaned ";
    int  hist_min = -250;
    cout <<"- initialising histogram with min value " <<hist_min<<endl;
    //run_histogram.set_min_value(hist_min);

    if (opts.imprint_file_time_stamp) {
        cout<<"Imprinting time stamps using approximate starting point of first file timestamp and millisecond exposure duration of: "
            <<camconf.get_calibration_info_for_CCD_id(
                frame_vec.front().ccd_id).milliseconds_per_frame_exposure<<endl;

    }

//=========================  Begin pipeline construction ========================================================

//==========================    Create Filters:   ====================================
    cout <<"Initialising pipeline with "<<opts.n_threads<<" threads..."<<endl;
    tbb::task_scheduler_init init_user(opts.n_threads);
    tbb::pipeline pipeline;

    using namespace single_CCD_filters;

    Sequential_File_Buffer_Filter<> input_filter(frame_vec, opts.n_frames_to_process);
    Decompress_Filter decompressor;

    pipeline.add_filter(input_filter);
    pipeline.add_filter(decompressor);

    Timestamp_Load_Filter timestamp_loader; //mainly for use with simulated data
    if (cdi.extension.find("fits")!=string::npos && first_header.key_exists("FRAMETIM")) {
        pipeline.add_filter(timestamp_loader);
    }

    File_Timestamp_Imprinting_Filter f_timestamper(frame_vec.front(),
            camconf.get_calibration_info_for_CCD_id(
                frame_vec.front().ccd_id).milliseconds_per_frame_exposure);
    if (opts.imprint_file_time_stamp) { pipeline.add_filter(f_timestamper); }

    Frame_Count_Display_Filter count_filter(true);
    pipeline.add_filter(count_filter);

    FrameCropFilter crop_filter(cropped_PixelRange);

    BiasDriftTracker t_db_filter(opts.temporal_debiasing_hist_box);

    UniformDebias uniform_debias_filter;
    Frame_Summation zeroed_avg_filter(first.pix.range(), first.CCD_grid.image_outline_);

    BiasFrameSubtractor static_col_db_filter;
    static_col_db_filter.set_bias_frame(combined_bias_frame);

    Frame_Summation debiased_sum_filter(first.pix.range(), first.CCD_grid.image_outline_);


    if (camconf.simulated_data==false) {
        pipeline.add_filter(
            crop_filter); //Usually a null op for simulated data, but it ensures correct initialization, etc.
        if (!opts.build_raw_histogram_only) {
            if (opts.frame_by_frame_temporal_debiasing_on) { pipeline.add_filter(t_db_filter); }
            pipeline.add_filter(uniform_debias_filter);
            pipeline.add_filter(zeroed_avg_filter);
            pipeline.add_filter(static_col_db_filter);
        }
    }

    pipeline.add_filter(debiased_sum_filter);


    Dynamic_Row_Debias row_db_filt;
    if (opts.dynamic_row_debiasing_on) { pipeline.add_filter(row_db_filt); }

    Frame_Output_Filter cleaned_output_filter(opts.cleaned_output_folder);
    if (opts.output_cleaned_files) { pipeline.add_filter(cleaned_output_filter); }

    Histogram14BitBuildFilter cleaned_hist_filter(gain_box, hist_min);
    if (opts.build_faint_histogram) { pipeline.add_filter(cleaned_hist_filter); }

    //To do - add a "STOP" limit, where if too many frames have candidates, the filter throws an exception and brings the program to a halt.
    // Otherwise, a threshold too low causes excessive diagnosis output and slowdown.
    // The limit should be implemented as a percentage of the frames processed so far, once say the first 1000 have been done.
    Cosmic_Ray_Candidate_Detection cosmic_ray_detector(opts.cosmic_ray_candidate_threshold);

    //NB All further cosmic ray detection depends on candidates, so turning this off effectively turns it all
    if (!camconf.simulated_data) {
        pipeline.add_filter(cosmic_ray_detector);
    }

    Cosmic_Ray_Confirmation cosmic_ray_checker(10.0, 6.0);
    pipeline.add_filter(cosmic_ray_checker);


    CrossCorrelator reg_filter(
        opts.guide_star_regions, kernel, opts.analysis_resample_factor, opts.convolution_thresh,
        kernel.psf_image.pix.range().x_dim()/1.5,
        opts.use_parabola_fit);
    if (opts.gs_registration_on) { pipeline.add_filter(reg_filter); }

    Region_Centroid centroid_filter(centroiding_box, frame_vec.size());
    if (opts.perform_centroiding) { pipeline.add_filter(centroid_filter); }



    Cosmic_Ray_Diagnosis_Output ray_output_filter(cdi.CCD_outputdir,
            opts.cosmic_ray_diagnostics_output_subfolder);
    pipeline.add_filter(ray_output_filter);

    Frame_Info_Collection_Filter list_collector(frame_vec.size());
    pipeline.add_filter(list_collector);


    Serial_Decommission_Filter end_filter;
    pipeline.add_filter(end_filter);





//==========================    Pipeline constructed; run it:   ====================================
    tick_count t0;
    cout <<"Pipelining "<<input_filter.n_frames_set_to_process()<<" of "<< frame_vec.size()
         <<" files..."<<endl;
    t0 = tick_count::now();
    pipeline.run(lucky_n_tokens);
    tick_count t1 = tick_count::now();
    double time =(t1-t0).seconds();
    cout <<"Averaged "<<count_filter.n_counted()<< " files in "<<
         time<<" seconds at an avg of "<< time/ count_filter.n_counted()<<" secs per frame or "<<
         count_filter.n_counted()/ time <<" fps"<<endl;
    cout <<"Loaded: "
         <<input_filter.bytes_loaded()/1e6<<"MiB at Averaged file access rate: "<<
         input_filter.bytes_loaded() / 1e6 / time <<" MiB per sec"<<endl;


//==========================    Process gain information if available ====================================
    gain_utils::gain_info hist_fit;
    if (opts.build_faint_histogram) {
        cout<<"Fitting gain hist..."<<endl;

//        bmp_histogram_container raw_hist=raw_hist_filter.histogram();

        HistogramContainer14bit clean_hist=cleaned_hist_filter.histogram();
        string histogram_output_full_path = cdi.CCD_outputdir+
                                            opts.faint_histogram_output_filename;
        clean_hist.write_to_file(histogram_output_full_path);

        try {
            hist_fit = gain_utils::fit_CCD_histogram(clean_hist.convert_to_value_count_map());
            hist_fit.write_to_file(string_utils::strip_file_extension(histogram_output_full_path)
                                   +"_fit_pars.txt");

            opts.dataset.most_recently_output_histogram = histogram_output_full_path;
            opts.dataset.write_to_file(opts.dataset_info_filename);

        } catch (runtime_error& e) {
            cerr<<string("Encountered an error while attempting to deal with gain info:\n")+e.what()
                <<endl;
            cerr<<"Will continue regardless..."<<endl;
            opts.build_faint_histogram=false;
        }
    }


//==========================    Gather average images at various levels of processing, write to files   ====================================
    if (opts.output_avg_images) {
        cout<<"Outputting average images..."<<flush;
        cout<<"dir: "<<cdi.dataset_output_base_dir<<endl;
        cout<<"To folder:\n"<<cdi.dataset_output_base_dir  +"/avgs/" <<endl;


        string avgs_filestem = cdi.dataset_output_base_dir  +"/avgs/" +cdi.out_filestem+"_"+
                               string_utils::itoa(cdi.ccd_id)+"_";
        int n_files=count_filter.n_counted();
        if (opts.n_frames_to_process!=0) {
            //(If we were just doing a quick snapshot)
            avgs_filestem = cdi.dataset_output_base_dir  +"/avgs/" +"First_"+ string_utils::itoa(
                                n_files)+"_"+cdi.out_filestem+"_"+ string_utils::itoa(cdi.ccd_id)+"_";
        }
        boost::filesystem::create_directories(cdi.dataset_output_base_dir  +"/avgs/");

        FitsHeader additional_hdr_info;
        additional_hdr_info.add_keyword("CCD_ID", string_utils::itoa(cdi.ccd_id),
                                        "CCD ID number");


        if (!camconf.simulated_data) {
            CCDImage<double> mask = CCDImage<double>(image_cleanup::get_CCD_default_weight_map(
                                        camconf.get_calibration_info_for_CCD_id(cdi.ccd_id)));
            mask.write_to_file(cdi.dataset_output_base_dir  +"/avgs/" +cdi.out_filestem+"_"+
                               string_utils::itoa(cdi.ccd_id)+"mask.fits");
            CCDImage<double> dc_zeroed = zeroed_avg_filter.Avg();
            //        dc_zeroed.add_comment("Run filestem:  "+ cs.in_filestem);
            //        dc_zeroed.set_keyword("CAMERAID", string_utils::itoa(cs.ccd_id));
            dc_zeroed.pix *= mask.pix;

            dc_zeroed.write_to_file(avgs_filestem   + "zeroed_avg.fits", additional_hdr_info);

            CCDImage<double> debiased_avg;
            if (camconf.simulated_data==false) { debiased_avg = debiased_sum_filter.Avg(); }
            else { debiased_avg = dc_zeroed; }

            debiased_avg.pix *= mask.pix;
            //        debiased_avg.add_comment("Run filestem:  "+ cs.in_filestem);
            //        debiased_avg.set_keyword("CAMERAID", string_utils::itoa(cs.ccd_id));

            string debias_method_str="_DB";
            debiased_avg.write_to_file(avgs_filestem  +debias_method_str+".fits",
                                       additional_hdr_info);


            CCDImage<double> normed_avg=debiased_avg;
            if (opts.build_faint_histogram) {
                gain_utils::normalise_CCD_with_uniform_gain(normed_avg, hist_fit);
                normed_avg.write_to_file(avgs_filestem + debias_method_str+"_normed.fits");

                const CCD_calibration_info& ccd_props= camconf.get_calibration_info_for_CCD_id(
                        cdi.ccd_id);
                if (ccd_props.dark_current_frames_available) {
                    cout<<"Subtracting dark current calibration frame: \""<<ccd_props.normalised_DC_frame_path<<"\""<<endl;
                    CCDImage<double> dark_current_frame =
                        CCDImage<double>::load_from_unknown_filetype(ccd_props.normalised_DC_frame_path);
                    //do dark current subtraction
                    normed_avg.pix -= dark_current_frame.pix;
                    normed_avg.write_to_file(avgs_filestem +debias_method_str+ "_normed_DC_subtracted.fits",
                                             additional_hdr_info);
                }
            }
        } else if (camconf.simulated_data) {
            CCDImage<double> sim_avg = debiased_sum_filter.Avg();
            sim_avg.write_to_file(avgs_filestem+"sim_avg.fits", additional_hdr_info);
        }

        cout<< " Done."<<endl;
    }


    //==========================  //Output frame list ====================================
    {
        string frame_list_filename= cdi.dataset_output_base_dir +"CCD_"+string_utils::itoa(
                                        cdi.ccd_id);
        if (opts.gs_registration_on) {
            //Construct a long and informative filename for the frame list...
            using string_utils::itoa;
            using string_utils::ftoa;
//            if (opts.using_log_kernel) {
//                frame_list_filename += "_log_psf_lock_";
//            } else frame_list_filename += "_reg_psf_lock_";
//            frame_list_filename+=ftoa(opts.analysis_resample_factor) + "x";
//            if (opts.use_parabola_fit)    frame_list_filename+="_parabola_fit";

            frame_list_filename+="_"+opts.xcorr_template_type;
        } else { frame_list_filename +="_noGS"; }

        if (opts.n_frames_to_process!=0) {frame_list_filename+="_first_"+string_utils::itoa(opts.n_frames_to_process);}

        frame_list_filename+=".txt";


        vector<FrameInfo> processed_frame_list = list_collector.collected_frames();
        FrameInfo::write_list_to_file(processed_frame_list, frame_list_filename);
        cout<<"Wrote frame list to "<<frame_list_filename<<endl;

        opts.dataset.uniform_bias_pedestal_estimate =
            FrameInfo::get_median_bias_pedestal_estimate(processed_frame_list);
        opts.dataset.most_recently_output_frame_list = frame_list_filename;

        if (opts.perform_centroiding) {
            using string_utils::itoa;
            string centroid_list_filename = cdi.dataset_output_base_dir +"CCD_" +itoa(
                                                cdi.ccd_id)+"_centroids.txt";
            FrameInfo::write_list_to_file(centroid_filter.collected_frames(), centroid_list_filename);
            opts.dataset.most_recently_output_frame_list = centroid_list_filename;
        }
        opts.dataset.write_to_file(opts.dataset_info_filename);
    }

    cout <<"Cleanup terminated normally, saved info for "<<list_collector.collected_frames().size()
         <<" frames"<<endl;
//          << bias_pedestal_flag_count_filter.n_counted()<<" frames marked for diagnosis"<<endl;

    if (cosmic_ray_detector.number_of_frames_with_pixel_above_threshold()) {
        cout<<cosmic_ray_detector.number_of_frames_with_pixel_above_threshold()
            <<" frames with pixels above ray analysis threshold, of which\n";
        cout<<cosmic_ray_checker.number_of_frames_with_confirmed_rays()
            <<" marked as cosmic rays "<<endl;
        FrameInfo::write_list_to_file(list_collector.collected_frames(),
                                      cdi.CCD_outputdir+"/frame_list_after_CR_scanning.txt");
    }

//    if (bias_pedestal_flag_count_filter.n_counted()){
//        cout<< ( (opts.output_diagnostic_files ) ?
//            "Dark level Marked frames written to folder "+opts.bias_pedestal_diagnostics_output_folder :
//            "Rerun with diagnosis output on to save frames")<<endl;
//    }

    if (ray_output_filter.number_of_frames_output()) {
        ofstream ray_vals(string(opts.cosmic_ray_diagnostics_output_subfolder
                                 +"/pix_vals_log.txt").c_str()) ;
        vector<double> pixel_vals = ray_output_filter.get_candidate_ray_pixel_peak_values();
        for (size_t i = 0; i!=pixel_vals.size(); ++i) {
            ray_vals<<pixel_vals[i]<<"\n";
        }
        ray_vals.close();
    }


    return (EXIT_SUCCESS);
}

