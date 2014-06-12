#include "../clean_and_register_subroutines.h"
#include <utility>
#include <iostream>
#include <coela_utility/src/lucky_math_funcs.h>
#include <coela_core/src/convolution.h>
#include "coela_core/src/pixel_array_routines.h"
using std::endl;
using std::flush;
using std::cout;
using tbb::tick_count;

namespace coela {
namespace clean_and_register_subroutines {


CCDImage<float> estimate_column_bias_pattern(
    const vector<FrameInfo>& frms,
    const string& output_dir,
    const PixelRange& crop_box,
    bool perform_temporal_debiasing,
    const PixelRange& temporal_debiasing_box,
    const size_t n_frames,
    const size_t n_threads)
{

    using namespace single_CCD_filters;

    CCDImage<float> bias_frm;
    bias_frm.pix = PixelArray2d<float>(crop_box.x_dim(), crop_box.y_dim(), 0.0);
    Sequential_File_Buffer_Filter<> input_filter(frms, n_frames);
    Decompress_Filter decompressor;
    Frame_Count_Display_Filter count_filter(true);
    FrameCropFilter crop_filter(crop_box); //Use these two if we have a good background region


    if (perform_temporal_debiasing) { assert(temporal_debiasing_box.is_valid()); }
    BiasDriftTracker drift_tracker(temporal_debiasing_box);

    UniformDebias db_filter;
    Col_Histogram_Gather col_hist_filter(crop_box, 50, crop_box.y_dim());
    col_hist_filter.set_hist_min_value(-300);

    Serial_Decommission_Filter end_filter;

    cout<<"Initializing scheduler with " << n_threads <<" threads."<<endl;

    tbb::task_scheduler_init scheduler(n_threads);
    tbb::pipeline pipeline;

    pipeline.add_filter(input_filter);
    pipeline.add_filter(decompressor);
    pipeline.add_filter(count_filter);
    pipeline.add_filter(crop_filter);
    if (perform_temporal_debiasing) { pipeline.add_filter(drift_tracker); }
    pipeline.add_filter(db_filter);
    pipeline.add_filter(col_hist_filter);
    pipeline.add_filter(end_filter);

    cout<<"Gathering histograms for column bias frame from first "<< n_frames
        <<" frames..."<<endl;


    tick_count t0;
    t0 = tick_count::now();

    pipeline.run(coela::single_CCD_filters::lucky_n_tokens);

    tick_count t1 = tick_count::now();
    double time =(t1-t0).seconds();


    string col_hists_folder=output_dir+"/col_hists/";
    boost::filesystem::create_directories(col_hists_folder);

    const vector<HistogramContainer14bit>& col_hists = col_hist_filter.get_hists();
    cout<<"\n";
    cout<<"Collected histograms of all columns in " <<n_frames<<
        " in "<<time<<" seconds or "<< time/ n_frames <<" secs per frame"<<endl;
    for (size_t x=0; x!=col_hists.size(); ++x) {
        cout<<"\r Fitting col " <<x+1 <<"\t\t"<<flush;
        col_hists[x].write_to_file(col_hists_folder+"/col"+string_utils::itoa(x)+".txt");
        double dark_estimate =
            image_cleanup::determine_histogram_bias_pedestal_via_thresholded_centroid(col_hists[x],
                    0.8);
//        double dark_estimate = gain_utils::fit_readout_gaussian(col_hists[x].convert_to_value_count_map(), false).front();
        pixel_array_routines::set_col_to_value(bias_frm.pix, x+1, (float)dark_estimate);
    }
    cout<<"\n";
    return bias_frm;
}


CCDImage<float> estimate_row_bias_pattern(
    const vector<FrameInfo>& frms,
    const std::string& output_dir,
    const PixelRange& crop_box,
    const CCDImage<float>& col_bias_pattern,
    const size_t n_frames,
    const size_t n_threads)
{
    using namespace single_CCD_filters;

    CCDImage<float> row_bias_frm;
    row_bias_frm.pix = PixelArray2d<float>(crop_box.x_dim(), crop_box.y_dim(), 0.0);
    Sequential_File_Buffer_Filter<> input_filter(frms, n_frames);
    Decompress_Filter decompressor;
    Frame_Count_Display_Filter count_filter(true);
    FrameCropFilter crop_filter(crop_box); //Use these two if we have a good background region

    UniformDebias uniform_db_filter;
    single_CCD_filters::BiasFrameSubtractor col_db_filt;
    col_db_filt.set_bias_frame(col_bias_pattern);
    single_CCD_filters::RowHistogram_Gather row_hist_filter(crop_box, 100);
    row_hist_filter.set_hist_min_value(-300);

    Serial_Decommission_Filter end_filter;

    cout<<"Initializing scheduler with " << n_threads <<" threads."<<endl;

    tbb::task_scheduler_init scheduler(n_threads);
    tbb::pipeline pipeline;

    pipeline.add_filter(input_filter);
    pipeline.add_filter(decompressor);
    pipeline.add_filter(count_filter);
    pipeline.add_filter(crop_filter);
    pipeline.add_filter(uniform_db_filter);
    pipeline.add_filter(col_db_filt);
    pipeline.add_filter(row_hist_filter);
    pipeline.add_filter(end_filter);

    cout<<"Gathering histograms for row bias frame from first "<< n_frames
        <<" frames..."<<endl;


    tick_count t0;
    t0 = tick_count::now();

    pipeline.run(coela::single_CCD_filters::lucky_n_tokens);

    tick_count t1 = tick_count::now();
    double time =(t1-t0).seconds();


    string row_hists_folder=output_dir+"/row_hists/";
    boost::filesystem::create_directories(row_hists_folder);

    const vector<HistogramContainer14bit>& row_hists = row_hist_filter.get_hists();
    cout<<"\n";
    cout<<"Collected histograms of all rows in " <<n_frames<< " in "<<time<<" seconds or "<<
        time/ n_frames <<" secs per frame"<<endl;
    for (size_t y=0; y!=row_hists.size(); ++y) {
        cout<<"\r Fitting row " <<y+1 <<"\t\t"<<flush;
        row_hists[y].write_to_file(row_hists_folder+"/row"+string_utils::itoa(y)+".txt");
        double dark_estimate =
            image_cleanup::determine_histogram_bias_pedestal_via_thresholded_centroid(row_hists[y],
                    0.8);
//        double dark_estimate = gain_utils::fit_readout_gaussian(col_hists[x].convert_to_value_count_map(), false).front();
        for (size_t x=1; x<=row_bias_frm.pix.range().x_dim(); ++x) {
            row_bias_frm.pix(x, y+1) = dark_estimate;
        }
    }
    cout<<"\n";
    return row_bias_frm;
}


CCDImage<double> create_debiased_average(
    const vector<FrameInfo>& frms,
    const PixelRange& crop_box,
    const CCD_BoxRegion crop_region,
    const CCDImage<float>& bias_frame,
    const size_t n_frames,
    const size_t n_threads)
{

    using namespace single_CCD_filters;
    Sequential_File_Buffer_Filter<> input_filter(frms, n_frames);
    Decompress_Filter decompressor;
    Frame_Count_Display_Filter count_filter(true);
    FrameCropFilter crop_filter(crop_box);

    UniformDebias uniform_db_filter;
    single_CCD_filters::BiasFrameSubtractor db_filt;
    db_filt.set_bias_frame(bias_frame);
    single_CCD_filters::Frame_Summation sum_filter(crop_box, crop_region);

    Serial_Decommission_Filter end_filter;

    cout<<"Initializing scheduler with " << n_threads <<" threads."<<endl;

    tbb::task_scheduler_init scheduler(n_threads);
    tbb::pipeline pipeline;

    pipeline.add_filter(input_filter);
    pipeline.add_filter(decompressor);
    pipeline.add_filter(count_filter);
    pipeline.add_filter(crop_filter);
    pipeline.add_filter(uniform_db_filter);
    pipeline.add_filter(db_filt);
    pipeline.add_filter(sum_filter);
    pipeline.add_filter(end_filter);

    cout<<"Creating debiased average from"<< n_frames <<" frames..."<<endl;

    tick_count t0;
    t0 = tick_count::now();

    pipeline.run(coela::single_CCD_filters::lucky_n_tokens);

    tick_count t1 = tick_count::now();
    double time =(t1-t0).seconds();

    cout<<"\n";
    cout<<"Averaged" <<n_frames<< " in "<<time<<" seconds or "<< time/ n_frames
        <<" secs per frame"<<endl;

    return sum_filter.Avg();
}


vector< std::pair<int, double> > estimate_per_column_EM_gain(
    const vector<FrameInfo>& frms,
    const std::string& output_dir,
//    const bool use_full_fit,
    const PixelRange& crop_box,
//    const CCD_BoxRegion crop_region,
    const CCDImage<float>& bias_frame,
    const size_t n_frames,
    const size_t n_threads)
{

    using namespace single_CCD_filters;
    Sequential_File_Buffer_Filter<> input_filter(frms, n_frames);
    Decompress_Filter decompressor;
    Frame_Count_Display_Filter count_filter(true);
    FrameCropFilter crop_filter(crop_box);

    UniformDebias uniform_db_filter;
    single_CCD_filters::BiasFrameSubtractor db_filt;
    db_filt.set_bias_frame(bias_frame);
    Col_Histogram_Gather col_hist_filter(crop_box, 50, crop_box.y_dim());
    col_hist_filter.set_hist_min_value(-300);
//    single_CCD_filters::Frame_Summation sum_filter(crop_box, crop_region);

    Serial_Decommission_Filter end_filter;

    cout<<"Initializing scheduler with " << n_threads <<" threads."<<endl;

    tbb::task_scheduler_init scheduler(n_threads);
    tbb::pipeline pipeline;

    pipeline.add_filter(input_filter);
    pipeline.add_filter(decompressor);
    pipeline.add_filter(count_filter);
    pipeline.add_filter(crop_filter);
    pipeline.add_filter(uniform_db_filter);
    pipeline.add_filter(db_filt);
    pipeline.add_filter(col_hist_filter);
    pipeline.add_filter(end_filter);

    tick_count t0;
    t0 = tick_count::now();

    pipeline.run(coela::single_CCD_filters::lucky_n_tokens);

    tick_count t1 = tick_count::now();
    double time =(t1-t0).seconds();

    string col_hists_folder=output_dir+"/col_hists/";
    boost::filesystem::create_directories(col_hists_folder);

    const vector<HistogramContainer14bit>& col_hists = col_hist_filter.get_hists();
    cout<<"\n";
    cout<<"Collected histograms of all columns in " <<n_frames<<
        " in "<<time<<" seconds or "<< time/ n_frames <<" secs per frame"<<endl;
    std::vector< std::pair<int, double> > col_gain_estimates;
    for (size_t x=0; x!=col_hists.size(); ++x) {
        cout<<"\r Fitting col " <<x+1 <<"\t\t"<<flush;
        col_hists[x].write_to_file(col_hists_folder+"/col"+string_utils::itoa(x)+".txt");
//        double dark_estimate = image_cleanup::determine_histogram_bias_pedestal_via_thresholded_centroid(col_hists[x], 0.8);
        gain_utils::gain_info col_fit = gain_utils::fit_CCD_histogram(
                                            col_hists[x].convert_to_value_count_map(),
                                            false
                                        );
        col_gain_estimates.push_back(std::pair<int,double>(int(x+1), col_fit.gain));

    }
    cout<<"\n";
    return col_gain_estimates;
}

//===================================================================================

psf_models::reference_psf generate_airy_core_template(
    const CameraConfigInfo& camconf,
    const int ccd_id,
    const double desired_pixel_scale_relative_to_CCD,
    const double radius_in_CCD_pix
)
{

    double rads_per_CCD_pixel =
        camconf.lens_inf.nominal_pixel_scale_in_mas*lucky_math::radians_per_milliarcsecond;
    psf_models::airy_psf_model airy_gen(
        100, rads_per_CCD_pixel,
        camconf.get_filter_info_for_ccd_id(ccd_id).central_wavelength_in_metres,
        camconf.lens_inf.telescope_outer_diameter,
        camconf.lens_inf.telescope_inner_diameter
    );

    double kernel_radius;
    if (radius_in_CCD_pix==0) {
        double minima_radius = psf_models::find_first_minima_in_CCD_pix(airy_gen);
        kernel_radius=minima_radius;
    } else {
        kernel_radius = radius_in_CCD_pix;
    }
//        2.5*minima_radius;

    psf_models::reference_psf airy_kernel =
        psf_models::generate_centred_psf(
            airy_gen,
            kernel_radius,
            desired_pixel_scale_relative_to_CCD,
            8
        );
    airy_kernel.apply_mask();

    airy_kernel.psf_image.pix /= airy_kernel.psf_image.pix.max_val();
//    psf_models::reference_psf unit_airy(airy_kernel);
    //peak is now 1.0
    double val_conv_with_unit_airy_core = spatial_convolution::detail::
                                          unsafe_convolve_with_kernel_at_Position(
                                                  airy_kernel.psf_image.pix, airy_kernel.central_pixel,
                                                  airy_kernel.psf_image.pix, airy_kernel.central_pixel
                                          );

    airy_kernel.psf_image.pix /= val_conv_with_unit_airy_core;

//    double new_conv = spatial_convolution::detail::
//        unsafe_convolve_with_kernel_at_Position(
//            unit_airy.psf_image.pix, unit_airy.central_pixel,
//            airy_kernel.psf_image.pix, airy_kernel.central_pixel
//            );
//
//    std::cerr<<"unit airy: "<<new_conv<<endl;
    return airy_kernel;
}

psf_models::reference_psf generate_normalised_core_template(
    const CameraConfigInfo& camconf,
    const int ccd_id,
    const double desired_pixel_scale_relative_to_CCD
)
{

    psf_models::reference_psf ref =
        generate_airy_core_template(camconf,
                                    ccd_id,
                                    desired_pixel_scale_relative_to_CCD);

    psf_models::reference_psf unit_airy(ref);
    unit_airy.psf_image.pix /= unit_airy.psf_image.pix.max_val();

    for (PixelIterator p(ref.psf_image.pix.range()); p!=p.end; ++p) {
        if (ref.psf_image.pix(p) > 0) {
            ref.psf_image.pix(p) /= sqrt(2*ref.psf_image.pix(p));
        }
    }

    ref.apply_mask();
    ref.psf_image.pix -= ref.psf_image.pix.min_val();

    double val_conv_with_unit_airy_core = spatial_convolution::detail::
                                          unsafe_convolve_with_kernel_at_Position(
                                                  unit_airy.psf_image.pix, unit_airy.central_pixel,
                                                  ref.psf_image.pix, ref.central_pixel
                                          );

    ref.psf_image.pix /= val_conv_with_unit_airy_core;

//    ref.psf_image.pix /= ref.psf_image.pix.max_val();
//    ref.psf_image.pix /= ref.psf_image.pix.sum();
    return ref;

}

//gaussian
//    double airy_fwhm_radians =
////        1.02* camconf.get_filter_info_for_ccd_id(ccd_id).central_wavelength_in_metres /
////            camconf.lens_inf.telescope_outer_diameter;
////
////    double airy_fwhm_CCD_pix = airy_fwhm_radians / rads_per_CCD_pixel;
////    double airy_gauss_approx_sigma = airy_fwhm_CCD_pix /
////        lucky_math::gaussian_fwhm_to_sigma_ratio;
////    double pixel_approx_sigma = 1.0 /
////        lucky_math::gaussian_fwhm_to_sigma_ratio;
////
////    double convolved_sigma = sqrt(pixel_approx_sigma*pixel_approx_sigma +
////        airy_gauss_approx_sigma + airy_gauss_approx_sigma);
////
////    psf_models::gaussian_psf_model g_model(0.75,
////            convolved_sigma
////            );


}//namespace
}//namespace
