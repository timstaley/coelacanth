/*
 * File:   clean_and_register_subroutines.h
 * Author: ts337
 *
 * Created on 12 July 2011, 13:15
 */

#ifndef COELA_CLEAN_AND_REGISTER_SUBROUTINES_H
#define COELA_CLEAN_AND_REGISTER_SUBROUTINES_H
#include "coela_utility/src/misc_math.h"

#include "coela_core/src/image_utils.h"
#include "coela_core/src/ds9_region.h"

#include "coela_luckypipe/src/ccd_dataset_info.h"
#include "coela_luckypipe/src/camera_config_info.h"

#include "coela_luckypipe/src/gain_utils.h"
#include "coela_luckypipe/src/single_ccd_filters.h"
#include "coela_luckypipe/src/image_cleanup.h"
#include "coela_luckypipe/src/frame_info.h"
#include "coela_analysis/src/psf_generation.h"

//#include "../threading/threading_tools.h"
#include <tbb/tick_count.h>
#include <tbb/task_scheduler_init.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <stdexcept>


namespace coela {
namespace clean_and_register_subroutines {

CcdImage<float> estimate_column_bias_pattern(
    const vector<FrameInfo>& frms,
    const std::string& output_dir,
    const PixelRange& crop_box,
    bool perform_temporal_debiasing,
    const PixelRange& temporal_debiasing_box,
    const size_t n_frames,
    const size_t n_threads);

CcdImage<float> estimate_row_bias_pattern(
    const vector<FrameInfo>& frms,
    const std::string& output_dir,
    const PixelRange& crop_box,
    const CcdImage<float>& col_bias_pattern,
    const size_t n_frames,
    const size_t n_threads);

CcdImage<double> create_debiased_average(
    const vector<FrameInfo>& frms,
    const PixelRange& crop_box,
    const CcdBoxRegion crop_region,
    const CcdImage<float>& bias_frame,
    const size_t n_frames,
    const size_t n_threads);

std::vector< std::pair<int, double> > estimate_per_column_EM_gain(
    const std::vector<FrameInfo>& frms,
    const std::string& output_dir,
//    const bool use_full_fit,
    const PixelRange& crop_box,
//    const CcdBoxRegion crop_region,
    const CcdImage<float>& bias_frame,
    const size_t n_frames,
    const size_t n_threads);

//===================================================================================
psf_models::ReferencePsf generate_airy_core_template(
    const CameraConfigInfo& camconf,
    const int ccd_id,
    const double desired_pixel_scale_relative_to_CCD,
    const double radius_in_CCD_pix=0 //defaults to first airy minima
);


psf_models::ReferencePsf generate_normalised_core_template(
    const CameraConfigInfo& camconf,
    const int ccd_id,
    const double desired_pixel_scale_relative_to_CCD
//    const double seeing_FWHM_C
);


}//namespace
}//namespace

#endif  /* CLEAN_AND_REGISTER_SUBROUTINES_H */

