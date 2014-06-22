/*
 * File:   image_cleanup.h
 * Author: ts337
 *
 * Created on 09 April 2009, 13:55
 */

#ifndef COELA_IMAGE_CLEANUP_H
#define COELA_IMAGE_CLEANUP_H
#include "coela_core/src/ccd_image.h"
#include "coela_utility/src/histogram_container.h"
#include "camera_config_info.h"
namespace coela {
//    typedef image<float> float_bitmap;
namespace image_cleanup {
//=====================================================================================================================

//        const CcdBoxRegion NOT09_flat_background(360,50,1050,1040);
////        const pixel_box NOT09_crop_region(49,17,49,17 ); //xmax, ymax will be reset to input_bmp.x_dim(), y_dim
//
//        const float cosmic_ray_avg_multiplier =5.0;
//        const float cosmic_ray_additive_events =100.f;
//        const int cosmic_ray_masking_pixel_radius = 3;



void col_10p_debias_preserving_bg_level(PixelArray2d<double>*,
                                        size_t y_min);
template<typename T>
void row_10p_debias_preserving_bg_level(PixelArray2d<T>* img_ptr);

template<typename T>
void row_median_debias(PixelArray2d<T>* img_ptr);

//NB clears and overwrites the "working_memory" histogram
template<class HistogramType>
double determine_bias_pedestal_from_box_in_raw_image(const PixelArray2d<float>& ,
        const PixelRange&,
        HistogramType& working_memory
                                                    );

template<class HistogramType>
double determine_histogram_bias_pedestal_via_thresholded_centroid(
    const HistogramType& hist,
    const double threshold_factor=0.8) ;

CcdImage<float> get_CCD_default_weight_map(const CcdCalibrationInfo& ccd_inf);

//=====================================================================================================================
}//end namespace coela::image_cleanup
}//end namespace coela
#endif  /* _IMAGE_CLEANUP_H */

