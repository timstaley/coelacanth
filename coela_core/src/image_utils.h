/*
 * File:   image_utils.h
 * Author: ts337
 *
 * Created on 26 August 2008, 16:33
 */
#ifndef COELA_IMAGE_UTILS_H
#define COELA_IMAGE_UTILS_H

#include "ccd_image.h"



namespace coela {

///forward declaration:
//struct DrizzleSettings;

namespace image_utils {
//=========================================================================================================

//All Flux preserving:
template<typename input_datatype>
CcdImage<double> bicubic_resample_non_edge_region(
    const CcdImage<input_datatype>& input_bmp,
    const CcdBoxRegion,
    const double resample_factor);

///Nb for coord x1 in orig--> x2 in resampled;      x2 = 2*(x1) - 1
CcdImage<double> bicubic_resample_2x(const CcdImage<double>& input_bmp);
/// Similarly x3 = 4*(x1) - 3
CcdImage<double> bicubic_resample_4x(const CcdImage<double>& input_bmp);

CcdImage<double> bin_image(const CcdImage<double>& input_image,
                           const int bin_pixel_width);


//
//        float_bitmap interpolate_non_uniform_sampling_points(
//                const pixel_box& image_shape,                 const ImageGrid<coordinate_type::CCD>& ref_grid,
//                const vector<CcdPosition>& sample_Positions, const vector<double>& sample_values);
//
//        float_bitmap interpolate_non_uniform_sampling_points(
//                const pixel_box& image_shape,
//                const vector<PixelPosition>& sample_Positions, const vector<double>& sample_values);
//
//        float_bitmap inverse_distance_weighting_map(const pixel_box& map_shape, const PixelPosition& sample_Position);




//
//        float_bitmap median_filter(const float_bitmap& input, const int neighbour_depth=1);
//
//        float_bitmap masked_median_filter(const float_bitmap& input, const float_bitmap& mask, const int neighbour_depth=1);
//        float_bitmap masked_box_smooth(const float_bitmap& input, const float_bitmap& mask, const int box_width=3);
//=========================================================================================================
} //end namespace image_utils
}//end namespace coela

#endif  /* _IMAGE_UTILS_H */

