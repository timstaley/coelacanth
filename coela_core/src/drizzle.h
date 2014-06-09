/*
 * File:   drizzle.h
 * Author: ts337
 *
 * Created on 09 April 2009, 12:08
 */

#ifndef COELA_DRIZZLE_H
#define COELA_DRIZZLE_H

#include "pixel_array2d.h"
#include "cartesian_coords.h"

//NB implemented for PixelArray2d rather than image, to make the code as re-usable as possible.
//(Easy to implement a wrapper that checks the input / output relative scales are correct for an image with CCD_grid.)

namespace coela {
namespace drizzle {
//=====================================================================================================================================================================

/// Drizzle algorithm, as described by Fruchter and Hook (1998).
/// This is a simple "translate" version which simply applies the same translation to all input pixels (no distortion correction).
/// drizzle_scale_factor relates the linear size of the output pixels to the input pixels (scale, s, in F&H paper).
/// drizzle_pixel_fraction controls the linear size of the pixel "drop" in relation to the input (pixfrac, in F&H).

/// so e.g. pixel_fraction =0.6, scale_factor =0.5 are good settings for attemping to double the resolution from a series of interlaced input images.

/// Null translation vector matches lower left bitmap corner (0,0) of input and output
/// translation vector (1,0) matches (0,0) in input to (1/scale,0) in output. (i.e. shift is applied positively to input before dropping onto output)
template<typename input_array_type, typename output_array_type>
void translate_and_drizzle_frame(const PixelArray2d<input_array_type>& input,
        const PixelArray2d<input_array_type>& input_weights,
        PixelArray2d<output_array_type>& output,
        PixelArray2d<output_array_type>& output_weights,
        const PixelShift& translation_vector_at_input_pixel_scale,
        const double drizzle_scale_factor, const double drizzle_pixel_fraction);

template<typename input_array_type>
PixelArray2d<input_array_type> unweight_drizzle_results(
        const PixelArray2d<input_array_type>& output,
        const PixelArray2d<input_array_type>& output_weights);

///Optimization for drizzling both regular and alternate (e.g. thresholded) images.
///NB weights will be the same for both regular and thresholded versions.
template<typename input_array_type, typename output_array_type>
void dual_translate_and_drizzle_frame(
        const PixelArray2d<input_array_type>& input,
        const PixelArray2d<input_array_type>& dual_input,
        const PixelArray2d<input_array_type>& input_weights,
        PixelArray2d<output_array_type>& output,
        PixelArray2d<output_array_type>& dual_output,
        PixelArray2d<output_array_type>& output_weights,
        const PixelShift& translation_vector_at_input_pixel_scale,
        const double drizzle_scale_factor, const double drizzle_pixel_fraction);

//=====================================================================================================================================================================
}//end namespace coela::drizzle
} //end namespace coela

#endif  /* _DRIZZLE_H */

