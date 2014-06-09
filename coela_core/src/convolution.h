/*
 * File:   frame_registration.h
 * Author: ts337
 *
 * Created on 02 March 2009, 16:53
 */

#ifndef COELA_CONVOLUTION_H
#define COELA_CONVOLUTION_H

#include "coela_core/src/pixel_array2d.h"
#include "coela_core/src/ccd_image.h"

//In fact, these are all applicable to a PixelArray2d

namespace coela {
namespace spatial_convolution {
//===================================================================================================================================================
///Creates an image the same size as the input.
//NB edge regions where the kernel overlay would extend beyond the input region
// are set equal to zero.
template<typename data_type>
PixelArray2d<data_type> convolve_with_kernel(
        const PixelArray2d<data_type>& input,
        const PixelArray2d<data_type>& kernel,
        const PixelIndex& kernel_key_pixel,
        bool output_progress_updates = false);

template<typename data_type>
PixelArray2d<data_type> convolve_with_kernel_at_pixels_above_threshold(
        const PixelArray2d<data_type>& input,
        const PixelArray2d<data_type>& kernel,
        const PixelIndex& kernel_key_pixel, const data_type threshold);

//===================================================================================================================================================
//===================================================================================================================================================
namespace detail {

/// Convolve a kernel with an input image, away from the edges where this is a well defined operation.
/// The "key_pixels" are the pixels that will be Positioned to overlap, when calculating which kernel pixels refer to which input pixels.
///
/// In debug mode, kernel and input overlap is asserted.
/// *WARNING* in release mode no such checking takes place, incorrect usage may result in segfault.
template<typename data_type>
data_type unsafe_convolve_with_kernel_at_Position(
        const PixelArray2d<data_type>& input, const PixelIndex& input_key_pixel,
        const PixelArray2d<data_type>& kernel,
        const PixelIndex& kernel_key_pixel);

PixelRange calculate_valid_region_for_convolution(
        const PixelRange& input_outline, const PixelRange& kernel_outline,
        const PixelIndex& kernel_key_pixel);

} //end namespace coela::convolution::detail
//===================================================================================================================================================
} //end namespace coela::convolution
} //end namespace coela

#endif  /* _CONVOLUTION_H */

