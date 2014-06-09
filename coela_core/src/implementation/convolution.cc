#include <math.h>
#include <vector>
#include "../convolution.h"

#include "coela_core/src/image_utils.h"


#include <stdexcept>
#include <iostream>

using std::cout; using std::cerr; using std::endl;
using std::vector;
using namespace std;

namespace coela {
namespace spatial_convolution {
//===================================================================================================================================================

template< typename data_type>
PixelArray2d<data_type> convolve_with_kernel(const PixelArray2d<data_type>& input,
        const PixelArray2d<data_type>& kernel, const PixelIndex& kernel_key_pixel,
        bool output_progress_updates)
{
    PixelRange output_region = detail::calculate_valid_region_for_convolution(
                                   input.range(),
                                   kernel.range(), kernel_key_pixel);

    PixelArray2d<data_type> output(input.range().x_dim(), input.range().y_dim(), 0.0);

    if (output_progress_updates) {
        cout<<"Convolving... completion at 1.0..."<<endl;
        size_t npix = output_region.n_pix();
        size_t counter=0;
        for (PixelIterator i(output_region); i!=i.end; ++i) {
            output(i) =
                detail::unsafe_convolve_with_kernel_at_Position(input, i,
                        kernel, kernel_key_pixel);
            counter++;
            if (counter%100==0) { cout<<"\r" <<double(counter)/double(npix)<<flush; }
        }
    } else {
        for (PixelIterator i(output_region); i!=i.end; ++i) {
            output(i) =
                detail::unsafe_convolve_with_kernel_at_Position(input, i,
                        kernel, kernel_key_pixel);
        }
    }
    cout<<endl;
    return output;
}

template PixelArray2d<double> convolve_with_kernel(
    const PixelArray2d<double>& input,
    const PixelArray2d<double>& kernel, const PixelIndex& kernel_key_pixel,
    bool);


template<typename data_type>
PixelArray2d<data_type> convolve_with_kernel_at_pixels_above_threshold(
    const PixelArray2d<data_type>& input,
    const PixelArray2d<data_type>& kernel, const PixelIndex& kernel_key_pixel,
    const data_type threshold)
{

    PixelRange output_region = detail::calculate_valid_region_for_convolution(
                                   input.range(),
                                   kernel.range(), kernel_key_pixel);

    PixelArray2d<data_type> output(input.range().x_dim(), input.range().y_dim(), 0.0);
    for (PixelIterator i(output_region); i!=i.end; ++i) {
        if (input(i)> threshold) {
            output(i) =
                detail::unsafe_convolve_with_kernel_at_Position(input, i,
                        kernel, kernel_key_pixel);
        }
    }
    return output;
}

template PixelArray2d<double> convolve_with_kernel_at_pixels_above_threshold(
    const PixelArray2d<double>& ,
    const PixelArray2d<double>& , const PixelIndex& , const double threshold);


//===================================================================================================================================================
namespace detail {

template<typename data_type>
data_type unsafe_convolve_with_kernel_at_Position(
    const PixelArray2d<data_type>& input, const PixelIndex& input_key_pixel,
    const PixelArray2d<data_type>& kernel, const PixelIndex& kernel_key_pixel)
{
    PixelIndex input_zero_index(input_key_pixel - kernel_key_pixel);

    assert(input.range().contains_range(
               PixelRange(input_zero_index + PixelIndex(1,1),
                          input_zero_index + kernel.range().high)
           ));

//compute convolution result
    data_type sum(0.0);
    for (PixelIterator k_index(kernel.range()); k_index!=k_index.end; ++k_index) {
        sum +=input(input_zero_index + k_index) * kernel(k_index); // * kernel.data_mask(k_index);
    }
    return sum;
}

template double unsafe_convolve_with_kernel_at_Position(
    const PixelArray2d<double>& input, const PixelIndex& input_key_pixel,
    const PixelArray2d<double>& kernel, const PixelIndex& kernel_key_pixel);



PixelRange calculate_valid_region_for_convolution(
        const PixelRange& input_outline,
        const PixelRange& kernel_outline, const PixelIndex& kernel_key_pixel)
{
    assert(input_outline.is_valid());
    assert(kernel_outline.is_valid());
    assert(kernel_key_pixel.is_valid());

    PixelIndex kernel_lower_corner_offset = kernel_key_pixel - kernel_outline.low;
    PixelIndex kernel_upper_corner_offset = kernel_outline.high - kernel_key_pixel;
    return PixelRange(input_outline.low  + kernel_lower_corner_offset,
                      input_outline.high - kernel_upper_corner_offset);
}

} //end namespace coela::convolution::detail
//===================================================================================================================================================


//===================================================================================================================================================
} //end namespace coela::convolution
}//end namespace coela
