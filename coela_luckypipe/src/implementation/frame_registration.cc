#include "../frame_registration.h"
#include "coela_core/src/image_utils.h"
#include "coela_core/src/convolution.h"
#include <utility>
#include <stdexcept>
#include <iostream>
#include <coela_utility/src/misc_math.h>
using std::cerr; using std::endl;
using std::pair;
namespace coela {
namespace frame_registration {

using namespace psf_models;
//==============================================================================================================

psf_models::reference_psf convert_to_log_values(const psf_models::reference_psf& in)
{
    psf_models::reference_psf copy(in);
    if (copy.psf_image.pix(copy.psf_image.pix.min_PixelIndex()) <0)
        throw std::logic_error("ref_psf::convert_to_log_vals - "
                               "Not recommended for reference PSF models which have negative values");
    copy.psf_image.pix+=1.0;
    for (PixelIterator i(copy.psf_image.pix.range()); i!=i.end; ++i) {
        copy.psf_image.pix(i)=log(copy.psf_image.pix(i));
    }
    return copy;
}
//==============================================================================================================
template<typename input_datatype>
gs_lock<CCD_Position> find_best_psf_match(
    const CCDImage<input_datatype>& input_image,
    const CCD_BoxRegion& gs_region,
    const psf_models::reference_psf& ref_psf,
    const double input_resample_factor,
    const double convolution_threshold_factor,
    const bool fit_parabola)
{

    CCDImage<double> resampled_sub_image;
    if (input_resample_factor!=1.0) {
        resampled_sub_image=
            image_utils::bicubic_resample_non_edge_region(
                input_image , gs_region, input_resample_factor
            );
    } else {
        resampled_sub_image=CCDImage<double>(
                                CCDImage<input_datatype>::sub_image(input_image, gs_region));
    }

    if (resampled_sub_image.CCD_grid.pixel_width_ !=
            ref_psf.psf_image.CCD_grid.pixel_width_) {
//        cerr<<"input image pixel scale:" << input_image.CCD_grid.pixel_width_<<endl;
//        cerr<<"resampled image pixel scale:" << resampled_sub_image.CCD_grid.pixel_width_<<endl;
//        cerr<<"ref psf pixel scale:" << ref_psf.psf_image.CCD_grid.pixel_width_<<endl;

        throw std::runtime_error("find_best_psf_match: - reference psf pixel scale does not match resampled input;");
    }

    CCDImage<double> convolution;
    if (convolution_threshold_factor==0.0) {
        convolution.pix=
            spatial_convolution::convolve_with_kernel(
                resampled_sub_image.pix, ref_psf.psf_image.pix, ref_psf.central_pixel) ;
    } else {
        double threshold(resampled_sub_image.pix(resampled_sub_image.pix.max_PixelIndex()) *
                         convolution_threshold_factor);
        convolution.pix =
            spatial_convolution::convolve_with_kernel_at_pixels_above_threshold(
                resampled_sub_image.pix, ref_psf.psf_image.pix, ref_psf.central_pixel,
                threshold) ;
    }

    PixelShift ref_psf_centre_offset =
        ref_psf.exact_centre - PixelPosition::centre_of_pixel(ref_psf.central_pixel);

    PixelIndex conv_peak_pixel = convolution.pix.max_PixelIndex();
    double conv_peak = convolution.pix(conv_peak_pixel);
    PixelPosition conv_peak_position = PixelPosition::centre_of_pixel(conv_peak_pixel);

    //This next bit gets a bit hairy and should probably be functioned off with its own unit test.
    if (fit_parabola) {
        //Check the max isn't on the edge of the valid region
        PixelRange valid_conv_output_region =
            spatial_convolution::detail::calculate_valid_region_for_convolution(
                resampled_sub_image.pix.range(),
                ref_psf.psf_image.pix.range(), ref_psf.central_pixel);
        if (PixelRange::pad(valid_conv_output_region, -1).
                contains_pixel(conv_peak_pixel)) {
            //ok, we're in business:
            //First, we make sure that the pixels above and below have been convolved:
            vector<PixelIndex> adjacent_pixels;
            adjacent_pixels.push_back(PixelIndex(conv_peak_pixel.x-1, conv_peak_pixel.y));
            adjacent_pixels.push_back(PixelIndex(conv_peak_pixel.x+1, conv_peak_pixel.y));
            adjacent_pixels.push_back(PixelIndex(conv_peak_pixel.x, conv_peak_pixel.y-1));
            adjacent_pixels.push_back(PixelIndex(conv_peak_pixel.x, conv_peak_pixel.y+1));

            bool bad_convolution_threshold=false;
            for (size_t i=0; i!=adjacent_pixels.size(); ++i) {
                double adjacent_pixel_conv_value =
                    spatial_convolution::detail::
                    unsafe_convolve_with_kernel_at_Position(
                        resampled_sub_image.pix,
                        adjacent_pixels[i],
                        ref_psf.psf_image.pix,
                        ref_psf.central_pixel
                    );
                if (adjacent_pixel_conv_value > conv_peak) {
                    bad_convolution_threshold=true;
                    conv_peak = adjacent_pixel_conv_value;
                    conv_peak_pixel = adjacent_pixels[i];
                }
                convolution.pix(adjacent_pixels[i]) = adjacent_pixel_conv_value;
            }

            if (bad_convolution_threshold) {

            } else {
                //Now we can perform parabolic interpolation, first in the x direction:
                pair<double, double> x_interp_pos_val =
                    misc_math::fit_1d_parabola_to_find_local_maxima(
                        conv_peak_position.x,
                        convolution.pix(PixelIndex(conv_peak_pixel.x-1, conv_peak_pixel.y)),
                        conv_peak,
                        convolution.pix(PixelIndex(conv_peak_pixel.x+1, conv_peak_pixel.y))
                    );

                //And then y direction
                pair<double, double> y_interp_pos_val =
                    misc_math::fit_1d_parabola_to_find_local_maxima(
                        conv_peak_position.y,
                        convolution.pix(PixelIndex(conv_peak_pixel.x, conv_peak_pixel.y-1)),
                        conv_peak,
                        convolution.pix(PixelIndex(conv_peak_pixel.x, conv_peak_pixel.y+1))
                    );
                PixelPosition parabola_peak_posn(x_interp_pos_val.first, y_interp_pos_val.first);
                double parabola_peak_val = 0.5*(x_interp_pos_val.second + y_interp_pos_val.second);

                if (parabola_peak_val > conv_peak &&
                        coord_distance(parabola_peak_posn, conv_peak_position)<1.5) {
                    conv_peak = parabola_peak_val;
                    conv_peak_position = parabola_peak_posn;
                }
            }
        }//Else... Peak pixel is at the edge of the valid convolution region, so ignore and move on

    }//end of fit parabola section


    PixelPosition subimage_psf_lock = conv_peak_position + ref_psf_centre_offset;

    CCD_Position psf_lock_posn = resampled_sub_image.CCD_grid.corresponding_grid_Position(
                                     subimage_psf_lock);

    return gs_lock<CCD_Position>(psf_lock_posn, conv_peak);
}


template
gs_lock<CCD_Position> find_best_psf_match(
    const CCDImage<float>& input_image,
    const CCD_BoxRegion& gs_region,
    const psf_models::reference_psf& ref_psf,
    const double input_resample_factor,
    const double convolution_threshold_factor,
    const bool fit_parabola);
template
gs_lock<CCD_Position> find_best_psf_match(
    const CCDImage<double>& input_image,
    const CCD_BoxRegion& gs_region,
    const psf_models::reference_psf& ref_psf,
    const double input_resample_factor,
    const double convolution_threshold_factor,
    const bool fit_parabola);
//==============================================================================================================
}
}
