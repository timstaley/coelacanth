/*
 * File:   psf_generators.cc
 * Author: tim
 *
 * Created on 06 January 2010, 21:25
 */

#include "../psf_generation.h"
#include "coela_utility/src/misc_math.h"
#include "coela_core/src/image_utils.h"
#include "coela_core/src/pixel_array_routines.h"
#include "coela_utility/src/string_utils.h"
#include <stdexcept>
#include <iostream>
#include <valarray>
using namespace std;

namespace coela {
namespace psf_models {

//==================================================================================================================
vector<double> reference_psf::get_unmasked_ref_psf_pixel_values() const
{
    return pixel_array_routines::get_region_unmasked_pixels_vector(psf_image.pix,
            mask.pix, psf_image.pix.range());
}
//==================================================================================================================

reference_psf generate_psf(const psf_model_interface& psf,
                           const PixelRange& output_frame_shape,
                           const CCD_Position& output_CCD_region_low_corner,
                           const CCD_Position& psf_centre_point,
                           const double desired_pixel_scale_relative_to_CCD_pix,
                           const double oversampling,
                           const bool radius_limited,
                           const double max_CCD_radius_reqd
                          )
{
    if (fmod(oversampling, 1.0) != 0.0) {
        throw runtime_error("Generate_psf_with_params: Please supply an integer oversample factor");
    }
    //NB - could branch and use either "shift_and_rebin" or "pre-aligned downsample" depending on whether this is integer or rational, if needed. int only for now.

    if (radius_limited && max_CCD_radius_reqd==0.0) {
        throw logic_error("generate_psf_with_params - Should supply a max radius if radius limit is switched on!");
    }

    double max_CCD_generation_radius = max_CCD_radius_reqd + 1.0;
    //we increase the limit by a full CCD pixel, in case the centre of a sub-pixel is outside the radius, but the centre of the output pixel it contributes to is not.

    CCDImage<double> construction_image;
    construction_image.pix = PixelArray2d<double>(ceil(output_frame_shape.x_dim() *
                             oversampling),
                             ceil(output_frame_shape.y_dim() * oversampling),
                             0.0);

    construction_image.initialize_CCD_grid_to_specific_offset_and_scale(
        output_CCD_region_low_corner,
        desired_pixel_scale_relative_to_CCD_pix / oversampling);

    PixelPosition construction_centre =
        construction_image.CCD_grid.corresponding_pixel_Position(psf_centre_point);

    double construction_CCD_scale =
        construction_image.CCD_grid.pixel_width_;

    if (radius_limited==false) {
        for (PixelIterator img_it(construction_image.pix.range()); img_it!=img_it.end; ++img_it) {
            PixelPosition pixel_centre=PixelPosition::centre_of_pixel(img_it);
            PixelShift pixel_offset_frm_coords(pixel_centre - construction_centre);
            CCD_PixelShift pixel_offset_CCD_coords(
                pixel_offset_frm_coords.x*construction_CCD_scale    ,
                pixel_offset_frm_coords.y*construction_CCD_scale);
            construction_image.pix(img_it) = psf(pixel_offset_CCD_coords);
        }
    } else { // if (radius_limited==true) //exactly the same loop with an if (dist<max_radius) added.
        for (PixelIterator img_it(construction_image.pix.range()); img_it!=img_it.end; ++img_it) {
            PixelPosition pixel_centre=PixelPosition::centre_of_pixel(img_it);
            PixelShift pixel_offset_frm_coords(pixel_centre - construction_centre);
            CCD_PixelShift pixel_offset_CCD_coords(
                pixel_offset_frm_coords.x*construction_CCD_scale    ,
                pixel_offset_frm_coords.y*construction_CCD_scale);
            if ((pixel_offset_CCD_coords.length()) < max_CCD_generation_radius) {
                construction_image.pix(img_it) = psf(pixel_offset_CCD_coords);
            }
        }

    }

    reference_psf output_psf;
    if (oversampling!=1.0) {
        output_psf.psf_image = image_utils::bin_image(construction_image, oversampling);
        output_psf.psf_image.pix /= (oversampling*oversampling);
    } else { output_psf.psf_image = construction_image; }
    assert(output_psf.psf_image.CCD_grid.image_outline_.low ==
           output_CCD_region_low_corner); //sanity check
    assert(output_psf.psf_image.CCD_grid.pixel_width_==
           desired_pixel_scale_relative_to_CCD_pix); //sanity check

    output_psf.mask=CCDImage<double>(output_psf.psf_image);
    output_psf.mask.pix.assign(1.0);//initialised to same size, all=1.0; (unmasked)

    output_psf.exact_centre = output_psf.psf_image.CCD_grid.corresponding_pixel_Position(
                                  psf_centre_point);
    output_psf.central_pixel = PixelPosition::pixel_containing_point(output_psf.exact_centre);

    return output_psf;
}

reference_psf generate_centred_psf(const axisymmetric_psf_model_interface& psf_func,
                                   const double outer_radius_in_CCD_pix,
                                   const double desired_pixel_scale_relative_to_CCD_pix,
                                   const int oversampling_factor
                                  )
{
    //This is crucial - always keep oversampling_factor outside the float bit.
    //as we're using a simple rebin fudge here.
    double zoom_factor =1.0 / desired_pixel_scale_relative_to_CCD_pix;
    size_t output_size = (ceil(outer_radius_in_CCD_pix) * zoom_factor*2 +1);
    //odd number of kernel_pixels (after rebinning) means that there is a high value central pixel with symmetric kernel around. (For better Xcorr locking)
    PixelRange output_outline(1,1,output_size,output_size);

//    ref_psf.psf_image = image<double>(output_size, output_size, 0.0);
//    ref_psf.psf_image.initialize_CCD_grid_to_specific_offset_and_scale(CCD_Position(0,0) , desired_pixel_scale_relative_to_CCD_pix);

    //NB output size is always odd, so this places the PSF centre at a pixel centre.
    CCD_Position psf_centre =
        CCD_Position(output_size/2.0 * desired_pixel_scale_relative_to_CCD_pix,
                     output_size/2.0 * desired_pixel_scale_relative_to_CCD_pix) ;


    reference_psf ref_psf=
        generate_psf(
            psf_func,
            output_outline,
            CCD_Position(0,0),
            psf_centre,
            desired_pixel_scale_relative_to_CCD_pix,
            oversampling_factor,
            true,
            outer_radius_in_CCD_pix);

    PixelPosition psf_pixel_centre =
        ref_psf.psf_image.CCD_grid.corresponding_pixel_Position(psf_centre);
    double outer_radius_in_scaled_pix =
        outer_radius_in_CCD_pix/desired_pixel_scale_relative_to_CCD_pix;

    for (PixelIterator it(ref_psf.mask.pix.range()); it!=it.end; ++it) {
        if (coord_distance(PixelPosition::centre_of_pixel(it), psf_pixel_centre) >
                outer_radius_in_scaled_pix) {
            ref_psf.mask.pix(it)=0.0;
        }
    }
    return ref_psf;
}



//==================================================================================


//
//gauss_moffat_psf_generator::gauss_moffat_psf_generator(const gauss_psf_generator& g_component, const moffat_psf_generator& moff_component){
//    gauss=g_component;
//    moffat=moff_component;
//    combined_peak_val=gauss.peak_val+moffat.peak_val;
//    gauss_weight=gauss.peak_val/combined_peak_val;
//}
//
//void gauss_moffat_psf_generator::set_peak_value(double combined_peak_value){
//    combined_peak_val=combined_peak_value;
//    derive_and_set_component_peaks();
//}
//
//void gauss_moffat_psf_generator::set_gauss_component_weight(double gauss_peak_weight){
//    assert(gauss_peak_weight>=0 && gauss_peak_weight <=1.0);
//    gauss_weight=gauss_peak_weight;
//    derive_and_set_component_peaks();
//}
//
//void gauss_moffat_psf_generator::derive_and_set_component_peaks(){
//    assert(gauss_weight>=0 && gauss_weight <=1.0);
//    gauss.peak_val=gauss_weight*combined_peak_val;
//    moffat.peak_val=(1.0-gauss_weight)*combined_peak_val;
//}
//
//double gauss_moffat_psf_generator::operator()(double radius_in_CCD_pix) const{
//    return moffat(radius_in_CCD_pix)+ gauss(radius_in_CCD_pix);
//}
//==================================================================================================================
} //end namespace coela::psf_models
}//end namespace coela
