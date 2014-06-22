/*
 * File:   PsfFitting.cc
 * Author: ts337
 *
 * Created on 13 May 2011, 16:34
 */

#include <memory>

#include "../psf_fitting.h"
#include "coela_analysis/src/psf_generation.h"

#include <iostream>
using std::cerr;
using std::endl;

namespace coela {
namespace psf_fitting {

//==================================================================================
AnalyticPsfFitCostFCN::AnalyticPsfFitCostFCN(const CcdImage<double>& data,
        const CcdImage<double>& fitting_mask,
        const Mn2_models::MinuitModelInterface& model,
        const variance_estimation_base& noise_model,
        const int desired_oversampling_level)
    : input_img_ref(data), mask(fitting_mask),
      psf_model(model), noise_estimator(noise_model),
      oversampling_level(desired_oversampling_level)
{}

//------------------------------------------------------------------------------------------
double AnalyticPsfFitCostFCN::operator()(
    const std::vector<double>& model_pars) const
{

    CcdImage<double> model_img = generate_model_image(model_pars);

    return reduced_chi_squared(input_img_ref, mask,
                               model_img,
                               noise_estimator,
                               model_pars.size());
}
//------------------------------------------------------------------------------------------
CcdImage<double> AnalyticPsfFitCostFCN::generate_model_image(
    const std::vector<double>& fitting_pars) const
{
    assert(fitting_pars.size()>=2);

    CcdPosition estimated_psf_centre(fitting_pars[0], fitting_pars[1]);
    vector<double> model_pars(fitting_pars.begin()+2, fitting_pars.end());
    Mn2_models::ModelParamsWrapper bound_model(psf_model, model_pars);

    return psf_models::generate_psf(bound_model, mask.pix.range(),
                                    mask.CCD_grid.image_outline_.low,
                                    estimated_psf_centre,
                                    mask.CCD_grid.pixel_width_,
                                    oversampling_level).psf_image;

}
//==================================================================================
double reduced_chi_squared(const CcdImage<double>& input,
                           const CcdImage<double>& input_mask,
                           const CcdImage<double>& fitting_model,
                           const variance_estimation_base& noise_estimator,
                           const int num_fitted_parameters)
{

    assert(input_mask.pix.range()==fitting_model.pix.range());
    assert(input_mask.CCD_grid.pixel_width_ ==
           fitting_model.CCD_grid.pixel_width_);

    //NB Such verbosity *is* necessary if we are dealing with unknown zoom and offset levels
    CcdPosition input_mask_low_pixel_centre =
        input_mask.CCD_grid.corresponding_grid_Position(
            PixelPosition::centre_of_pixel(PixelIndex(1,1))
        );
    PixelIndex input_pixel_matching_mask_low_corner = PixelPosition::pixel_centred_at(
                input.CCD_grid.corresponding_pixel_Position(input_mask_low_pixel_centre));
//    cerr<<"Matching input pixel: " <<input_pixel_matching_mask_low_corner<<endl;

    PixelIndex input_pixels_offset = input_pixel_matching_mask_low_corner - PixelIndex(1,1);

    double sum=0.0;
    int n_fitted_data_points=0;

    for (PixelIterator model_index(fitting_model.pix.range());
            model_index!=model_index.end; ++model_index) {
        if (input_mask.pix(model_index)!=0.0) {
            double pixvar = noise_estimator(model_index + input_pixels_offset);
            double d = input.pix(model_index + input_pixels_offset) - fitting_model.pix(model_index);
            d*=d; //diff_sqrd
            d/=pixvar;
            sum+=d;
            ++n_fitted_data_points;
        }
    }

    int degs_of_freedom = n_fitted_data_points - num_fitted_parameters - 1;
    return sum / degs_of_freedom;
}


} //end namespace coela::PsfFitting
}//end namespace coela

