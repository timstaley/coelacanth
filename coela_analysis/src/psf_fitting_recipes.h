/*
 * File:   psf_fitting_recipes.h
 * Author: ts337
 *
 * Created on 16 May 2011, 12:01
 */

#ifndef COELA_PSF_FITTING_RECIPES_H
#define COELA_PSF_FITTING_RECIPES_H

#include "psf_fitting.h"

namespace coela {
namespace psf_fitting {

Mn2_models::PSF_Fit<psf_models::gaussian_psf_model> fit_gaussian_model(
    const CCDImage<double>& background_subtracted_image,
    const CCD_Position& estimated_centre,
    const variance_estimation_base& noise_model,
    const double fitting_radius_in_CCD_pix,
    const int oversampling,
    bool verbose_display);

} //end namespace coela::psf_fitting
}//end namespace coela


#endif  /* PSF_FITTING_RECIPES_H */

