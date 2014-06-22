/*
 * File:   PsfFitting_recipes.h
 * Author: ts337
 *
 * Created on 16 May 2011, 12:01
 */

#ifndef COELA_PsfFitTING_RECIPES_H
#define COELA_PsfFitTING_RECIPES_H

#include "psf_fitting.h"

namespace coela {
namespace psf_fitting {

Mn2_models::PsfFit<psf_models::GaussianPsfModel> fit_gaussian_model(
    const CcdImage<double>& background_subtracted_image,
    const CcdPosition& estimated_centre,
    const variance_estimation_base& noise_model,
    const double fitting_radius_in_CCD_pix,
    const int oversampling,
    bool verbose_display);

} //end namespace coela::PsfFitting
}//end namespace coela


#endif  /* PsfFitTING_RECIPES_H */

