/*
 * File:   psf_fitting.h
 * Author: ts337
 *
 * Created on 13 May 2011, 16:34
 */

#ifndef COELA_PsfFitTING_H
#define COELA_PsfFitTING_H

#include "coela_analysis/src/psf_generation.h"

#include "minuit_psf_models.h"
#include "fitting_noise_models.h"
#include <Minuit2/FCNBase.h>

namespace coela {
namespace psf_fitting {

//To do - implement image<bool> for masking purposes

//==================================================================================
//Cost function to be fed into Minuit minimizer routines.

class AnalyticPsfFitCostFCN: public ROOT::Minuit2::FCNBase {
public:

    ///Takes a reference to input data, data model, and noise_model.
    ///NB does not take ownership, so ALL references must remain valid whilst in use!
    ///Therefore, don't instantiate e.g. a temporary least_squares_noise_model() on the parameter input line.

    /// The cost function is agnostic about how the fit region is described.
    /// You simply pass it a "Fitting mask" image.
    /// The model psf is generated over the same area as the fitting_mask pixel array, so the mask should not be unnecessarily large or performance will be slow.
    /// Then the model and data are compared over all pixels for which the corresponding mask values are non-zero.

    //p.s. Oversampling levels above 1 will slow down proceedings considerably.
    AnalyticPsfFitCostFCN(const CcdImage<double>& data,
                              const CcdImage<double>& fitting_mask,
                              const Mn2_models::MinuitModelInterface& model,
                              const variance_estimation_base& noise_model,
                              const int desired_oversampling_level=1
                             );
    ///<TO DO - change these args to pointers, since this better fits / implies their usage.

    double operator()(const std::vector<double>& model_pars) const;
    double Up() const {return 1.0;} ///<Required by Minuit, see docs re: "ErrorDef"

private:
    CcdImage<double> generate_model_image(const std::vector<double>& fitting_pars) const;

    //TO DO - change these to pointers, since this better fits / implies their usage.
    const CcdImage<double>& input_img_ref;
    const CcdImage<double>&  mask;
    const Mn2_models::MinuitModelInterface& psf_model;
    const variance_estimation_base& noise_estimator;
    const int oversampling_level;

};

//==================================================================================
//Sub-routines:
double reduced_chi_squared(const CcdImage<double>& input,
                           const CcdImage<double>& input_mask,
                           const CcdImage<double>& fitting_model,
                           const variance_estimation_base& noise_estimator,
                           const int num_fitted_parameters);


//==================================================================================




} //end namespace coela::PsfFitting
}//end namespace coela

#endif  /* PsfFitTING_H */

