/*
 * File:   minuit_psf_models.cc
 * Author: ts337
 *
 * Created on 13 May 2011, 18:26
 */

#include "../minuit_psf_models.h"
#include "coela_utility/src/misc_math.h"
namespace coela {
namespace psf_fitting {
namespace Mn2_models {
using std::vector;

//==================================================================================
ModelParamsWrapper::ModelParamsWrapper(
    const MinuitModelInterface& parameterised_model,
    const vector<double>& params)
    : model_ref(parameterised_model),
      pars_copy(params)
{}

double ModelParamsWrapper::operator()(const CcdPixelShift& offset) const
{
    return model_ref(offset, pars_copy);
}

//==================================================================================
const std::string PositionFit::centre_x_suffix("_centre_CCD_x"),
      PositionFit::centre_y_suffix("_centre_CCD_y");

void PositionFit::set_MnPars(const PositionFit & p,
                             ROOT::Minuit2::MnUserParameters * pars_ptr,
                             const std::string& prefix)
{
    pars_ptr->Add(prefix+centre_x_suffix, p.centre.x, p.error_estimate.x);
    pars_ptr->Add(prefix+centre_y_suffix, p.centre.y, p.error_estimate.y);
    return;
}

//------------------------------------------------------------
PositionFit PositionFit::pull_from_MnPars(const ROOT::Minuit2::MnUserParameters& pars,
        const std::string& prefix)
{
    PositionFit fit;
    fit.centre.x = pars.Value(prefix+centre_x_suffix);
    fit.centre.y = pars.Value(prefix+centre_y_suffix);
    fit.error_estimate.x = pars.Error(prefix+centre_x_suffix);
    fit.error_estimate.y = pars.Error(prefix+centre_y_suffix);
    return fit;
}


//==================================================================================

const std::string GaussianPsf::model_name("Gauss"),
      GaussianPsf::peak_suffix("_peak"),
      GaussianPsf::sigma_suffix("_sigma");
//------------------------------------------------------------
double GaussianPsf::operator()(const CcdPixelShift& offset,
                                const std::vector<double>& psf_params) const
{
    assert(psf_params.size()==2);
    double radius = offset.length();
    return misc_math::gaussian_1d_function(radius, psf_params[0], psf_params[1]);
}
//-------------------------------------------------------------------------------
void GaussianPsf::set_MnPars_from_PsfFit(
    const PsfFit<psf_models::GaussianPsfModel> initial_fit,
    ROOT::Minuit2::MnUserParameters* pars_ptr,
    const std::string& prefix)
{
    PositionFit::set_MnPars(initial_fit.Position,
                            pars_ptr,
                            prefix+model_name);
    pars_ptr->Add(prefix + model_name + peak_suffix,
                  initial_fit.model.peak_val,
                  initial_fit.model.peak_val * 0.25);
    pars_ptr->SetLowerLimit(prefix + model_name + peak_suffix, 0.0);

    pars_ptr->Add(prefix + model_name + sigma_suffix,
                  initial_fit.model.sigma_in_CCD_pix,
                  initial_fit.model.sigma_in_CCD_pix * 0.5);
    pars_ptr->SetLowerLimit(prefix + model_name + sigma_suffix, 1e-3);

}
//-------------------------------------------------------------------------------
PsfFit<psf_models::GaussianPsfModel> GaussianPsf::pull_PsfFit_from_MnPars(
    const ROOT::Minuit2::MnUserParameters& pars,
    const std::string& prefix)
{
    PsfFit<psf_models::GaussianPsfModel> fit;

    fit.Position = PositionFit::pull_from_MnPars(pars, prefix+model_name);

    fit.model.peak_val = pars.Value(prefix+model_name+peak_suffix);
    fit.model.sigma_in_CCD_pix = pars.Value(prefix+model_name+sigma_suffix);

    return fit;
}
//==================================================================================


}//end namespace coela::PsfFitting::minuit_psf_models
} //end namespace coela::PsfFitting
}//end namespace coela
