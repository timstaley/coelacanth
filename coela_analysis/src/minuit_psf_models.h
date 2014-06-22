/*
 * File:   minuit_psf_models.h
 * Author: ts337
 *
 * Created on 13 May 2011, 18:26
 */

#ifndef COELA_MINUIT_PSF_MODELS_H
#define COELA_MINUIT_PSF_MODELS_H
#include "coela_analysis/src/psf_models.h"
#include <Minuit2/MnUserParameters.h>
namespace coela {

namespace psf_fitting {
namespace Mn2_models {
//==================================================================================
//interface for stateless PSF model functions for generation with parameters varying

class MinuitModelInterface {
public:
    virtual double operator()(const CcdPixelShift& offset_from_psf_centre,
                              const std::vector<double>& psf_params) const = 0;
};

//==================================================================================
//Takes a reference of a minuit psf model, and the current set of params,
//And wraps them together for compatibility with psf_model_base (so we can re-use psf generation code)
class ModelParamsWrapper: public psf_models::PsfModelInterface {
public:
    ModelParamsWrapper(const MinuitModelInterface& parameterised_model,
                         const std::vector<double>& psf_params);

    //Override the psf_model_base function:
    double operator()(const CcdPixelShift& offset_from_psf_centre) const;
private:
    const MinuitModelInterface& model_ref;
    const std::vector<double> pars_copy;
};

//==================================================================================
struct PositionFit {
    CcdPosition centre;
    CcdPixelShift error_estimate;

    static void set_MnPars(const PositionFit &,
                           ROOT::Minuit2::MnUserParameters * ,
                           const std::string& prefix="unnamed_model");
    //------------------------------------------------------------
    static PositionFit pull_from_MnPars(const ROOT::Minuit2::MnUserParameters& ,
                                        const std::string& prefix="unnamed_model");
private:
    const static std::string centre_x_suffix, centre_y_suffix;
};

template<class PsfType>
struct PsfFit {
    PsfType model;
    PositionFit Position;
};

//==================================================================================

class GaussianPsf: public MinuitModelInterface {
public:
    double operator()(const CcdPixelShift& offset_from_psf_centre,
                      const std::vector<double>& psf_params) const;

    static void set_MnPars_from_PsfFit(
        const PsfFit<psf_models::GaussianPsfModel> initial_fit ,
        ROOT::Minuit2::MnUserParameters * ,
        const std::string& prefix);

    static PsfFit<psf_models::GaussianPsfModel> pull_PsfFit_from_MnPars(
        const ROOT::Minuit2::MnUserParameters& fitted_pars,
        const std::string& prefix);

private:
    const static std::string model_name, peak_suffix, sigma_suffix;
    //Params: [peak, sigma]
};



}//end namespace coela::PsfFitting::minuit_psf_models
} //end namespace coela::PsfFitting
}//end namespace coela

#endif  /* MINUIT_PSF_MODELS_H */

