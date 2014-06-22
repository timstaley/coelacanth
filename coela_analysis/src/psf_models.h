/*
 * File:   psf_models.h
 * Author: ts337
 *
 * Created on 13 May 2010, 13:30
 */

#ifndef COELA_PSF_MODELS_H
#define COELA_PSF_MODELS_H

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/version.hpp>

#include "coela_core/src/cartesian_coords.h"
#include "coela_analysis/src/psf_characterisation.h"

namespace coela {
namespace psf_models {


//==================================================================================
//Interface classes-
//By specifying a common model for all our PSF model classes, we can use  common PSF image generation code,
//and common psf fitting code.

class PsfModelInterface {
public:
    virtual double operator()(const CcdPixelShift& offset_to_psf_centre) const=0;
    virtual ~PsfModelInterface() {}
};

//Converts pixel_offsets into simple radial offsets for axisymmetric models
class AxisymmetricPsfModelInterface: public PsfModelInterface {
public:
    double operator()(const CcdPixelShift& offset_to_psf_centre) const { return (*this)(offset_to_psf_centre.length()); }
    virtual double operator()(double radius_in_CCD_pix) const=
        0; //now only need to define the radial function, and it will get wrapped into the Positional call.
    virtual ~AxisymmetricPsfModelInterface() {}
};
//==================================================================================
//A few analysis routines

///Uses a very crude numerical integration routine - may update to something a bit more sophisticated in the future, but it'll do for now.
double calculate_total_flux(const AxisymmetricPsfModelInterface& psf_gen,
                            const double step_size_in_CCD_pixels=1.0/50.0,
                            const size_t CCD_pixel_limit = 1e4,
                            const double precision=1e-10
                           );

double calculate_flux_enclosed_at_radius(
    const AxisymmetricPsfModelInterface& psf_gen,
    const double radius_in_CCD_pixels,
    const double step_size_in_CCD_pixels=1.0/50.0
);

double fwhm_in_CCD_pix(const AxisymmetricPsfModelInterface& psf_gen,
                       const double CCD_pixel_width_precision=1e-5);
double find_first_minima_in_CCD_pix(const AxisymmetricPsfModelInterface& psf_gen,
                                    const double CCD_pixel_stepsize=0.01);

double calculate_strehl(const AxisymmetricPsfModelInterface& psf_gen,
                        const double unit_peak_airy_total_flux);



//==================================================================================

class UniformPsfModel: public AxisymmetricPsfModelInterface {
public:
    UniformPsfModel();

    UniformPsfModel(double peak_value):peak_val(peak_value) {}

    double operator()(double /*radius_in_CCD_pix*/) const { return peak_val;}

    //data members
    double peak_val;
};

//==================================================================================

class AiryPsfModel: public AxisymmetricPsfModelInterface {
public:
    AiryPsfModel();

    AiryPsfModel(double peak_value, double rads_per_CCD_pixel,
                   double wavelength,
                   double aperture_diameter, double central_obscuration_diameter);

    double operator()(double radius_in_CCD_pix) const;

    //data members
    double peak_val;
    double rads_per_CCD_pix, wavelength, aperture_diameter, obscuration_diameter;
};

//==================================================================================

class GaussianPsfModel: public AxisymmetricPsfModelInterface {
public:
    GaussianPsfModel();
    GaussianPsfModel(double peak_value, double sigma_measured_in_CCD_pixels);

    double operator()(double radius_in_CCD_pix) const;

    //data members
    double peak_val, sigma_in_CCD_pix;

    //--------------------------------------------------------------------------------
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */) {
        ar  & BOOST_SERIALIZATION_NVP(peak_val)
        & BOOST_SERIALIZATION_NVP(sigma_in_CCD_pix);
    }
};
//==================================================================================

///NB Beware - currently assumes the radial values of the profile points correspond to CCD pixels!
class EmpiricalPsfModel: public AxisymmetricPsfModelInterface {
public:
    EmpiricalPsfModel(const std::vector<psf_characterisation::detail::psf_profile_point>&);

    double operator()(double radius_in_CCD_pix) const; //Returns (interpolated) median profile

    psf_characterisation::detail::psf_profile_point get_interpolated_point_at_radius(
        double radius_in_CCD_pix) const;

private:
    std::vector<psf_characterisation::detail::psf_profile_point> datapoints;
};

class HybridRadialSwitchoverPsfModel: public AxisymmetricPsfModelInterface {
public:
    HybridRadialSwitchoverPsfModel(const AxisymmetricPsfModelInterface*
                                       inner_model_ptr,
                                       const AxisymmetricPsfModelInterface* outer_model_ptr,
                                       const double switchover_radius,
                                       const double normalization_factor);
    double operator()(double radius_in_CCD_pix) const;

private:
    const AxisymmetricPsfModelInterface * const inner_model_ptr_;
    const AxisymmetricPsfModelInterface * const outer_model_ptr_;
    const double switchover_radius_;
    const double normalization_factor_;
};

//struct moffat_psf_model: public axisymmetric_psf_base{
//    moffat_psf_model():peak_val(-1.0),sigma_in_CCD_pix(-1.0), beta(-1.0){}
//    moffat_psf_model(double peak_value, double sigma_in_CCD_pixels, double beta_):
//        peak_val(peak_value), sigma_in_CCD_pix(sigma_in_CCD_pixels), beta(beta_){}
//
//    double peak_val, sigma_in_CCD_pix, beta;
//
//
//private:
//    friend class boost::serialization::access;
//    template<class Archive>
//    void serialize(Archive & ar, const unsigned int /* file_version */){
//        ar  & BOOST_SERIALIZATION_NVP(peak_val)
//            & BOOST_SERIALIZATION_NVP(sigma_in_CCD_pix)
//            & BOOST_SERIALIZATION_NVP(beta);
//    }
//};


} //end namespace coela::psf_models
}//end namespace coela

#endif  /* _PSF_MODELS_H */

