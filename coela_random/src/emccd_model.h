/*
 * File:   CCD_model.h
 * Author: ts337
 *
 * Created on 22 March 2011, 18:17
 */

#ifndef COELA_EMCCD_MODEL_H
#define COELA_EMCCD_MODEL_H

#include <unuran.h>
#include "random.h"
#include "boost_serialization_header_includes.h"
#include <boost/noncopyable.hpp>

namespace coela {
namespace emccd_detail {


unsigned long integer_factorial(int n);
double approx_factorial_double(double number);
double integer_power(double x, int n);

//------------------------------------------------------------------------------------
///Probability mass function for the stochastic EMCCD gain.
///Specifies the probability of an ouput value, given n_photons and gain parameters.
///for high flux gain, use a gaussian ~N( mean = n*g,  var=n*g*g)
/// See Basden et al 2003 for the derivation:
double low_flux_gain_pmf(const int n_photons, const double gain,
                                      const int output_variate_value);


double cicir_pmf(const int N_EM_serial_register_stages,
                                const double gain, const int output_variate_value);

//=======================================================================================
class CicirVariate: private boost::noncopyable_::noncopyable {
public:
    CicirVariate(const int N_EM_serial_register_stages,
                        const double gain,
                        unuran::StreamWrapper&);

    ~CicirVariate();

    int operator()();

private:
    UNUR_GEN *gen;
};
//=======================================================================================

///Helper class; employed by the (variable flux) EMCCD gain variate class.
///Cannot be copied safely ,
/// because we may copy a pointer to a gen which is then destroyed as the owner
/// class goes out of scope.
class FixedLowFluxGainVariate:  private boost::noncopyable_::noncopyable {
public:
    const static int n_photons_limit_for_low_flux=12; //Cannot be more than 20 or so,
    //(factorial becomes too large for unsigned int)

    const static int max_variate_value=10000;


    FixedLowFluxGainVariate(int n_photons, double gain,
                                  unuran::StreamWrapper&);

    ~FixedLowFluxGainVariate();

    int operator()();

private:
//    double lower, upper;
    UNUR_GEN *gen;

};
//=======================================================================================

//Non-copyable again, due to implicitly inherited "EMCCD_fixed_low_flux_gain_variate" objects
class GainVariate: private boost::noncopyable_::noncopyable {
public:
    GainVariate(double gain_init,
                      unuran::StreamWrapper& rns_init);

    ~GainVariate();

    int operator()(int n_photons_input);

private:
    //Bit kludgy  - use pointers to avoid copying / assigning the actual object. It works.
    //Vector index corresponds to n_photons parameter.
    std::vector < FixedLowFluxGainVariate* > fixed_flux_generators;

    unuran::StreamWrapper * rns;
    double gain;
};

//=======================================================================================
} //end namespace coela::CCD_models::detail

//=======================================================================================
//( namespace is coela::CCD_models)

struct EmccdModelParams {
    double bias_pedestal, readout_sigma;
    double photon_gain;
    double serial_CIC_rate;
    int N_EM_serial_register_stages;

    EmccdModelParams() {}
    EmccdModelParams(const std::string& filename);

    void write_to_file(const std::string& filename);

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */) {
        ar
        & BOOST_SERIALIZATION_NVP(bias_pedestal)
        & BOOST_SERIALIZATION_NVP(readout_sigma)
        & BOOST_SERIALIZATION_NVP(photon_gain)
        & BOOST_SERIALIZATION_NVP(serial_CIC_rate)
        & BOOST_SERIALIZATION_NVP(N_EM_serial_register_stages)
        ;
    }
};
std::ostream& operator<<(std::ostream& os, const EmccdModelParams&);


//=======================================================================================


/// A variate to simulate Data Number output from an EMCCD, with stochastic multiplication, CIC, and readout noise.
///  NB the EMCCD_model does not have a function to "Simulate one pixel value" because this would be quite inefficient.
///  Instead, the different sub-process simulators are called per-image from the function:
///  "image_noise_simulation::EMCCD_simulated_image(...)"

class EmccdModel {
public:

    EmccdModel(const EmccdModelParams& parameter_init,
                unuran::StreamWrapper& rns_in);

    const EmccdModelParams params;

    int stochastically_multiply_photons(int n_photons_input) {
        return photon_gain_variate(n_photons_input);}
    int CICIR_variate() {return CICIR_variate_();}
    int readout_noise_variate();

private:
    emccd_detail::GainVariate photon_gain_variate;
    unuran::GaussianRandomVariate readout_variate;
    emccd_detail::CicirVariate CICIR_variate_;
};



}//end namespace coela


#endif  /* EMCCD_MODEL_H */

