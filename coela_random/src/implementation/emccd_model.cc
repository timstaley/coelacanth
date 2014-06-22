
#include "../emccd_model.h"
#include "../boost_serialization_body_includes.h"
#include <cassert>

using namespace std;
using namespace coela;

namespace coela {
namespace emccd_detail {

unsigned long integer_factorial(int n)
{
    assert(n<20);
    if (n==0) { return 1ul; }
    return n*integer_factorial(n-1);
}


///Sidesteps the integer limits, but as such is not exact. Used for calculating PDFs for an EMCCD (which are already approximations anyway).
double approx_factorial_double(double number)
{
    if (number <= 1.0) { return 1.0; }
    return number * approx_factorial_double(number - 1.0);
}

double integer_power(double x, int n)
{
    if (n==0) { return 1.0; }
    if (n<0) { return 1.0/integer_power(x,-1*n); }

    double result=x;
    while (n!=1) {
        result*=x;
        n--;
    }
    return result;
}



double low_flux_gain_pmf(const int n_photons, const double gain,
                                      const int output_variate_value)
{
    assert(n_photons<12);
    assert(n_photons>0);

    if (output_variate_value < n_photons) { return 0.0; }

    double denom =
        integer_power(gain, n_photons)* integer_factorial(n_photons-1);

    double numerator = integer_power(output_variate_value, n_photons-1)
                       * exp(-output_variate_value /gain);

    return numerator/denom;
}

double cicir_pmf(const int num_serial_register_stages,
                                const double gain, const int output_variate_value)
{
    if (output_variate_value < 0) { return 0; }

    const int & N_stages = num_serial_register_stages;
    const double prob_per_stage = pow(gain, 1.0/N_stages);
    const double norm = 1.0 / N_stages;
    double output_prob=0.0;
    for (int CICIR_stages_passed = 1; CICIR_stages_passed<=N_stages; ++CICIR_stages_passed) {
        const double & g = pow(prob_per_stage, CICIR_stages_passed);
        output_prob +=  exp((0.5-output_variate_value)/g) / g;
    }
    return norm * output_prob;
}
//map<int,double> FullEmccdHistogramFitFCN::get_CICIR_histogram(
//    const double gain,
//    const double N_CICIR_pix,
//    const int N_EM_stages,
//    const int max_val_to_generate) const{
//
//    int min_N_stages=1;
//    double prob_per_stage = pow(gain, 1.0/N_EM_stages);
//
//
//    map<int,double> gain_for_N_stages;
//    for (int stages_passed = 1; stages_passed<=N_EM_stages; ++stages_passed){
//        gain_for_N_stages[stages_passed] = pow(prob_per_stage, stages_passed);
//    }
//
////    cerr<<"gain at 604:"<< gain_for_N_stages[604]<<endl;
////    const double frequency_normalisation =  N_CICIR_pix ;
//    const double frequency_normalisation = N_CICIR_pix / (N_EM_stages - min_N_stages);
//
//    map<int,double> CICIR_histogram;
//    for (int register_output_val=1; register_output_val<=max_val_to_generate; ++register_output_val){
//        double output_frequency=0.0;
//        for (int CICIR_stages_passed = min_N_stages; CICIR_stages_passed<=N_EM_stages; ++CICIR_stages_passed){
//            const double & g = gain_for_N_stages[CICIR_stages_passed];
//
//            output_frequency += frequency_normalisation * exp((0.5-register_output_val)/g) / g;
//        }
//        CICIR_histogram[register_output_val]=output_frequency;
//    }
//    return CICIR_histogram;
//}

//=======================================================================================
//EMCCD_fixed_low_flux_gain_variate::EMCCD_fixed_low_flux_gain_variate()
//:gen(NULL), init_state(false){}

FixedLowFluxGainVariate::FixedLowFluxGainVariate(int n_photons, double gain,
        unuran::StreamWrapper& rns)
    :gen(NULL)
{

    assert(n_photons<FixedLowFluxGainVariate::n_photons_limit_for_low_flux);

    UNUR_DISTR *distr;    /* distribution object                   */
    UNUR_PAR   *par;      /* parameter object                      */

    vector<double> prob_vec;
    prob_vec.resize(FixedLowFluxGainVariate::max_variate_value, 0.0);

    for (size_t i=1; i!=prob_vec.size(); ++i) {
        prob_vec[i] = low_flux_gain_pmf(n_photons, gain, i);
    }

    distr = unur_distr_discr_new();

    unur_distr_discr_set_pv(distr, &prob_vec[0] , prob_vec.size());


//    par = unur_auto_new(distr);
    par = unur_dgt_new(distr);
//    par = unur_dau_new(distr);

    unur_set_urng(par, rns.get_stream_pointer());
    gen = unur_init(par);
    unur_distr_free(distr);

    if (gen==NULL) {
        throw std::runtime_error(
                "FixedLowFluxGainVariate - could not construct generator.");
    }
}

FixedLowFluxGainVariate::~FixedLowFluxGainVariate()
{
    unur_free(gen);
}


int FixedLowFluxGainVariate::operator()()
{
    return unur_sample_discr(gen);
}

//=======================================================================================
CicirVariate::CicirVariate(const int num_serial_register_stages,
        const double gain,
        unuran::StreamWrapper& rns)
    :gen(NULL)
{
    UNUR_DISTR *distr;    /* distribution object                   */
    UNUR_PAR   *par;      /* parameter object                      */

    vector<double> prob_vec;
    prob_vec.resize(FixedLowFluxGainVariate::max_variate_value, 0.0);

    for (size_t i=1; i!=prob_vec.size(); ++i) {
        prob_vec[i] = cicir_pmf(num_serial_register_stages,  gain, i);
    }

    distr = unur_distr_discr_new();

    unur_distr_discr_set_pv(distr, &prob_vec[0] , prob_vec.size());


//    par = unur_auto_new(distr);
    par = unur_dgt_new(distr);
//    par = unur_dau_new(distr);

    unur_set_urng(par, rns.get_stream_pointer());
    gen = unur_init(par);
    unur_distr_free(distr);

    if (gen==NULL) { throw std::runtime_error("EMCCD_low_flux_gain_variate::EMCCD_low_flux_gain_variate - could not construct generator."); }
}

CicirVariate::~CicirVariate()
{
    unur_free(gen);
}


int CicirVariate::operator()()
{
    return unur_sample_discr(gen);
}

//=======================================================================================
GainVariate::GainVariate(double gain_init,
                                     unuran::StreamWrapper& rns_init):
    fixed_flux_generators(FixedLowFluxGainVariate::n_photons_limit_for_low_flux, NULL),
    rns(&rns_init),
    gain(gain_init)
{}

int GainVariate::operator()(int n_photons_input)
{
    //Low flux case:
    if (n_photons_input<(int)fixed_flux_generators.size()) {
        if (n_photons_input) {
            //non-zero case
            if (!fixed_flux_generators[n_photons_input]) {
                //if NULL
                fixed_flux_generators[n_photons_input] =
                    new FixedLowFluxGainVariate(n_photons_input, gain, *rns);
            }

            return (fixed_flux_generators[n_photons_input])->operator()();
        } else { return 0; } //if zero photons input.
    }
    //high flux case
    unuran::GaussianRandomVariate high_flux_grv(
        n_photons_input*gain,
        sqrt(n_photons_input*gain*gain),
        *rns
    );
    return rint(high_flux_grv());
}

GainVariate::~GainVariate()
{
    for (size_t i=0; i!=fixed_flux_generators.size(); ++i) { delete(fixed_flux_generators[i]); }
}

//=======================================================================================
} //end namespace coela::emccd_detail
//=======================================================================================

EmccdModelParams::EmccdModelParams(const std::string& filename) {
        using namespace boost_serialization;
        (*this) = load_serializable_object_from_file<EmccdModelParams>(filename);
    }

void EmccdModelParams::write_to_file(const std::string& filename) {
        using namespace boost_serialization;
        save_serializable_object_to_file(filename, *this);
    }

std::ostream& operator<<(std::ostream& os, const EmccdModelParams& pars)
{
    os<<"#EMCCD model characteristics file\n";
    os<<"#Version 1.0\n";
    os<<"Bias pedestal: "<<pars.bias_pedestal<<"\n";
    os<<"Readout sigma: "<<pars.readout_sigma<<"\n";
    os<<"Photon gain: "<<pars.photon_gain<<"\n";
    os<<"CIC rate: "<<pars.serial_CIC_rate<<"\n";
    os<<"N EM register stages: "<<pars.N_EM_serial_register_stages<<"\n";
//    os<<"CIC gain: "<<pars.serial_CIC_gain<<"\n";
//    os<<"Background flux: "<<pars.background_flux_per_Pixel<<"\n";
    return os;
}

//EMCCD_model_parameters::EMCCD_model_parameters(const std::string& filename){
//    *this = serialization::load_serializable_object_from_file<EMCCD_model_parameters>(filename);
//}

//void EMCCD_model_parameters::write_to_file(const std::string& filename){
//
//}

EmccdModel::EmccdModel(const EmccdModelParams& parameter_init,
                         unuran::StreamWrapper& rns_in):
    params(parameter_init),
    photon_gain_variate(parameter_init.photon_gain, rns_in),
    readout_variate(parameter_init.bias_pedestal, parameter_init.readout_sigma, rns_in),
    CICIR_variate_(parameter_init.N_EM_serial_register_stages, parameter_init.photon_gain,
                   rns_in)
{}

int EmccdModel::readout_noise_variate()
{
    return nearbyint(readout_variate());
}

//=======================================================================================

}//end namespace coela
