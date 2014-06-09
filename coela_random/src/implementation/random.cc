/*
 * File:   random.cc
 * Author: ts337
 *
 * Created on 28 February 2011, 11:58
 */

#include "../random.h"
#include <boost/random/poisson_distribution.hpp>

#include <unuran.h>
#include <unuran_urng_rngstreams.h>
#include <vector>
#include <stdexcept>
#include <string>

using std::string;
using std::logic_error;
using std::to_string;

namespace coela {
//=========================================================================================================
namespace boost_random {
UniformRandomVariate::UniformRandomVariate(unsigned int seed):
    rng(seed), dist(0,1), rv(rng,dist) {}
} //end namespace coela::boost_random
//=========================================================================================================
//=========================================================================================================
namespace unuran {

int StreamWrapper::num_streams_created=0;
bool StreamWrapper::package_seed_set=false;

void StreamWrapper::set_unuran_package_seed(
    std::vector<unsigned long>& six_seed_numbers,
    bool check_if_already_set
)
{
    if (six_seed_numbers.size()!=6) {
        throw logic_error("random_number_stream::set_unuran_package_seed - "
                            "please supply vector of six numbers");
    }

    if (check_if_already_set && unuran_package_has_been_seeded()) {
        throw logic_error("random_number_stream::set_unuran_package_seed - "
                            "package seed already set");
    }

    RngStream_SetPackageSeed(& six_seed_numbers[0]);

    package_seed_set=true;

}

void StreamWrapper::advance_package_seed(const size_t n_steps)
{
    for (size_t i=0; i!=n_steps; ++i) {
        StreamWrapper temp_instantiation_to_force_seed_advance;
    }
    return;
}

StreamWrapper::StreamWrapper()
{
    if (!package_seed_set) throw logic_error(
            "random_number_stream::random_number_stream() - please set package seed first");

    rn_stream_ob_ptr = unur_urng_rngstream_new(string("urng-"+
            to_string(num_streams_created)).c_str());

    if (rn_stream_ob_ptr == NULL) throw std::runtime_error(
            "random_number_stream::random_number_stream - could not create new rng stream");
    num_streams_created++;
}

StreamWrapper::~StreamWrapper()
{
    unur_urng_free(rn_stream_ob_ptr);
}


//=========================================================================================================
UniformRandomVariate::UniformRandomVariate(double lower_bound, double upper_bound,
        StreamWrapper& rns)
    :gen(NULL)
{

    assert(lower_bound < upper_bound);
    UNUR_DISTR *distr;    /* distribution object                   */
    UNUR_PAR   *par;      /* parameter object                      */
    double urv_params[2];
    urv_params[0]=lower_bound;
    urv_params[1]=upper_bound;
    distr = unur_distr_uniform(urv_params, 2);
    par = unur_cstd_new(distr);
    unur_set_urng(par, rns.get_stream_pointer());
    gen = unur_init(par);
    unur_distr_free(distr);

    if (gen==NULL) { throw std::runtime_error("uniform_random_variate (constructor) - could not construct generator."); }

}
UniformRandomVariate::~UniformRandomVariate()
{
    unur_free(gen);
}


double UniformRandomVariate::operator()()
{
    return unur_sample_cont(gen);
}


//=========================================================================================================

GaussianRandomVariate::GaussianRandomVariate(double mean, double sigma,
        StreamWrapper& wrapped_rng)
    :gen(NULL)
{
    UNUR_DISTR *distr;    /* distribution object                   */
    UNUR_PAR   *par;      /* parameter object                      */
    double gaussian_params[2];
    gaussian_params[0]=mean;
    gaussian_params[1]=sigma;
    distr = unur_distr_normal(gaussian_params, 2);
    par = unur_cstd_new(distr);

    // Choose generation method: http://statmath.wu.ac.at/unuran/doc/unuran.html#normal
    //TO DO: Profile which is faster for likely use case of single generation per instance.
//    unur_cstd_set_variant(par, 2);


    unur_set_urng(par, wrapped_rng.get_stream_pointer());
    gen = unur_init(par);
    unur_distr_free(distr);

    if (gen==NULL) { throw std::runtime_error("gaussian_random_variate::gaussian_random_variate - could not construct generator."); }
}

GaussianRandomVariate::~GaussianRandomVariate()
{
    unur_free(gen);
}

//void gaussian_random_variate::update_params(double mean, double sigma){
//    double gaussian_params[2];
//    gaussian_params[0]=mean;
//    gaussian_params[1]=sigma;
//
//    unur_distr_cont_set_pdfparams(unur_get_distr(gen), gaussian_params, 2);
//    unur_distr_cont_upd_mode(unur_get_distr(gen));
//    unur_distr_cont_upd_pdfarea(unur_get_distr(gen));
//    if ( unur_reinit(gen) ) throw runtime_error("Gaussian_random_variate::update_params - threw an error code");
//}

double GaussianRandomVariate::operator()()
{
    return unur_sample_cont(gen);
}

//=========================================================================================================
//=========================================================================================================

// Use default random number stream- deprecated to encourage good use of seeds
//poisson_random_variate::poisson_random_variate(double lambda)
//:gen(NULL)
//{
//    UNUR_DISTR *distr;    /* distribution object                   */
//    UNUR_PAR   *par;      /* parameter object                      */
//    double poisson_params[1];
//    poisson_params[0]=lambda;
//    distr = unur_distr_poisson(poisson_params, 1);
//    par = unur_auto_new(distr);
//    gen = unur_init(par);
//    unur_distr_free(distr);
//
//    if (gen==NULL) throw std::runtime_error("simple_poisson_random_variate::simple_poisson_random_variate - could not construct generator.");
// }

PoissonRandomVariate::PoissonRandomVariate(double lambda,
        StreamWrapper& wrapped_rng)
    :gen(NULL)
{
    UNUR_DISTR *distr;    /* distribution object                   */
    UNUR_PAR   *par;      /* parameter object                      */
    double poisson_params[1];
    poisson_params[0]=lambda;
    distr = unur_distr_poisson(poisson_params, 1);
    par = unur_dstd_new(distr);
    unur_set_urng(par, wrapped_rng.get_stream_pointer());
    gen = unur_init(par);
    unur_distr_free(distr);

    if (gen==NULL) {
        throw std::runtime_error(
            "simple_poisson_random_variate - could not construct generator.");
    }
}

PoissonRandomVariate::~PoissonRandomVariate()
{
    unur_free(gen);
}

int PoissonRandomVariate::operator()()
{
    return unur_sample_discr(gen);
}

//=========================================================================================================
}//end namespace coela::unuran
}//end namespace coela
