/*
 * File:   random.h
 * Author: ts337
 *
 * Created on 28 February 2011, 11:58
 */

#ifndef COELA_RANDOM_H
#define COELA_RANDOM_H


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/noncopyable.hpp>
#include <vector>

#include <unuran.h>

namespace coela {
//=========================================================================================================
namespace boost_random{

///Convenient wrapper class to boost URV.
///Returns a generator class that returns operator() in the same way.
class UniformRandomVariate {
public:
    UniformRandomVariate(unsigned int seed);
    double operator()() {return rv();} //returns a random double between 0 and 1

private:
    boost::mt19937 rng;
    boost::uniform_real<> dist;
    boost::variate_generator< boost::mt19937, boost::uniform_real<> > rv;
};

} //end namespace coela::boost_wrappers

//=========================================================================================================
//=========================================================================================================
namespace unuran {
//-----------------------------------------------------------------------
//Alternative:
// Based around unuran since bug found in boost poisson variate generator (too many zeroes):
class StreamWrapper: private boost::noncopyable_::noncopyable {
public:
    StreamWrapper();

    ~StreamWrapper();

    static void set_unuran_package_seed(
        std::vector<unsigned long>& six_seed_numbers,
        bool check_if_already_set=true);

    static void advance_package_seed(const size_t n_steps);

    static bool unuran_package_has_been_seeded() {return package_seed_set;}

    UNUR_URNG * get_stream_pointer() {return rn_stream_ob_ptr;}

//    void reset()

private:
    static int num_streams_created;
    static bool package_seed_set;
    UNUR_URNG * rn_stream_ob_ptr;

};

//=========================================================================================================

class UniformRandomVariate {
public:
    /// Init with default RNG stream -- deprecated to encourage good use of seeds
//    unuran_poisson_random_variate(double lambda);

    /// Init with explicit RNG stream
    ///NB the stream must not be destructed while variate is in use!
    UniformRandomVariate(double lower_bound, double upper_bound,
                           StreamWrapper&);

    double operator()();
    ~UniformRandomVariate();
private:
//    double lower, upper;
    UNUR_GEN *gen;
};

//=========================================================================================================

class GaussianRandomVariate {
public:

    /// Init with explicit RNG stream
    ///NB the stream must not be destructed while variate is in use!
    GaussianRandomVariate(double mean, double sigma, StreamWrapper&);


    //Deprecated -- seems to work but produces weird warning messages from UNURAN library.
//    void update_params(double mean, double sigma);

    double operator()();
    ~GaussianRandomVariate();
private:
    UNUR_GEN *gen;
};

//=========================================================================================================
class PoissonRandomVariate {
public:
    /// Init with default RNG stream -- deprecated to encourage good use of seeds
//    unuran_poisson_random_variate(double lambda);

    /// Init with explicit RNG stream
    ///NB the stream must not be destructed while poisson variate is in use!
    PoissonRandomVariate(double lambda, StreamWrapper&);

    int operator()();
    ~PoissonRandomVariate();
private:
    UNUR_GEN *gen;
};

//=========================================================================================================
}//end namespace coela::unuran
}//end namespace coela


#endif  /* RANDOM_H */

