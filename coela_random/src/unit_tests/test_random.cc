#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include "../random.h"

#include <boost/filesystem.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <vector>

#include<unuran.h>

#include <iostream>
#include <fstream>

using namespace coela;
using namespace std;

SUITE(random)
{

    string test_suite_output_dir = "random_tests/";

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"random number gen\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(Unuran_Seed) {
        //Need to do this within a function, hence declaration within a TEST macro.
        std::vector < unsigned long > seed;

        seed.push_back(1111);
        seed.push_back(222);
        seed.push_back(333);
        seed.push_back(444);
        seed.push_back(555);
        seed.push_back(666);

        unuran::
        StreamWrapper::set_unuran_package_seed(seed, false);

    }

    TEST(Unuran_Seed_Manipulation) {
        using namespace coela::unuran;

        std::vector < unsigned long > seed;
        seed.push_back(1111);
        seed.push_back(222);
        seed.push_back(333);
        seed.push_back(444);
        seed.push_back(555);
        seed.push_back(666);

        int n_advancements = 3;

        StreamWrapper::set_unuran_package_seed(seed, false);
        StreamWrapper rns1;
        StreamWrapper::advance_package_seed(n_advancements);
        StreamWrapper alternate_rns1;

        UniformRandomVariate urv1(0,1.0,rns1);
        UniformRandomVariate alternate_urv1(0,1.0,alternate_rns1);


        StreamWrapper::set_unuran_package_seed(seed, false);
        StreamWrapper rns2;
        StreamWrapper::advance_package_seed(n_advancements);
        StreamWrapper alternate_rns2;

        UniformRandomVariate urv2(0,1.0,rns2);
        UniformRandomVariate alternate_urv2(0,1.0,alternate_rns2);

        vector<double> results1, results2;

        for (size_t i=0; i!=10; ++i) {
            results1.push_back(urv1());
            results2.push_back(urv2());
        }

        CHECK(results1==results2);

        vector<double> alt_results1, alt_results2;

        for (size_t i=0; i!=10; ++i) {
            alt_results1.push_back(alternate_urv1());
            alt_results2.push_back(alternate_urv2());
        }

        CHECK(alt_results1==alt_results2);

        CHECK(alt_results1!=results1); //(sanity check)
        CHECK(alt_results2!=results2);
    }

//
//    TEST(Unuran_Uniform_Distribution) {
//        using namespace coela::unuran;
//
//        StreamWrapper rns;
//
//        uniform_random_variate rv_0_1(0, 1.0, rns);
//
//        coela::HistogramContainer14bit hist;
//
//        size_t n_its = 1e6;
//
//        int n_bins=10;
//
//        //visual check
////        for (size_t i=0; i!=10;++i){
////            cout<<"rv_0_1: " <<rv_0_1()<<endl;
////        }
//
//        for (size_t i=0; i!=n_its; ++i) {
//            hist.count_int_value(floor(rv_0_1()*n_bins));
//        }
//
//        int n_per_bin = n_its / n_bins;
//
//        for (int ind=0; ind<n_bins; ++ind) {
//            CHECK_CLOSE(n_per_bin, hist[ind], n_per_bin/100);
//        }
//
//        hist.write_to_file(test_suite_output_dir+"uniform_dist_tenths.txt");
//    }
//


//
//    TEST(Boost_Poisson_Distribution_mean_1){
//        double mean =1.0;
////        coela::simple_01_uniform_random_variate rv(172);
//
//
//        unuran_poisson_random_variate poisson_gen(mean);
//
//
//        coela::HistogramContainer14bit hist;
//
//        size_t n_its = 1e6;
//
//        boost::math::poisson_distribution<> poisson_pdf(mean);
//
//        for (size_t i=0; i!=n_its; ++i){
//            hist.count_int_value((int)poisson_gen());
//
//        }
//
//        for (int ind=0; ind<100; ++ind){
//            double predicted_num = boost::math::pdf(poisson_pdf,ind)*n_its;
//            CHECK_CLOSE(predicted_num, (double)hist[ind] , max(predicted_num/10, 5.0));
//        }
//
//        hist.write_to_file(test_suite_output_dir+"poisson_dist_mean"+string_utils::ftoa(mean,2)+".txt");
//
//    }


//    TEST(Boost_Poisson_Distribution_mean_11_and_a_quarter){
//        double mean =11.25;
//        coela::simple_01_uniform_random_variate rv(123);
//
//
//        boost::poisson_distribution<> poisson_gen(mean);
//
//
//        coela::HistogramContainer14bit hist;
//
//        size_t n_its = 1e6;
//
//        boost::math::poisson_distribution<> poisson_pdf(mean);
//
//        int zeroes=0;
//        for (size_t i=0; i!=n_its; ++i){
//            assert(fmod(poisson_gen(rv),1.0) == 0.0);
//            int result = poisson_gen(rv);
//            if (result==0) ++zeroes;
//            hist.count_int_value(result);
//
//        }
//
//        cerr<<"zeros "<<zeroes<<endl;
//
//        for (int ind=0; ind<100; ++ind){
//            double predicted_num = boost::math::pdf(poisson_pdf,ind)*n_its;
//            CHECK_CLOSE(predicted_num, (double)hist[ind] , max(predicted_num/10, 10.0));
//        }
//
//        hist.write_to_file(test_suite_output_dir+"boost_poisson_dist_mean"+string_utils::ftoa(mean,4)+".txt");
//
//    }


//    TEST(unuran_poisson_11_quarter) {
//
//        std::vector < unsigned long > seed;
//        seed.push_back(1111);
//        seed.push_back(222);
//        seed.push_back(333);
//        seed.push_back(444);
//        seed.push_back(555);
//        seed.push_back(666);
//
//        unuran::
//        StreamWrapper::set_unuran_package_seed(seed, false);
//
//        using namespace random;
//        double mean =11.25;
//
//        unuran::StreamWrapper rns1;
//
//        unuran::PoissonRandomVariate gen(mean, rns1);
//        coela::HistogramContainer14bit hist;
//
//        size_t n_its = 1e6;
//
//        int zeroes=0;
//        for (size_t i=0; i!=n_its; ++i) {
//            assert(fmod(gen(),1.0) == 0.0);
//            int result = gen();
//            if (result==0) { ++zeroes; }
//            hist.count_int_value(result);
//
//        }
//
////        cerr<<"zeros "<<zeroes<<endl; //turns out the boost method produced far too many zero counts
//
//        boost::math::poisson_distribution<> poisson_pdf(mean);
//
//        for (int ind=0; ind<100; ++ind) {
//            double predicted_num = boost::math::pdf(poisson_pdf,ind)*n_its;
//            CHECK_CLOSE(predicted_num, (double)hist[ind] , max(predicted_num/10, 10.0));
//        }
//
//        hist.write_to_file(test_suite_output_dir+"unuran_poisson_dist_mean"+string_utils::ftoa(
//                               mean,4)+".txt");
//
//    }

    TEST(unuran_poisson_seed_varies) {
        //Confirm that different poisson variates using the same stream will create different outputs:


        double mean = 23.1234;

        unuran::StreamWrapper rns1, rns2;

        //NB changing the seed may cause some tests to fail due to overly tight / crude statistical checks.

        //Here we can explicitly use different streams, or use default behaviour
        //(which  is to feed both from the default random number stream)

        //The advantage to explicitly specifying different streams? MULTI-THREADING.
        //(Each stream object looks after its own state)
        //But also, using a specific stream object makes the user aware of the need to supply a seed
        //(and possibly vary said seed between simulations)

        unuran::PoissonRandomVariate p_gen1(mean, rns1);
        unuran::PoissonRandomVariate p_gen2(mean, rns1);
//        unuran_poisson_random_variate p_gen1(mean);
//        unuran_poisson_random_variate p_gen2(mean);


        vector<int> p1, p2;

        for (size_t i=0; i!=10; ++i) {
            p1.push_back(p_gen1());
            p2.push_back(p_gen2());

            //visual check:
//            cout<<i<<" - rns poisson 1: "<<gen1()<<"; rns poisson 2: "<<gen2()<<endl;
        }

        CHECK(p1!=p2);
    }

    TEST(unuran_gaussian) {

        unuran::StreamWrapper rns;

        double mean(5.0), sigma(2.5);
        unuran::GaussianRandomVariate grv(mean, sigma ,rns);

        size_t n_its=1e5;
        vector<double> vals; vals.reserve(n_its);

        ofstream datfile(string(test_suite_output_dir+"gaussian_vals.txt").c_str());

        for (size_t i=0; i!=n_its; ++i) {
            vals.push_back(grv());
            datfile<<vals.back()<<"\n";
        }

        datfile.close();

//        cout<<"gaussian limits: "<< endl;
//        cout<<"mean: "<< 2* sigma/sqrt((double)n_its)<<endl;
//        cout<<"sigma: "<< 2* sigma*sigma/sqrt((double)n_its) <<endl;

//        CHECK_CLOSE(mean, vector_mean(vals),  2*sigma/sqrt((double)n_its));
//        CHECK_CLOSE(sigma, vector_std_dev(vals),  sigma*sigma/sqrt((double)n_its));

        //Now test the parameter update: (produces weird warning messages, stuff it.)
//        vals.clear();
//        mean = 55;
//        sigma = 15;
//
//        grv.update_params(mean,sigma);
//        for (size_t i=0; i!=n_its;++i){
//            vals.push_back( grv());
//        }
//        CHECK_CLOSE(mean, vector_mean(vals),  sigma/sqrt((double)n_its));
//        CHECK_CLOSE(sigma, vector_std_dev(vals),  sigma*sigma/sqrt((double)n_its));
    }

}
