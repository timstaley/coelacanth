#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>
#include <string>
#include <boost/filesystem.hpp>
#include <UnitTest++/UnitTest++.h>


#include "../emccd_model.h"

using namespace coela;
using namespace std;

SUITE(CCD_model)
{
    using namespace emccd_detail;


    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }


    TEST(integer_factorial) {
        CHECK_EQUAL(1ul, integer_factorial(0));
        CHECK_EQUAL(1ul, integer_factorial(1));
        CHECK_EQUAL(2ul, integer_factorial(2));
        CHECK_EQUAL(6ul, integer_factorial(3));
        CHECK_EQUAL(24ul, integer_factorial(4));

//        cout<<"Unsigned long max: " << numeric_limits<unsigned long>::max() <<endl;
//
//        double max_int_fac =1;
//        while (approx_factorial_double(max_int_fac) < numeric_limits<unsigned long>::max())
//            max_int_fac++;
//
//        cout<<"Implies max input for integer_factorial:"<<max_int_fac<<endl;

    }


    TEST(EMCCD_pars_write_to_file) {
        EmccdModelParams pars;

        pars.bias_pedestal=1.1;
        pars.readout_sigma=2.1;
        pars.photon_gain=3.1;
        pars.serial_CIC_rate=4.1;
        pars.N_EM_serial_register_stages=604;

        pars.write_to_file(test_suite_output_dir+"sample_EMCCD_model_pars.xml");
        EmccdModelParams loaded(test_suite_output_dir+"sample_EMCCD_model_pars.xml");
        CHECK_EQUAL(pars.photon_gain, loaded.photon_gain);
    }

    TEST(EMCCD_gain_PMF) {
        double gain = 100;

        CHECK_CLOSE(exp(-1/gain)/gain , low_flux_gain_pmf(1, gain, 1) , 1e-6);

        for (int n_photons = 1; n_photons<=10; ++n_photons) {
            double sum_probability=0;
            for (int x_out=n_photons; x_out<10000; ++x_out) {
                double PMF_val = low_flux_gain_pmf(n_photons, gain, x_out);
                CHECK(PMF_val>0);
                CHECK(PMF_val<1);
                sum_probability+=PMF_val;
            }
            //Any more stringent tests (than 1 percent) fail due to numerical roundoff in sum.
            CHECK_CLOSE(1.0, sum_probability,
                        0.01);
        }

    }

//
//    TEST(Write_gain_PMFs_to_file_for_inspection){
//        double gain =100;
//        for (int n_photons = 1; n_photons<=10; ++n_photons){
//            ofstream datfile(string(test_suite_output_dir+"EMCCD_gain_PMF_for_"+string_utils::itoa(n_photons)+"_photons.txt").c_str());
//            for (int x_out=n_photons; x_out<10000; ++x_out){
//                double PMF_val = EMCCD_gain_PMF_at_value(n_photons, gain, x_out);
//                datfile<<x_out<<" "<<PMF_val<<"\n";
//            }
//            datfile.close();
//        }
//    }

    int hist_max = 10000;
    void count_result(vector<unsigned int>& hist, int result){
        if (fixed < 0) {
            throw logic_error(
                "emccd gain variate should not produce negative values");
        }
        if (result < hist_max) {
            ++hist[result];
        }
    }

    TEST(gain_variates) {
        vector<unsigned long> seed;
        seed.push_back(11);
        seed.push_back(25);
        seed.push_back(35);
        seed.push_back(44);
        seed.push_back(55);
        seed.push_back(66);
        unuran::StreamWrapper::set_unuran_package_seed(seed,
                false);

        unuran::StreamWrapper rns1;
        double gain=100.0;

        size_t n_its_per_hist = 1e6;

        GainVariate var_gen(gain, rns1);

        //Test for low flux:
        for (int n_photons=1;
                n_photons!=FixedLowFluxGainVariate::n_photons_limit_for_low_flux; ++n_photons) {
            vector<unsigned int> fixed_results(hist_max, 0);
            vector<unsigned int> var_results(hist_max, 0);

            FixedLowFluxGainVariate fixed_gen(n_photons, gain, rns1);

////            check default constructor followed by assignment:
//            EMCCD_fixed_low_flux_gain_variate gen;
//            gen = EMCCD_fixed_low_flux_gain_variate(n_photons, gain, rns1);

            for (size_t i=0; i!=n_its_per_hist; ++i) {
                count_result(fixed_results, fixed_gen());
                count_result(var_results, var_gen(n_photons));
            }


            if (n_photons==0) {
                CHECK_EQUAL(n_its_per_hist, fixed_results[0]);
                CHECK_EQUAL(n_its_per_hist, var_results[0]);
            }


            for (int ind=0; ind<100; ++ind) {
                double predicted_num = n_its_per_hist *
                                        low_flux_gain_pmf(n_photons, gain, ind);

                CHECK_CLOSE(predicted_num, (double)fixed_results[ind] ,
                        max(5*sqrt(predicted_num), 5.0));
                CHECK_CLOSE(predicted_num, (double)var_results[ind] ,
                                        max(5*sqrt(predicted_num), 5.0));
            }

        }

        //Test for high flux:

        n_its_per_hist = 1e5;
        int n_photons =55; //test high flux
        {
            vector<unsigned int> high_flux_results(hist_max, 0);
            for (size_t i=0; i!=n_its_per_hist; ++i) {
                count_result(high_flux_results, var_gen(n_photons));
            }

//            actual_hist.write_to_file(
//                    test_suite_output_dir+"EMCCD_gain_var_"+
//                        to_string(n_photons)+"_actual_hist.txt");
        }
    }





}//end of test suite
