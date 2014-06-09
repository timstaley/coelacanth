
#include "../gain_utils.h"

#include <UnitTest++/UnitTest++.h>
#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>

#include <cmath>
using namespace coela;
using namespace std;
/*
SUITE(gain_utils){
    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup){
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(histogram_convolution){
        map<int, double> single_point_hist;
        double first_point_count=10000;

        double sigma = 10.0;
        int peak_val = 100;
        for (int i = -10*sigma+peak_val; i< 10*sigma+peak_val; ++i){
            single_point_hist[i] = 0;
        }
        single_point_hist[peak_val] = first_point_count;

        map<int, double> simple_conv =
            gain_utils::convolve_histogram_with_gaussian(single_point_hist, 10.0, 50);


        //The convolved output is only valid in regions where we can apply a full kernel
        // i.e. a subset of the original histogram
        int original_min = gain_utils::get_min_value_pair(single_point_hist).first;

        int convolved_min = gain_utils::get_min_value_pair(simple_conv).first;
        CHECK(convolved_min > original_min);

        double convolved_counts = gain_utils::sum_counts_in_histogram(simple_conv);
        CHECK_CLOSE(first_point_count, convolved_counts, first_point_count*0.005);

        gain_utils::write_histogram_data_to_file(single_point_hist,
                test_suite_output_dir+"single_point_hist.txt");
        gain_utils::write_histogram_data_to_file(simple_conv,
                test_suite_output_dir+"single_point_conv_hist.txt");


        map<int, double> exp_hist;
        for (int x=101; x!=16000; ++x){
            exp_hist[x] = 0.05 * exp(-x/10.0);
        }

        map<int, double> exp_conv=
            gain_utils::convolve_histogram_with_gaussian(exp_hist, 10.0, 50);

        gain_utils::write_histogram_data_to_file(exp_hist,
                test_suite_output_dir+"exp_hist.txt");
        gain_utils::write_histogram_data_to_file(exp_conv,
                test_suite_output_dir+"exp_conv_hist.txt");

    }

    TEST(CICIR_histogram){
        map<int,long> blank_hist;
        for (int value=0; value!=1500; value++){
            blank_hist[value] = 0;
        }

        gain_utils::full_EMCCD_histogram_fit_FCN full_model(blank_hist);
        double gain= 100;
        int N_CICIR_pix =  500000;
        int max_val_to_generate = 12000;
        int N_EM_stages = 604;

        map<int, double> CICIR_hist =
            full_model.get_CICIR_model_histogram(gain, N_CICIR_pix ,
                                            N_EM_stages,
                                            max_val_to_generate);

        double hist_counts = gain_utils::sum_counts_in_histogram(CICIR_hist);
//        cerr<<"hist counts " << hist_counts / N_CICIR_pix <<endl;
        CHECK_CLOSE(N_CICIR_pix , hist_counts, N_CICIR_pix*0.01);

        double mean_CICIR_gain =0.0;
        for (map<int,double>::const_iterator it=CICIR_hist.begin();
                it!=CICIR_hist.end(); ++it){
            mean_CICIR_gain += it->first*it->second;
        }
        mean_CICIR_gain /=hist_counts;
        double expected_mean_gain = gain / log(gain);
        CHECK_CLOSE(expected_mean_gain, mean_CICIR_gain, expected_mean_gain*0.05);

//        cout<<" Mean gain: " <<mean_gain<<"; Expected mean gain: "<< expected_mean_gain<<endl;

    }

    TEST(calculate_thresholded_SNR_runs){
//        const double RO=10;
//        const double gain= 100;
//        const double light_level=0.1;
//        const double dark_current_level=0.05;
//        const double CICIR_rate=0.05;
//        const double threshold=0.3;
//
//        gain_utils::Thresholded_SNR_Calculator calc1(
//                RO,
//                gain,
//                light_level,
//                dark_current_level,
//                CICIR_rate
//                );
//
//        double SNR = calc1.calc_SNR_at_threshold(threshold);
    }
    TEST(thresholded_SNR_converges_to_ideal_case){

        const double RO=1;
        const double gain= 100;
        const double light_level=0.1;
        const double CICIR_rate=0.00;
        const double threshold=0.1;

        gain_utils::Thresholded_SNR_Calculator ideal_case(
                RO,
                gain,
                light_level,
                CICIR_rate
                );

        double SNR = ideal_case.calc_SNR_at_threshold(threshold);

        double ideal_SNR = sqrt(light_level);
//        cerr<<"ideal SNR" <<ideal_SNR<<endl;
//        cerr<<"calc SNR" <<SNR<<endl;

        CHECK(SNR < ideal_SNR);
        CHECK_CLOSE(ideal_SNR, SNR, ideal_SNR*0.1);
    }


//
//    TEST(readout_model_consistency){
//        gain_utils::gain_info model_pars;
//        model_pars.bias_pedestal = 25;
//        model_pars.readout_sigma = 12.5;
//        model_pars.N_dark_pix = 100000;
//
//        model_pars.N_light_pix=0;
//        model_pars.gain=100;
//        model_pars.N_serial_CIC_pix=0;
//
//        vector<double> model_pars_vec;
//        model_pars_vec.push_back(model_pars.bias_pedestal);
//        model_pars_vec.push_back(model_pars.N_dark_pix);
//        model_pars_vec.push_back(model_pars.readout_sigma);
//        model_pars_vec.push_back(model_pars.N_light_pix);
//        model_pars_vec.push_back(model_pars.gain);
//        model_pars_vec.push_back(model_pars.N_serial_CIC_pix);
//
//        map<int,long> blank_hist;
//        for (int value=0; value!=15000; value++){
//            blank_hist[value] = 0;
//        }
//        gain_utils::full_EMCCD_histogram_fit_FCN full_model(blank_hist);
//
//        map<int, double> model_hist = full_model.get_model_histogram(model_pars_vec);
//
//        map<int, long> approx_model_hist;
//        for (map<int,double>::const_iterator it = model_hist.begin();
//                it!=model_hist.end(); ++it){
//            approx_model_hist[ it->first ] = rint(it->second);
//        }
//
//
//        vector<double> gpars = gain_utils::fit_readout_gaussian(approx_model_hist);
//
//        CHECK_CLOSE(model_pars.bias_pedestal, gpars[0], 1.0);
//        CHECK_CLOSE(model_pars.readout_sigma, gpars[2], 1.0);
//
//        cerr<<"model / fit bias:" << model_pars.bias_pedestal <<" , "<< gpars[0]<<endl;
//        cerr<<"model / fit ROsig:" << model_pars.readout_sigma <<" , "<< gpars[2]<<endl;
//
//        //Now check the full fit works too:
//    }

}
 */