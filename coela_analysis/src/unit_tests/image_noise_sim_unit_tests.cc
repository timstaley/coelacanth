/*
 * File:   FitsHeader_tests.cc
 * Author: ts337
 *
 * Created on 07 February 2011, 12:41
 */
#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include "../image_noise_sim.h"
#include "coela_utility/src/histogram_container.h"
//#include <tbb/tick_count.h>

#include <boost/filesystem.hpp>
#include <iostream>

using namespace coela;
using namespace std;

//SLOW!
///*
SUITE(Image_Noise_Sim)
{
    string test_suite_output_dir = "image_noise_sim_tests/";

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"image_noise_sim_tests\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(Init_Unuran_Seed) {
        //Need to do this within a function, hence declaration within a TEST macro.
        std::vector < unsigned long > seed;

        seed.push_back(111);
        seed.push_back(222);
        seed.push_back(333);
        seed.push_back(444);
        seed.push_back(555);
        seed.push_back(666);

        //NB changing the seed may cause some tests to fail due to overly tight / crude statistical checks.

        if (! unuran::
                StreamWrapper::unuran_package_has_been_seeded()) {

            unuran::
            StreamWrapper::set_unuran_package_seed(seed);
        }


    }

    TEST(uniform_random_pixel_variate) {
        PixelArray2d<int> img(5, 5, 0);


        unuran::StreamWrapper rns;

        image_noise_simulation::uniform_random_pixel_variate pixel_picker(
            img.range(), rns);

        size_t n_picks=1e5;
        for (size_t i=0; i!=n_picks; ++i) {
            img(pixel_picker())+=1;
        }

        double prob = 1.0/(double)img.range().n_pix();
        double expected_val = (double)n_picks * prob;
        double expected_SD = sqrt(expected_val*(1.0 - prob));
        for (PixelIterator pix(img.range()); pix!=pix.end; ++pix) {
            CHECK_CLOSE(expected_val, img(pix), expected_SD*3);
        }

//        img.write_to_file(test_suite_output_dir+"uniformly_picked_pixels.fits");

    }
//

    TEST(intensity_map_random_pixel_variate_setup_time) {

    }

    TEST(intensity_map_random_pixel_variate_pick_time) {

    }

    TEST(Photon_Shot_Noise) { //Slow!
        PixelArray2d<double> intensity_map(50,50, 0.5);

        intensity_map(2,2)=3.14;
        intensity_map(3,3)=15.555;
        intensity_map(4,4)=5.14;

        vector< HistogramContainer14bit > hists;

        unuran::StreamWrapper rns;

        size_t n_its =2e3;

        PixelArray2d<double> avg_img(intensity_map);
        avg_img.assign(0);

//        tbb::tick_count t0;
//        t0 = tbb::tick_count::now();

        for (size_t i=0; i!=n_its; ++i) {
            PixelArray2d<int> photon_img =
                image_noise_simulation::generate_photon_arrival_map(
                    intensity_map,
                    intensity_map.sum(),
                    0.0,
                    rns);
            avg_img +=photon_img;
        }

//        tbb::tick_count t1 = tbb::tick_count::now();
//        double time =(t1-t0).seconds();
//
//        cout<<"Performed photon_shot_noise sim in "<<time<<" seconds "<<endl;

//        double expected_phots_per_it = intensity_map.sum();
//        cout <<time / n_its / expected_phots_per_it *1e9<< " nanoseconds per photon" <<endl;
//        cout<<"( Expected phots per frame: " << expected_phots_per_it<< " )" << endl;


        intensity_map.write_to_file(test_suite_output_dir+"photon_shot_test_intensity_map.fits");
        image_noise_simulation::generate_photon_arrival_map(
            intensity_map,
            intensity_map.sum(),
            0.0,
            rns
        ).write_to_file(
            test_suite_output_dir+"photon_shot_test_sample_instance.fits");

        avg_img/=n_its;
        avg_img.write_to_file(test_suite_output_dir+"photon_shot_test_avg_map.fits");

        for (PixelIterator p_it(avg_img.range()); p_it!=p_it.end; ++p_it) {
            CHECK_CLOSE(intensity_map(p_it), avg_img(p_it), sqrt(intensity_map(p_it))*0.5);
        }
    }

///*

//    TEST(EMCCD_image_sim_dark_image_zero_cic){
//        image<int> photon_map(5,5, 0);
//
//        CCD_models::EMCCD_model_parameters model_pars;
//        model_pars.photon_gain=100;
//        model_pars.serial_CIC_gain=10;
//        model_pars.CIC_rate_per_pix=0.0;
//        model_pars.readout_sigma=10.5;
//        model_pars.bias_pedestal=500;
//
//
////        photon_map(2,2)=1;
////        photon_map(3,3)=2;
////        photon_map(4,4)=3;
//
//        HistogramContainer14bit  hist;
//
//                std::vector < unsigned long > seed;
//
//        seed.push_back(111);
//        seed.push_back(222);
//        seed.push_back(333);
//        seed.push_back(444);
//        seed.push_back(555);
//        seed.push_back(666);
//
//        unuran::
//        StreamWrapper::set_unuran_package_seed(seed, false);
//
//
//
//
//        unuran::StreamWrapper rns;
//
//        CCD_models::EMCCD_model model(model_pars, rns);
//
//
//        size_t n_its =1e5;
//
//        image<double> avg_img(photon_map.range().x_dim(), photon_map.range().y_dim(), 0);
//
//        for (size_t i=0; i!=n_its;++i){
//            image<int> dark_frame =
//                image_noise_simulation::EMCCD_simulated_image(
//                    photon_map, model, rns);
//            avg_img +=dark_frame;
//
//            for (PixelIterator pix(dark_frame.range()); pix!=pix.end; ++pix){
//                hist.count_int_value(dark_frame(pix));
//            }
//        }
//
//        avg_img/=(double)n_its;
//
//        for (PixelIterator pix(avg_img.range()); pix!=pix.end; ++pix){
//            CHECK_CLOSE( model_pars.bias_pedestal, avg_img(pix), 3 * model_pars.readout_sigma / sqrt(n_its) );
//        }
//
//        avg_img.write_to_file(test_suite_output_dir+"EMCCD_dark_run_sim_CIC_0.fits");
//
//        hist.write_to_file(test_suite_output_dir+"EMCCD_dark_run_sim_hist_CIC_0.txt");
//
//
//    }



//    TEST(EMCCD_image_sim_dark_image_nonzero_CIC){
//        image<int> photon_map(5,5, 0);
//
//        CCD_models::EMCCD_model_parameters model_pars;
//        model_pars.photon_gain=100;
//        model_pars.serial_CIC_gain=10;
//        model_pars.CIC_rate_per_pix=0.04;
//        model_pars.readout_sigma=10.5;
//        model_pars.bias_pedestal=500;
//
//
////        photon_map(2,2)=1;
////        photon_map(3,3)=2;
////        photon_map(4,4)=3;
//
//        HistogramContainer14bit  hist;
//
//                std::vector < unsigned long > seed;
//
//        seed.push_back(111);
//        seed.push_back(222);
//        seed.push_back(333);
//        seed.push_back(444);
//        seed.push_back(555);
//        seed.push_back(666);
//
//        unuran::
//        StreamWrapper::set_unuran_package_seed(seed, false);
//
//
//
//
//        unuran::StreamWrapper rns;
//
//        CCD_models::EMCCD_model model(model_pars, rns);
//
//
//        size_t n_its =1e5;
//
//        image<double> avg_img(photon_map.range().x_dim(), photon_map.range().y_dim(), 0);
//
//        for (size_t i=0; i!=n_its;++i){
//            image<int> dark_frame =
//                image_noise_simulation::EMCCD_simulated_image(
//                    photon_map, model, rns);
//            avg_img +=dark_frame;
//
//            for (PixelIterator pix(dark_frame.range()); pix!=pix.end; ++pix){
//                hist.count_int_value(dark_frame(pix));
//            }
//        }
//
//        avg_img/=(double)n_its;
//
//        for (PixelIterator pix(avg_img.range()); pix!=pix.end; ++pix){
//            CHECK_CLOSE( model_pars.bias_pedestal + model_pars.CIC_rate_per_pix*model_pars.serial_CIC_gain,
//                    avg_img(pix), 4 * model_pars.readout_sigma / sqrt(n_its) );
//        }
//
//        avg_img.write_to_file(
//            test_suite_output_dir+"EMCCD_dark_run_sim_CIC_"+string_utils::ftoa(model_pars.CIC_rate_per_pix)+ ".fits");
//
//        hist.write_to_file(test_suite_output_dir+"EMCCD_dark_run_sim_hist_CIC_" +string_utils::ftoa(model_pars.CIC_rate_per_pix)+  ".txt");
//    }

    /*
        TEST(EMCCD_image_sim_uniform_photon_image_zero_CIC){
            int uniform_n_photons =1;
            image<int> photon_map(5,5, uniform_n_photons);

            CCD_models::EMCCD_model_parameters model_pars;
            model_pars.photon_gain=100;
            model_pars.serial_CIC_gain=10;
            model_pars.CIC_rate_per_pix=0.00;
            model_pars.readout_sigma=10.5;
            model_pars.bias_pedestal=500;


    //        photon_map(2,2)=1;
    //        photon_map(3,3)=2;
    //        photon_map(4,4)=3;

            HistogramContainer14bit  hist;

                    std::vector < unsigned long > seed;

            seed.push_back(111);
            seed.push_back(222);
            seed.push_back(333);
            seed.push_back(444);
            seed.push_back(555);
            seed.push_back(666);

            unuran::
            StreamWrapper::set_unuran_package_seed(seed, false);




            unuran::StreamWrapper rns;

            CCD_models::EMCCD_model model(model_pars, rns);


            size_t n_its =1e5;

            image<double> avg_img(photon_map.range().x_dim(), photon_map.range().y_dim(), 0);

            for (size_t i=0; i!=n_its;++i){
                image<int> dark_frame =
                    image_noise_simulation::EMCCD_simulated_image(
                        photon_map, model, rns);
                avg_img +=dark_frame;

                for (PixelIterator pix(dark_frame.range()); pix!=pix.end; ++pix){
                    hist.count_int_value(dark_frame(pix));
                }
            }

            avg_img/=(double)n_its;

            for (PixelIterator pix(avg_img.range()); pix!=pix.end; ++pix){
                CHECK_CLOSE( model_pars.bias_pedestal + uniform_n_photons*model_pars.photon_gain,
                        avg_img(pix),
                        3 * sqrt(model_pars.readout_sigma*model_pars.readout_sigma + uniform_n_photons*model_pars.photon_gain*model_pars.photon_gain )/ sqrt(n_its) );
            }

            avg_img.write_to_file(
                test_suite_output_dir+"EMCCD_uniform_photon_image_zero_CIC.fits");

            hist.write_to_file(test_suite_output_dir+"EMCCD_uniform_photon_image_zero_CIC_hist.txt");
        }

        TEST(EMCCD_image_sim_varying_photon_image_nonzero_CIC){

            image<int> photon_map(5,5, 0.0);

            CCD_models::EMCCD_model_parameters model_pars;
            model_pars.photon_gain=100;
            model_pars.serial_CIC_gain=10;
            model_pars.CIC_rate_per_pix=0.04;
            model_pars.readout_sigma=10.5;
            model_pars.bias_pedestal=500;


            photon_map(2,2)=1;
            photon_map(3,3)=2;
            photon_map(4,4)=3;

            HistogramContainer14bit  hist;

                    std::vector < unsigned long > seed;

            seed.push_back(111);
            seed.push_back(222);
            seed.push_back(333);
            seed.push_back(444);
            seed.push_back(555);
            seed.push_back(666);

            unuran::
            StreamWrapper::set_unuran_package_seed(seed, false);




            unuran::StreamWrapper rns;

            CCD_models::EMCCD_model model(model_pars, rns);


            size_t n_its =1e5;

            image<double> avg_img(photon_map.range().x_dim(), photon_map.range().y_dim(), 0);

            for (size_t i=0; i!=n_its;++i){
                image<int> dark_frame =
                    image_noise_simulation::EMCCD_simulated_image(
                        photon_map, model, rns);
                avg_img +=dark_frame;

                for (PixelIterator pix(dark_frame.range()); pix!=pix.end; ++pix){
                    hist.count_int_value(dark_frame(pix));
                }
            }

            avg_img/=(double)n_its;

            avg_img.write_to_file(
                test_suite_output_dir+"EMCCD_image_sim_varying_photon_image_nonzero_CIC.fits");

            hist.write_to_file(test_suite_output_dir+"EMCCD_image_sim_varying_photon_image_nonzero_CIC_hist.txt");
        }

        TEST(EMCCD_image_sim_stochastic_photon_image_nonzero_CIC){
            double uniform_prob_photons =0.1;
            image<float> intensity_map(5,5, uniform_prob_photons);

            CCD_models::EMCCD_model_parameters model_pars;
            model_pars.photon_gain=100;
            model_pars.serial_CIC_gain=10;
            model_pars.CIC_rate_per_pix=0.1;
            model_pars.readout_sigma=10.5;
            model_pars.bias_pedestal=500;


    //        photon_map(2,2)=1;
    //        photon_map(3,3)=2;
    //        photon_map(4,4)=3;

            HistogramContainer14bit  hist;

                    std::vector < unsigned long > seed;

            seed.push_back(111);
            seed.push_back(222);
            seed.push_back(333);
            seed.push_back(444);
            seed.push_back(555);
            seed.push_back(666);

            unuran::
            StreamWrapper::set_unuran_package_seed(seed, false);




            unuran::StreamWrapper rns;

            CCD_models::EMCCD_model model(model_pars, rns);


            size_t n_its =1e5;

            image<double> avg_img(intensity_map.range().x_dim(), intensity_map.range().y_dim(), 0);

            for (size_t i=0; i!=n_its;++i){
                image<int> photon_map = image_noise_simulation::photon_arrival_map(intensity_map, 1.0, rns);
                image<int> gain_cal_frame =
                    image_noise_simulation::EMCCD_simulated_image(
                        photon_map, model, rns);
                avg_img +=gain_cal_frame;

                for (PixelIterator pix(gain_cal_frame.range()); pix!=pix.end; ++pix){
                    hist.count_int_value(gain_cal_frame(pix));
                }
            }

            avg_img/=(double)n_its;

        avg_img.write_to_file(
                test_suite_output_dir+ string(UnitTest::CurrentTest::Details()->testName) +".fits");


            hist.write_to_file(test_suite_output_dir+string(UnitTest::CurrentTest::Details()->testName)+"_hist.txt");
        }


    */
}

