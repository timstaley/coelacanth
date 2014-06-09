#include <UnitTest++/UnitTest++.h>
#include "../drizzle.h"



#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>

#include <cmath>
using namespace coela;
using namespace std;

SUITE(drizzle)
{

    //run tests at float precision
//    typedef PixelArray2d<float> PixelArray2d_of_real;
//    double precision_limit = 1e-3; //Float precision results in errors around the 1e-3 level
    // (OK for lucky imaging single exposure purposes, but not great)

    //run tests at double precision
    typedef PixelArray2d<double> PixelArray2d_of_real;
    double precision_limit = 1e-11; //Double precision is accurate to much better tolerances.

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }



    TEST(single_input_identity_drizzle) {
        PixelArray2d_of_real input(10,10, 3.14);
        PixelArray2d_of_real input_weights(10,10, 1.0);

        //Add a hot pixel, downweight it:
        input(5,5)=50;
        double downweighted_input_sum = input(5,5);
        input_weights(5,5)=0;


        PixelArray2d_of_real output(10,10,0.0);
        PixelArray2d_of_real output_weights(10,10,0.0);

        PixelShift null_translation(0,0);

        double drizzle_scale =1.0;
        double pixfrac = 1.0;

        drizzle::translate_and_drizzle_frame(
            input, input_weights,
            output, output_weights,
            null_translation,
            drizzle_scale,
            pixfrac);


        PixelArray2d_of_real drz_img = drizzle::unweight_drizzle_results(output,output_weights);

        CHECK_CLOSE(input.sum() - downweighted_input_sum, output.sum(), precision_limit);
        CHECK_CLOSE(input.sum() - downweighted_input_sum,
                    drz_img.sum() *drizzle_scale*drizzle_scale, precision_limit);
        CHECK_CLOSE(input_weights.sum(), output_weights.sum(), precision_limit);


        for (PixelIterator i(drz_img.range()); i!=i.end; ++i) {
            if (input_weights(i)!=0.0) { CHECK_CLOSE(input(i), drz_img(i), precision_limit); }
            if (input_weights(i)==0.0) { CHECK_CLOSE(0.0 , drz_img(i), precision_limit); }
        }
    }

    TEST(single_input_null_shift_resample_drizzle) {
        PixelArray2d_of_real input(10,10, 3.14);
        PixelArray2d_of_real input_weights(10,10, 1.0);

        double drizzle_scale =0.5;
        double pixfrac = 0.9;
        double zoom = 1.0 / drizzle_scale;

        PixelArray2d_of_real output(input.range().x_dim()*zoom ,
                                    input.range().y_dim()*zoom,
                                    0.0);

        PixelArray2d_of_real output_weights(output);

        PixelShift null_translation(0,0);

        drizzle::translate_and_drizzle_frame(
            input, input_weights,
            output, output_weights,
            null_translation,
            drizzle_scale,
            pixfrac);



        PixelArray2d_of_real drz_img = drizzle::unweight_drizzle_results(output,output_weights);

        CHECK_CLOSE(input.sum(), output.sum(), precision_limit);
        CHECK_CLOSE(input.sum(), drz_img.sum() *drizzle_scale*drizzle_scale, precision_limit);
        CHECK_CLOSE(input_weights.sum(), output_weights.sum(), precision_limit);

        for (PixelIterator i(input.range()); i!=i.end; ++i) {
            CHECK_CLOSE(input(i),  drz_img(i.x*2, i.y*2)   , precision_limit);
        }
        for (PixelIterator i(drz_img.range()); i!=i.end; ++i) {
            CHECK_CLOSE(3.14,  drz_img(i)   , precision_limit);
        }

    }

    TEST(drizzle_output_weighting_correct_for_non_int_shift) {

        PixelArray2d_of_real input(5,5, 0.0);
        PixelArray2d_of_real input_weights(input);
        input_weights.assign(1.0);

        input(3,3)=50;

        double drizzle_scale =1.0;
        double pixfrac = 0.6;
        double zoom = 1.0 / drizzle_scale;

        PixelShift input_translation(0.45,0.45);

        PixelArray2d_of_real output((input.range().x_dim() + ceil(input_translation.x)) *zoom ,
                                    (input.range().y_dim() + ceil(input_translation.y)) *zoom ,
                                    0.0);

        PixelArray2d_of_real output_weights(output);

        drizzle::translate_and_drizzle_frame(
            input, input_weights,
            output, output_weights,
            input_translation,
            drizzle_scale,
            pixfrac);

        PixelArray2d_of_real drz_img = drizzle::unweight_drizzle_results(output,output_weights);

        CHECK_CLOSE(input.sum(), output.sum(), precision_limit);
        CHECK_CLOSE(input.sum(), drz_img.sum() *drizzle_scale*drizzle_scale, precision_limit);
        CHECK_CLOSE(input_weights.sum(), output_weights.sum(), precision_limit);



    }


    TEST(generalized_drizzle_check_no_bounds_overflow) {
        PixelArray2d_of_real input(15,15, 0.0);
        PixelArray2d_of_real input_weights(input);
        input_weights.assign(1.0);

        input(5,5)=50;

        input_weights(2,2)=0;
        input_weights(7,7)=0.5;

        double drizzle_scale =0.5;
        double pixfrac = 1.0;
        double zoom = 1.0/drizzle_scale;

        PixelShift input_translation(2.0,-1.25);

        PixelArray2d_of_real output(input.range().x_dim()*zoom,
                                    input.range().y_dim()*zoom,
                                    0.0);
        PixelArray2d_of_real output_weights(output);

        drizzle::translate_and_drizzle_frame(
            input, input_weights,
            output, output_weights,
            input_translation,
            drizzle_scale,
            pixfrac);

        PixelArray2d_of_real drz_img = drizzle::unweight_drizzle_results(output,output_weights);

        CHECK_CLOSE(input.sum(), drz_img.sum()*drizzle_scale*drizzle_scale, precision_limit);

        input.write_to_file("drz_input.fits");
        input_weights.write_to_file("drz_input_weights.fits");
        output.write_to_file("drz_output_raw.fits");
        output_weights.write_to_file("drz_output_weights.fits");
        drz_img.write_to_file("drz_output_unweighted.fits");
    }

    TEST(dual_drizzle_test1) {
        PixelArray2d_of_real input(15,15, 3.14);
        PixelArray2d_of_real input_weights(input);
        input_weights.assign(1.0);

        input(5,5)=50;

        input_weights(2,2)=0;
        input_weights(7,7)=0.5;

        double drizzle_scale =0.5;
        double pixfrac = 1.0;
        double zoom = 1.0/drizzle_scale;

        PixelShift input_translation(2.0,-1.25);

        PixelArray2d_of_real output(input.range().x_dim()*zoom,
                                    input.range().y_dim()*zoom,
                                    0.0);
        PixelArray2d_of_real output_weights(output);

        PixelArray2d_of_real dual_output(output);

        drizzle::dual_translate_and_drizzle_frame( //test with same input to both outputs
            input, input, input_weights,
            output, dual_output, output_weights,
            input_translation,
            drizzle_scale,
            pixfrac);

        PixelArray2d_of_real drz_img1 = drizzle::unweight_drizzle_results(output,output_weights);
        PixelArray2d_of_real drz_img2 = drizzle::unweight_drizzle_results(dual_output,
                                        output_weights);

        CHECK(drz_img1==drz_img2);
    }

//    TEST(combine_drizzles){
//        PixelArray2d_of_real input1(15,15, 3.14);
//        PixelArray2d_of_real input1_weights(input1);
//        input1_weights.assign(1.0);
//
//        PixelArray2d_of_real input2(15,15, 0);
//        PixelArray2d_of_real input2_weights(input2);
//        input2_weights.assign(0.5);
//
//        input1(5,5)=50;
//
//        input2(6,8)=13;
//
//        input1_weights(2,2)=0;
//        input2_weights(7,7)=5.0;
//
//
//
//        double drizzle_scale =0.5;
//        double pixfrac = 0.9;
//        double zoom = 1.0/drizzle_scale;
//
//        pixel_shift input_translation1(2.0,-1.25), input_translation2(1.0,0.666);
//
//        PixelArray2d_of_real output1(   input1.range().x_dim()*zoom,
//                                                input1.range().y_dim()*zoom,
//                                                    0.0 );
//        PixelArray2d_of_real output2( output1 );
//
//        PixelArray2d_of_real output1_weights(output1);
//        PixelArray2d_of_real output2_weights(output2);
//
//        PixelArray2d_of_real output_double_drizzle( output1 );
//        PixelArray2d_of_real output_double_drizzle(output1);
//
//
//        drizzle::translate_and_drizzle_frame(
//                input1,  input1_weights,
//                output1, output1_weights,
//                input_translation1,
//                drizzle_scale,
//                pixfrac);
//
//        drizzle::translate_and_drizzle_frame(
//                input2,  input2_weights,
//                output2, output2_weights,
//                input_translation2,
//                drizzle_scale,
//                pixfrac);
//
//        drizzle::translate_and_drizzle_frame(
//                input1,  input1_weights,
//                output_double_drizzle, output_double_drizzle,
//                input_translation1,
//                drizzle_scale,
//                pixfrac);
//
//        drizzle::translate_and_drizzle_frame(
//                input2,  input2_weights,
//                output_double_drizzle, output_double_drizzle,
//                input_translation2,
//                drizzle_scale,
//                pixfrac);
//    }
}