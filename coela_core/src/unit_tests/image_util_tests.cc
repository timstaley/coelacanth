#include <UnitTest++/UnitTest++.h>
#include "../pixel_array_routines.h"
#include "../image_utils.h"

#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <cmath>
using namespace coela;
using namespace std;

SUITE(image_utils)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";
    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }



    TEST(bin_image) {
        CCDImage<double> input_image;
        input_image.pix = PixelArray2d<double>(42,42, 0.0);

        input_image.pix(1,1)=5;
        input_image.pix(1,2)=5;
        input_image.pix(2,1)=5;
        input_image.pix(2,2)=5;

        input_image.pix(10,10)=20;
        input_image.pix(20,20)=55;

        input_image.initialize_CCD_grid_for_raw_data();

        CCDImage<double> bin_1x = image_utils::bin_image(input_image, 1);
        CCDImage<double> bin_2x = image_utils::bin_image(input_image, 2);
        CCDImage<double> bin_3x = image_utils::bin_image(input_image, 3);
        CCDImage<double> bin_42x = image_utils::bin_image(input_image, 42);

        CHECK(bin_1x.pix==input_image.pix);
//        CHECK(bin_1x==input_image);

        CHECK_EQUAL(input_image.CCD_grid.image_outline_, bin_1x.CCD_grid.image_outline_);
        CHECK_EQUAL(input_image.CCD_grid.image_outline_, bin_2x.CCD_grid.image_outline_);
        CHECK_EQUAL(input_image.CCD_grid.image_outline_, bin_3x.CCD_grid.image_outline_);

        CHECK_EQUAL(input_image.pix.sum(), bin_2x.pix.sum());
        CHECK_EQUAL(input_image.pix.sum(), bin_3x.pix.sum());

        CHECK_EQUAL(20.0, bin_2x.pix(1,1));

//        input_image.write_to_file("rebin_input.fits");
//        bin_2x.write_to_file("rebin_2x.fits");
//        bin_3x.write_to_file("rebin_3x.fits");
//        bin_42x.write_to_file("rebin_42x.fits");
    }

    TEST(resample_image) {
        CCDImage<double> input_image;
        input_image.pix = PixelArray2d<double>(40,40, 0.0);
        input_image.initialize_CCD_grid_for_raw_data();

        input_image.pix(20,20)=5;
        input_image.pix(20,21)=5;
        input_image.pix(21,20)=5;
        input_image.pix(21,21)=5;

        CCDImage<double> downsample4x = image_utils::bin_image(input_image, 4);
        CCDImage<double> downsample2x = image_utils::bin_image(input_image, 2);


        CCDImage<double> interp4x = image_utils::bicubic_resample_4x(downsample4x);
        CCDImage<double> interp2x = image_utils::bicubic_resample_2x(downsample2x);

        CHECK_EQUAL(input_image.pix.sum(), interp2x.pix.sum());
        CHECK_EQUAL(input_image.pix.sum(), interp4x.pix.sum());

        CCD_Position orig_centroid = input_image.CCD_grid.corresponding_grid_Position(
                                         pixel_array_routines::centroid(
                                             input_image.pix, input_image.pix.range())
                                     );

        CCD_Position centroid_4x = interp4x.CCD_grid.corresponding_grid_Position(
                                       pixel_array_routines::centroid(interp4x.pix, interp4x.pix.range())
                                   );

        CHECK_CLOSE(orig_centroid.x, centroid_4x.x, 1e-6);
        CHECK_CLOSE(orig_centroid.y, centroid_4x.y, 1e-6);

//        input_image.write_to_file("resample_input.fits");
//        interp4x.write_to_file("resample_4x.fits");
//        interp2x.write_to_file("resample_2x.fits");
    }


}
