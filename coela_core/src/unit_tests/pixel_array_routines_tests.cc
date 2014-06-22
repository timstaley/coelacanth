#include <UnitTest++/UnitTest++.h>
#include "../pixel_array_routines.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <sstream>
#include <cmath>
using namespace coela;
using namespace std;

SUITE(pixel_array_routines)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";
    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(RegionCentroid) {
        PixelArray2d<double> img(20,20, 0.0);

        PixelIndex first_pixel(5,10);
        img(first_pixel) = 1.0;

        CHECK(PixelPosition::centre_of_pixel(first_pixel) ==
              pixel_array_routines::centroid(img, img.range())) ;

        PixelIndex second_pixel(15,10);
        img(second_pixel) = 1.0;

        PixelPosition avg_posn = PixelPosition::mid_point(
                                     PixelPosition::centre_of_pixel(first_pixel) ,
                                     PixelPosition::centre_of_pixel(second_pixel));

        CHECK(avg_posn == pixel_array_routines::centroid(img, img.range()));

        PixelRange cutoff_box = img.range();
        cutoff_box.high.x = second_pixel.x -1;

        CHECK(PixelPosition::centre_of_pixel(first_pixel) ==
              pixel_array_routines::centroid(img,cutoff_box));

    }

    TEST(bin_image) {
        //Tested complete with co-ord system under image_utils tests.
    }




}
