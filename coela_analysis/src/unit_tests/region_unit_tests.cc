#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>
#include <boost/filesystem.hpp>
#include <fstream>

#include "coela_utility/src/string_utils.h"
#include "coela_core/src/ccd_image.h"
#include "../regions.h"

using namespace coela;
using namespace std;

SUITE(regions)
{
    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";
    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(basic_functionality) {
        PixelPosition centre(5.25,5.25);
        double radius=2.5;
        regions::circular_aperture<coordinate_types::pixels> ap1(centre, radius);

        CCDImage<double> eight_pix;
        eight_pix.pix = PixelArray2d<double>(8,8,0.0);
        eight_pix.initialize_CCD_grid_for_raw_data();
        CCDImage<double> eight_pix_mask =
            ap1.generate_mask_for_covered_portion_of_image(eight_pix);
        string filename = test_suite_output_dir + "eight_pix_img.fits";
        eight_pix.write_to_file(filename);
        filename = test_suite_output_dir + "eight_pix_mask.fits";
        eight_pix_mask.write_to_file(filename);

        CCDImage<double> five_pix;
        five_pix.pix = PixelArray2d<double>(5,5,0.0);
        five_pix.initialize_CCD_grid_for_raw_data();
        CCDImage<double> five_pix_mask = ap1.generate_mask_for_covered_portion_of_image(five_pix);
        filename = test_suite_output_dir + "five_pix_img.fits";
        five_pix.write_to_file(filename);
        filename = test_suite_output_dir + "five_pix_overlapping_mask.fits";
        five_pix_mask.write_to_file(filename);
    }

}



