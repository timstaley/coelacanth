#include <UnitTest++/UnitTest++.h>

#include "../drizzle_settings.h"

#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>


using namespace coela;
using namespace std;
SUITE(DrizzleSettings)
{

    string test_suite_output_dir = "DrizzleSettings_unit_tests/";

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"DrizzleSettings\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }



    TEST(file_save_and_load) {
        DrizzleSettings ds1;

        ds1.output_base_folder="output/";

        ds1.CCD_ids_present.push_back(0);
        ds1.CCD_output_sub_folders.push_back("CCD_0");

        ds1.CCD_ids_present.push_back(2);
        ds1.CCD_output_sub_folders.push_back("CCD_2");

        ds1.normalisation_on=true;
        ds1.thresholding_on=false;
        ds1.post_drizzle_col_debias_on=true;
        ds1.create_drizzle_mosaic=true;

        ds1.drizzle_scale_factor=1.0;
        ds1.drizzle_pixel_fraction=0.6;

        ds1.output_percentiles.push_back(0.025);
        ds1.output_percentiles.push_back(0.05);
        ds1.output_percentiles.push_back(0.1);
        ds1.output_percentiles.push_back(0.25);
        ds1.output_percentiles.push_back(0.5);
        ds1.output_percentiles.push_back(0.75);
        ds1.output_percentiles.push_back(1.0);

        ds1.write_to_file(test_suite_output_dir+"ds1test.txt");
        DrizzleSettings ds2(test_suite_output_dir+"ds1test.txt");
        ds2.write_to_file(test_suite_output_dir+"ds2test.txt");

    }
}
