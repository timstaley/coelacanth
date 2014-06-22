
#include <UnitTest++/UnitTest++.h>
#include "../ccd_dataset_info.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <exception>
#include <boost/filesystem.hpp>

using namespace coela;
using namespace std;


extern string lucky_lib_test_resources_dir;

SUITE(CCD_dataset_info)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }


    TEST(serialization_to_stream) {

        CCD_DatasetInfo cdi1;
        cdi1.ccd_id=42;
        cdi1.CCD_inputdir = "/home/test/ccd42";
        cdi1.CCD_outputdir = "/home/test_out/ccd42";
        cdi1.in_filestem = "run123";
        cdi1.out_filestem = "run123";
        cdi1.extension = ".fits";
        cdi1.dataset_output_base_dir = "/home/test_out";

        cdi1.default_camera_config_file= "/home/settings/camconf";
        cdi1.default_guide_star_region_file= "/home/test_out/gs.txt";
        cdi1.default_faint_histogram_region_file= "/home/test_out/gain.txt";

        cdi1.most_recently_output_frame_list= "/home/test_out/frame_list.txt";

        std::stringstream ss1, ss2;
        string str1, str2;

        ss1<<cdi1;
        str1=ss1.str();

        CCD_DatasetInfo cdi2;
        ss1>>cdi2;

        ss2 << cdi2;
        str2=ss2.str();

        CHECK_EQUAL(str1,str2);


    }


}