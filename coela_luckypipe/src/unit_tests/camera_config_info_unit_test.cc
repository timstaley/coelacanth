#include <UnitTest++/UnitTest++.h>

#include "../camera_config_info.h"

#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>


using namespace coela;
using namespace std;
SUITE(CameraConfigInfo)
{

    string test_suite_output_dir = "camera_config_unit_tests/";

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"camera_config_info\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }


    TEST(file_save_and_load) {
        CameraConfigInfo cci;

        CcdCalibrationInfo ccd_props;
        ccd_props.ccd_id=0;
        ccd_props.precal_row_bias_frame_available=false;
        ccd_props.precal_row_bias_frame_path="none";
        ccd_props.cropped_PixelRange = PixelRange(10,10,500,500);
        ccd_props.crop_region = CcdBoxRegion(9.0,9.0,500,500);
        ccd_props.default_temporal_debiasing_histogram_region =
            CcdBoxRegion(250,250,500,500);
        ccd_props.dark_current_frames_available=false;
        ccd_props.normalised_DC_frame_path="none";
        ccd_props.thresholded_DC_frame_path="none";

        cci.CCD_vec.push_back(ccd_props);


        cci.lens_inf.nominal_pixel_scale_in_mas=55.3;
        cci.lens_inf.telescope_inner_diameter=0.5;
        cci.lens_inf.telescope_outer_diameter=2.5;
        cci.lens_inf.ccd_mosaic_offsets.push_back(
            pair<int, MosaicPixelShift>(int(0), MosaicPixelShift(50,50)));

        OpticalFilterInfo filter;
        filter.ccd_id=0;
        filter.central_wavelength_in_metres=770e-9;
        filter.name="I band";
        cci.filters.push_back(filter);

        cci.simulated_data=false;

        string filename = test_suite_output_dir+"camconf.txt";
        cci.write_to_files(filename);

        //Does it load?
        CameraConfigInfo cci2(filename);
        cci2.write_to_files(test_suite_output_dir+"camconf2.txt");

        //TO DO - check that loaded version actually matches saved version


    }

}

