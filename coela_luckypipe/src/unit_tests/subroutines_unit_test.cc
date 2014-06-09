#include <UnitTest++/UnitTest++.h>

#include "../clean_and_register_subroutines.h"
#include "coela_core/src/convolution.h"
#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>


using namespace coela;
using namespace std;

SUITE(Clean_And_Register_Subroutines)
{
    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    struct camconf_fixture {
        camconf_fixture() {
            camconf.lens_inf.telescope_inner_diameter=0.5;
            camconf.lens_inf.telescope_outer_diameter=2.5;
            camconf.lens_inf.nominal_pixel_scale_in_mas=32.5;
            OpticalFilterInfo f;
            f.ccd_id=0;
            f.central_wavelength_in_metres=7.7e-7;
            f.name="i";
            camconf.filters.push_back(f);
        }
        CameraConfigInfo camconf;

    };



    TEST_FIXTURE(camconf_fixture, Airy_template) {
        double resample_factor=4;


        psf_models::reference_psf kernel =
            clean_and_register_subroutines::generate_airy_core_template(
                camconf,
                0,
                1.0/resample_factor
            );

        kernel.psf_image.write_to_file(test_suite_output_dir+"airy_core.fits");

        double fwhm_est = psf_characterisation::estimate_fwhm_in_image_pix(kernel.psf_image.pix,
                          kernel.exact_centre,
                          kernel.psf_image.pix.max_val(),
                          kernel.psf_image.pix.range().x_dim()/2.0);

        cout<<"Airy fwhm in pix:"<<fwhm_est / resample_factor << endl;



        CHECK_CLOSE(2.0, fwhm_est / resample_factor,
                    0.05); //FWHM approx 1.02 lambda / D, (aperture inner diameter=0)

        PixelArray2d<double> airy_autoconvolution =
            spatial_convolution::convolve_with_kernel(
                kernel.psf_image.pix,
                kernel.psf_image.pix,
                kernel.central_pixel);

        airy_autoconvolution.write_to_file(test_suite_output_dir+"airy_autoconv.fits");
        double autocon_fwhm_est = psf_characterisation::estimate_fwhm_in_image_pix(
                                      airy_autoconvolution,
                                      kernel.exact_centre,
                                      airy_autoconvolution.max_val(),
                                      airy_autoconvolution.range().x_dim()/2.0);
        cout<<"Autoconv fwhm: "<< autocon_fwhm_est / resample_factor<<endl;

    }

    TEST_FIXTURE(camconf_fixture, custom_template) {
        double resample_factor=4;
        psf_models::reference_psf kernel =
            clean_and_register_subroutines::generate_normalised_core_template(
                camconf,
                0,
                1.0/resample_factor
            );

        kernel.psf_image.write_to_file(
            test_suite_output_dir+"normalised_core.fits");


    }
}
