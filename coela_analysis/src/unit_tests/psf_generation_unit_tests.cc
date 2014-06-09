#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>
#include <boost/filesystem.hpp>
#include <fstream>

#include "coela_core/src/pixel_array_routines.h"
#include "coela_analysis/src/psf_generation.h"
//#include "../coela_analysis/src/psf_characterisation.h"

using namespace coela;
using namespace std;

SUITE(psf_generation)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(gaussian_generation) {
        PixelRange img_size(1,1,50,50);
        double peak=5.0, sigma=1.2;
        psf_models::gaussian_psf_model g_model(peak, sigma);

        CCD_Position true_centre(img_size.x_dim()/2.0 + 0.5, img_size.y_dim()/2.0 +0.5);

        psf_models::reference_psf ref_psf =
            psf_models::generate_psf(g_model, img_size,
                                     CCD_Position(0,0),
                                     true_centre,
                                     1.0,
                                     2);

        string filename = test_suite_output_dir + "example_gauss_gen.fits";
        ref_psf.psf_image.write_to_file(filename);

        PixelPosition centroid = pixel_array_routines::centroid(ref_psf.psf_image.pix,
                                 ref_psf.psf_image.pix.range());

        CCD_Position CCD_centroid = ref_psf.psf_image.CCD_grid.
                                    corresponding_grid_Position(centroid);

        CHECK_CLOSE(CCD_centroid.x, true_centre.x, 0.01);
        CHECK_CLOSE(CCD_centroid.y, true_centre.y, 0.01);

        double max_pixel = ref_psf.psf_image.pix(ref_psf.psf_image.pix.max_PixelIndex());

        CHECK_CLOSE(max_pixel, g_model.peak_val, 0.25);

        double est_flux = 2*M_PI*g_model.sigma_in_CCD_pix*g_model.sigma_in_CCD_pix *
                          g_model.peak_val;
        double actual_flux = ref_psf.psf_image.pix.sum();
//        cerr<<"Est flux:"<<est_flux<<"; actual: "<<actual_flux<<endl;
        CHECK_CLOSE(est_flux, actual_flux, est_flux*0.01);
    }

    TEST(airy_generation) {
        double lambda = 7.7e-7;
        double d=2.5;
        double nyquist_pixel_width_rads = lambda / d / 2.0;

        double central_obscuration=0.0;
        psf_models::airy_psf_model airy_gen(
            1,
            nyquist_pixel_width_rads,
            lambda,
            d,
            central_obscuration
        );

        double first_minima_in_CCD_pix =
            psf_models::find_first_minima_in_CCD_pix(airy_gen);
//        cout<<"minima " <<first_minima<<endl;
        double unit_peak_flux = psf_models::calculate_total_flux(airy_gen);

        double first_minima_flux =
            psf_models::calculate_flux_enclosed_at_radius(airy_gen,
                    first_minima_in_CCD_pix);

        double ratio = first_minima_flux  / unit_peak_flux;
//        cout<<"Flux ratio" << ratio <<endl;

        if (central_obscuration==0.0) {
//            cout<<"Circular aperture..."<<endl;
            CHECK_CLOSE(1.22*2.0, first_minima_in_CCD_pix, 0.01);
            CHECK_CLOSE(0.838, ratio, 0.01); //83.8% flux inside first minima (Wikipedia)
        } else if (central_obscuration==0.5) {
            CHECK_CLOSE(1.16*2.0, first_minima_in_CCD_pix, 0.01);
        }


    }
}
