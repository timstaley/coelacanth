#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>
#include <boost/filesystem.hpp>
#include <fstream>

#include "coela_utility/src/string_utils.h"
#include "coela_core/src/ccd_image.h"
#include "coela_analysis/src/psf_generation.h"
#include "coela_analysis/src/psf_fitting_recipes.h"
#include "coela_random/src/random.h"

using namespace coela;
using namespace std;

SUITE(PsfFitting_recipes)
{


    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";
    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(Gaussian_fitting_recipe) {

        bool output_to_screen=false;
        PixelRange img_size(1,1,50,50);

        psf_models::GaussianPsfModel g_model(5.0, 1.2);

        vector<unsigned long> seed;
        seed.push_back(111);
        seed.push_back(122);
        seed.push_back(133);
        seed.push_back(144);
        seed.push_back(155);
        seed.push_back(166);

        unuran::StreamWrapper::
        set_unuran_package_seed(seed, false);
        unuran::StreamWrapper rns;
        unuran::UniformRandomVariate urv(10.,40., rns);

        CcdPosition g_model_true_centre(urv(), urv());

        cerr<<"CCD centre Position: " <<g_model_true_centre<<endl;

        psf_models::ReferencePsf test_psf_ref_img =
            psf_models::generate_psf(g_model, img_size,
                                     CcdPosition(0.0,0.0),
                                     g_model_true_centre,
                                     1.0,
                                     2,
                                     true, 10.0);

        string filename = test_suite_output_dir +"gauss_model.fits";
        test_psf_ref_img.psf_image.write_to_file(filename);

        PixelIndex peak_pixel = test_psf_ref_img.psf_image.pix.max_PixelIndex();
        CcdPosition initial_centre_estimate =
            test_psf_ref_img.psf_image.CCD_grid.corresponding_grid_Position(
                PixelPosition::centre_of_pixel(peak_pixel));


        psf_fitting::least_squares_weighting lsq_weight;

        psf_fitting::Mn2_models::PsfFit<psf_models::GaussianPsfModel>
        fit_result =
            psf_fitting::fit_gaussian_model(test_psf_ref_img.psf_image,
                                            initial_centre_estimate,
                                            lsq_weight,
                                            1.8,
                                            4,
                                            output_to_screen);

        CHECK_CLOSE(g_model.peak_val, fit_result.model.peak_val, g_model.peak_val * 0.02);
        CHECK_CLOSE(g_model.sigma_in_CCD_pix, fit_result.model.sigma_in_CCD_pix,
                    g_model.sigma_in_CCD_pix * 0.025);
        CHECK_CLOSE(g_model_true_centre.x, fit_result.Position.centre.x, 0.05);
        CHECK_CLOSE(g_model_true_centre.y, fit_result.Position.centre.y, 0.05);



    }

}
