#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>
#include <boost/filesystem.hpp>
#include <fstream>

#include "coela_utility/src/string_utils.h"
#include "coela_utility/src/histogram_container.h"

#include "coela_core/src/ccd_image.h"
#include "coela_analysis/src/psf_models.h"
#include "coela_analysis/src/minuit_psf_models.h"

#include <Minuit2/MnPrint.h>

using namespace coela;
using namespace std;

SUITE(minuit_psf_models)
{


    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";
    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(Position_MnPar_interface) {
        psf_fitting::Mn2_models::PositionFit posn;
        posn.centre=CCD_Position(3.14,2.55);
        posn.error_estimate=CCD_PixelShift(0.04,0.03);

        ROOT::Minuit2::MnUserParameters pars;
        string prefix = "Testing_PositionFit_Interface_";
        psf_fitting::Mn2_models::PositionFit::set_MnPars(posn, &pars, prefix);

        psf_fitting::Mn2_models::PositionFit recovered_posn =
            psf_fitting::Mn2_models::PositionFit::pull_from_MnPars(pars, prefix);

        CHECK_EQUAL(posn.centre, recovered_posn.centre);
        CHECK_EQUAL(posn.error_estimate, recovered_posn.error_estimate);

    }


    TEST(gaussian_MnPar_interface) {
        CCD_Position centre(3.14,2.718);
        CCD_PixelShift error_est(0.4,0.5);

        psf_models::gaussian_psf_model g_model(5.5, 2.3);
        psf_fitting::Mn2_models::PSF_Fit<psf_models::gaussian_psf_model> g_fit;
        g_fit.model = g_model;
        g_fit.Position.centre=centre;
        g_fit.Position.error_estimate=error_est;

        ROOT::Minuit2::MnUserParameters pars;

        string prefix = "testing_mnpars_interface_";

        psf_fitting::Mn2_models::gaussian_psf::
        set_MnPars_from_PSF_Fit(g_fit, &pars, prefix);

//        cout<<pars<<endl;

        psf_fitting::Mn2_models::
        PSF_Fit<psf_models::gaussian_psf_model> pulled_fit =
            psf_fitting::Mn2_models::gaussian_psf::pull_PSF_Fit_from_MnPars(pars,  prefix);

        CHECK_EQUAL(g_fit.model.sigma_in_CCD_pix, pulled_fit.model.sigma_in_CCD_pix);
        //etc. (to do)

    }

}
