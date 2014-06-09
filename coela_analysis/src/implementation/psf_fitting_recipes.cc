#include "../psf_fitting_recipes.h"
#include "../psf_characterisation.h"
#include "../regions.h"
#include <Minuit2/MnPrint.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMinimize.h>
#include <iostream>


namespace coela {
namespace psf_fitting {

Mn2_models::PSF_Fit<psf_models::gaussian_psf_model> fit_gaussian_model(
    const CCDImage<double>& background_subtracted_image,
    const CCD_Position& estimated_centre,
    const variance_estimation_base& noise_model,
    const double fitting_radius_in_CCD_pix,
    const int oversampling,
    bool verbose_display)
{

    using std::cout;
    using std::cerr;
    using std::endl;
    const CCDImage<double> & img=background_subtracted_image;

    Mn2_models::PSF_Fit<psf_models::gaussian_psf_model>
    initial_fit;
    initial_fit.Position.centre = estimated_centre;
    initial_fit.Position.error_estimate = CCD_PixelShift(0.75,0.75);

    PixelPosition init_centre_posn =
        img.CCD_grid.corresponding_pixel_Position(estimated_centre);

    PixelIndex init_centre_pix =
        PixelPosition::pixel_containing_point(init_centre_posn);

    initial_fit.model.peak_val =img.pix(init_centre_pix);

    double max_radius_for_fwhm_estimation = 20.0;

    double init_fwhm_img_pix =
        psf_characterisation::estimate_fwhm_in_image_pix(img.pix,
                init_centre_posn, initial_fit.model.peak_val,
                max_radius_for_fwhm_estimation);

    initial_fit.model.sigma_in_CCD_pix =
        init_fwhm_img_pix * img.CCD_grid.pixel_width_;

    ROOT::Minuit2::MnUserParameters init_pars;
    Mn2_models::gaussian_psf g_model;
    string prefix = "Test_Model_";

    g_model.set_MnPars_from_PSF_Fit(initial_fit, &init_pars, prefix);

    if (verbose_display) { cerr<<"Init Gauss pars:\n"<<init_pars<<endl; }

    regions::circular_aperture<coordinate_types::CCD> fitting_aperture(
        estimated_centre, fitting_radius_in_CCD_pix);

    CCDImage<double> fitting_mask =
        fitting_aperture.generate_mask_for_covered_portion_of_image(img);

//    fitting_mask.write_to_file("fitting_mask.fits");

    psf_fitting::analytic_psf_fit_cost_fcn psf_metric(
        img,
        fitting_mask,
        g_model,
        noise_model,
        oversampling);

    ROOT::Minuit2::MnMinimize minimizer(psf_metric, init_pars);
    if (verbose_display) { cerr<<"Running minimisation..."<<endl; }
    ROOT::Minuit2::FunctionMinimum results = minimizer();

    if (verbose_display) { cerr<<results<<endl; }

    return g_model.pull_PSF_Fit_from_MnPars(results.UserParameters(), prefix);
}

} //end namespace coela::psf_fitting
}//end namespace coela
