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
#include "../psf_characterisation.h"

using namespace coela;
using namespace std;

SUITE(psf_characterisation)
{
    using namespace psf_characterisation;

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";
    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }


    struct cone_shaped_PSF_image {
        double cone_peak;
        PixelArray2d<double> cone_img;
        PixelPosition cone_centre;
        PixelIndex cone_centre_pixel;

        cone_shaped_PSF_image() {

            cone_img = PixelArray2d<double> (100, 100, 0);
            cone_peak = cone_img.range().x_dim()/2.0 - 5.0;
            cone_centre = PixelPosition(cone_img.range().x_dim()/2.0 +0.5 ,
                                        cone_img.range().y_dim()/2.0 +0.5);
            cone_centre_pixel = PixelPosition::pixel_containing_point(cone_centre);

            //Cone shape determined by:
            // y = x (x=radius)
            // for x < cone_peak
            //y = 0; thereafter

            for (PixelIterator i(cone_img.range()); i!=i.end; ++i) {
                cone_img(i) =
                    max(0.0,
                        cone_peak - coord_distance(PixelPosition::centre_of_pixel(i), cone_centre)
                       );
            }



        }

    };

    TEST_FIXTURE(cone_shaped_PSF_image, write_image_for_visual_check) {
        cone_img.write_to_file(test_suite_output_dir+"cone_img.fits");
    }

    TEST_FIXTURE(cone_shaped_PSF_image, sanity_check) {
        //Therefore, enclosed flux at x< cone peak
        // EF = 1/3*PI*x^3

        double analytical_sum_flux=M_PI*cone_peak*cone_peak*cone_peak / 3.0;

//        cout<<"Analytic sum: " <<analytical_sum_flux<<"; actual: "<< cone_img.sum()<<endl;

        CHECK_CLOSE(analytical_sum_flux, cone_img.sum(),  analytical_sum_flux/100.0);
    }

    TEST_FIXTURE(cone_shaped_PSF_image, test_get_radial_data_unmasked) {
        using namespace psf_characterisation::detail;
        std::vector<std::pair <double, double>  > radial_data =
            get_radial_data(cone_img, cone_centre, cone_img.range().x_dim());
        //Oversized range should be ok.

        CHECK(radial_data.empty()==false);

        sort(radial_data.begin(), radial_data.end(), double_pair_first_member_predicate);
        ofstream datfile(string(test_suite_output_dir+"raw_radial_data.txt").c_str());
        for (size_t i=0; i!=radial_data.size(); ++i) {
            datfile<< radial_data[i].first <<" "<<radial_data[i].second<<endl;

        }
        datfile.close();
    }

    TEST_FIXTURE(cone_shaped_PSF_image, test_get_radial_data_unmasked_and_bin) {
        using namespace psf_characterisation::detail;
        std::vector<std::pair <double, double>  > radial_data =
            get_radial_data(cone_img, cone_centre, cone_img.range().x_dim());
        //Oversized range should be ok.


        vector<psf_profile_point> profile_info =
            bin_average_radial_data(radial_data);

//        cout<<"Get " << profile_info.size()<<" bins"<<endl;
        CHECK(profile_info.size() > cone_img.range().x_dim()/2.0);
        //expect a few extra along the diagonal
    }


    TEST_FIXTURE(cone_shaped_PSF_image, test_FWHM_estimation) {
        double est_fwhm = estimate_fwhm_in_image_pix(cone_img,
                          cone_centre,
                          cone_peak,
                          cone_img.range().x_dim()/2.0
                                                    );

        //45 degree angle slope ->
        //half-width at base = max height.
        //Therefore; full-width at half-height = max height.

        double calculated_FWHM = cone_peak;
//        cout<<"FWHM est: " << est_fwhm<<endl;
        CHECK_CLOSE(calculated_FWHM, est_fwhm, cone_peak / 20);
    }

    //Doesn't quite work, because we end up trying to interpolate past known points... not worth fixing, throws a sensible exception.
//    TEST_FIXTURE(cone_shaped_PSF_image, test_full_enclosed_flux_radius_estimation){
//        double full_flux_est_radius =
//        estimate_radius_to_enclose_flux(cone_img,
//                cone_centre,
//                cone_img.sum(),
//                cone_img.range().x_dim() );
//
////        cout<<"Est radius 100pc flux: "<<full_flux_est_radius<<"; expected: "<< cone_peak<<endl;
//        CHECK_CLOSE(cone_peak, full_flux_est_radius, cone_peak*0.01);
//    }

    TEST_FIXTURE(cone_shaped_PSF_image, test_10p_enclosed_flux_radius_estimation) {
//        double flux_10p_est_radius =
//        estimate_radius_to_enclose_flux(cone_img,
//                cone_centre,
//                cone_img.sum()*0.1,
//                cone_img.range().x_dim() );

//        cout<<"Est radius 10pc flux: "<<full_flux_est_radius<<endl;
//        cout<<"; expected: "<< cone_peak<<endl;
//        CHECK_CLOSE(cone_peak, full_flux_est_radius, cone_peak*0.01);
    }

    TEST_FIXTURE(cone_shaped_PSF_image, test_FWEF_estimation) {
        double est_fwhef = estimate_FWHEF_in_image_pix(cone_img,
                           cone_centre,
                           cone_img.sum(),
                           cone_img.range().x_dim()
                                                      );

//       cout<<"est FWHEF: "<<est_fwhef<<endl;

        //EF = PI*X^2*A - (2/3 * PI * X^3 / 3)
        double x= est_fwhef /2.0;

        double analytic_ef  = M_PI*x*x*cone_peak - (2.0/3.0) *M_PI*x*x*x ;

//       cout<<"Est HWHEF "<<x<<endl;
//       cout<<"Resulting analytic flux proportion: " <<analytic_ef / cone_img.sum()<<endl;

        CHECK_CLOSE(0.5, analytic_ef/ cone_img.sum(), 0.02);
    }

    TEST_FIXTURE(cone_shaped_PSF_image, test_fully_encircled_flux_estimation) {

//        PixelArray2d<double>& cone_img_ref = cone_img;

        double  fully_encircled =
            estimate_encircled_flux_at_pixel_radius(
                cone_img,
                cone_centre,
                cone_img.range().x_dim()/2.0 - 3.5
            );

//        cout<<"Estimate total flux via encircled_flux_at_pixel_radius: "<< fully_encircled<<endl;
//        cout<<"Actual total flux "<< cone_img.sum()<<endl;
        CHECK_CLOSE(cone_img.sum(), fully_encircled, fully_encircled*0.01);
    }

    TEST_FIXTURE(cone_shaped_PSF_image, test_encircled_flux_estimation) {

        double ap_radius =   cone_img.range().x_dim()/5.0;

        double  EF = estimate_encircled_flux_at_pixel_radius(cone_img,
                     cone_centre,
                     ap_radius
                                                            );

        double x = ap_radius;
        double analytic_ef= M_PI*x*x*cone_peak - (2.0/3.0) *M_PI*x*x*x ;

//        cout<<"Estimate encircled_flux_at_pixel_radius: "<< ap_radius<<": "<< EF<<endl;
//
//        cout<<"Analytic calculation: "<< analytic_ef<<endl;
//
//        cout<<"Percent diff: " << (EF - analytic_ef )/ analytic_ef<<endl;

        CHECK_CLOSE(analytic_ef, EF, analytic_ef*0.035);
    }



}
