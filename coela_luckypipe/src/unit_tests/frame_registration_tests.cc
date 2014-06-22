#include <UnitTest++/UnitTest++.h>
#include "../frame_registration.h"
#include "coela_core/src/pixel_array_routines.h"
#include "coela_core/src/image_utils.h"
#include <boost/filesystem.hpp>

#include <iostream>
#include <sstream>


using namespace coela;
using namespace std;
SUITE(frame_registration)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(find_best_psf_match) {

        CcdImage<double> construction_image;
        construction_image.pix = PixelArray2d<double>(200,200,0.0);

        construction_image.initialize_CCD_grid_to_specific_offset_and_scale(
            CcdPosition(0,0), 0.25);

        psf_models::ReferencePsf ref_psf4x;

        bool output_to_screen=false;

        ref_psf4x.psf_image.pix =PixelArray2d<double>(40, 40, 0);
        ref_psf4x.exact_centre = PixelPosition(20.61,20.61);
        ref_psf4x.central_pixel = PixelPosition::pixel_containing_point(ref_psf4x.exact_centre);
        ref_psf4x.psf_image.initialize_CCD_grid_to_specific_offset_and_scale(CcdPosition(0,0),
                0.25);
        for (PixelIterator i(ref_psf4x.psf_image.pix.range()); i!=i.end; ++i) {
            ref_psf4x.psf_image.pix(i) =
                max(0.0,
                    10 - coord_distance(PixelPosition::centre_of_pixel(i), ref_psf4x.exact_centre)
                   );
        }

        PixelIndex output_offset(83,71);
        PixelRange output_range(ref_psf4x.psf_image.pix.range().low + output_offset,
                                ref_psf4x.psf_image.pix.range().high + output_offset);
        for (PixelIterator in(ref_psf4x.psf_image.pix.range()), out(output_range);
                in!=in.end && out!=out.end;
                ++in, ++out) {
            construction_image.pix(out)+=ref_psf4x.psf_image.pix(in);
        }

        construction_image.write_to_file(test_suite_output_dir+"construction_img.fits");
        ref_psf4x.psf_image.write_to_file(test_suite_output_dir+"ref_psf4x.fits");

        psf_models::ReferencePsf ref_psf1x = ref_psf4x;
        ref_psf1x.exact_centre.x/=4.0;  ref_psf1x.exact_centre.y/=4.0;
        ref_psf1x.central_pixel = PixelPosition::pixel_containing_point(ref_psf1x.exact_centre);
        ref_psf1x.psf_image = image_utils::bin_image(ref_psf4x.psf_image, 4);

//
        PixelPosition expected_img_centre = ((ref_psf4x.exact_centre - PixelPosition::origin)
                                             + PixelShift(output_offset.x, output_offset.y)) + PixelPosition::origin;

        CHECK_CLOSE(pixel_array_routines::centroid(construction_image.pix,
                    construction_image.pix.range()).x ,
                    expected_img_centre.x,
                    0.01);

        CcdPosition known_image_centre =
            construction_image.CCD_grid.corresponding_grid_Position(expected_img_centre);


        CcdImage<double> test_image = image_utils::bin_image(construction_image, 4);
        test_image.write_to_file(test_suite_output_dir+"binned_img.fits");

        CHECK_CLOSE(1.0, test_image.CCD_grid.pixel_width_, 1e-6);

        CcdBoxRegion search_region(5,5,45,45);

        frame_registration::GuideStarLock<CcdPosition> lock1x =
            frame_registration::find_best_psf_match(test_image,
                    search_region,
                    ref_psf1x,
                    1.0,
                    0.25
                                                   );

        CHECK_CLOSE(known_image_centre.x, lock1x.Position.x, 0.5);
        CHECK_CLOSE(known_image_centre.y, lock1x.Position.y, 0.5);

        if (output_to_screen) {
            cerr<<"1x error: " << coord_distance(known_image_centre  ,lock1x.Position)<<endl;
            cerr<<"1x error xy: " << known_image_centre  - lock1x.Position<<endl;
        }
        frame_registration::GuideStarLock<CcdPosition> lock1x_parabolic =
            frame_registration::find_best_psf_match(test_image,
                    search_region,
                    ref_psf1x,
                    1.0,
                    0.25,
                    true
                                                   );

        CHECK_CLOSE(known_image_centre.x, lock1x_parabolic.Position.x, 0.5);
        CHECK_CLOSE(known_image_centre.y, lock1x_parabolic.Position.y, 0.5);
        if (output_to_screen) {
            cerr<<"1x parabolic error: " << coord_distance(known_image_centre  ,
                    lock1x_parabolic.Position)<<endl;
            cerr<<"1x parabolic error xy: " << known_image_centre  - lock1x_parabolic.Position<<endl;
        }

        frame_registration::GuideStarLock<CcdPosition> lock4x =
            frame_registration::find_best_psf_match(test_image,
                    search_region,
                    ref_psf4x,
                    4,
                    0.25
                                                   );

        CHECK_CLOSE(known_image_centre.x, lock4x.Position.x, 0.13);
        CHECK_CLOSE(known_image_centre.y, lock4x.Position.y, 0.13);
        if (output_to_screen) {
            cerr<<"4x error: " << coord_distance(known_image_centre ,lock4x.Position)<<endl;
            cerr<<"4x error xy: " << known_image_centre - lock4x.Position<<endl;
        }
        frame_registration::GuideStarLock<CcdPosition> lock4x_parabolic =
            frame_registration::find_best_psf_match(test_image,
                    search_region,
                    ref_psf4x,
                    4,
                    0.25,
                    true
                                                   );
        if (output_to_screen) {
            cerr<<"4x parabolic error: " << coord_distance(known_image_centre ,
                    lock4x_parabolic.Position)<<endl;
            cerr<<"4x paraboloic error xy: " << known_image_centre - lock4x_parabolic.Position<<endl;
        }
//        CHECK_CLOSE(known_image_centre.x, lock4x_parabolic.Position.x, 0.05);
//        CHECK_CLOSE(known_image_centre.y, lock4x_parabolic.Position.y, 0.05);


    }


}