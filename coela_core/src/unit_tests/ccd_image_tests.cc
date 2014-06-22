
#include <UnitTest++/UnitTest++.h>
#include "../ccd_image.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <exception>

#include <boost/filesystem.hpp>

using namespace coela;
using namespace std;

extern string lucky_lib_test_resources_dir;

SUITE(CcdImage)
{
    string test_suite_output_dir = "image_class_tests/";

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"image\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
//        CHECK(false);
    }

    TEST(sub_image_ds9_coord_display) {
        CcdImage<double> img1;
        img1.pix = PixelArray2d<double> (100,100,55.5);
        img1.initialize_CCD_grid_for_raw_data();

        for (PixelIterator p(img1.pix.range()); p!=p.end; ++p) {
            img1.pix(p)=p.x *cos(p.y/10.0) ;
        }
        img1.write_to_file("subimage_full.fits");

        PixelRange box1(20,20,50,50);

        CcdImage<double> img2 = CcdImage<double>::sub_image(img1,box1);
        CHECK_EQUAL(box1.x_dim(), img2.pix.range().x_dim());
        CHECK_EQUAL(box1.y_dim(), img2.pix.range().y_dim());
//        cerr<<"image 2 outline: " <<img2.array.outline()<<endl;
        img2.write_to_file("subimage_sub.fits");
    }

    TEST(load_from_generic_file_type) {
        CcdImage<float> flt_img;
        flt_img.pix = PixelArray2d<float>(10,10,3.14);
        string filename = test_suite_output_dir+"float_img.fits";
        flt_img.write_to_file(filename);

        CcdImage<double> dbl_img =
            CcdImage<double>::load_from_unknown_filetype(filename);

        for (PixelIterator i(dbl_img.pix.range()); i!=i.end; ++i) {
            CHECK_CLOSE(flt_img.pix(i), dbl_img.pix(i), 1e-7);
        }

    }

    TEST(load_from_float_FITS_cube) {
        string fits_cube_filename = lucky_lib_test_resources_dir+"cube_float_file.fits";
        FitsHeader fht(fits_cube_filename);
        PixelCubeHeader cube_hdr(fht);

        CHECK_EQUAL(size_t(99), cube_hdr.z_dim());

        CcdImage<float> plane1_from_cube =
            CcdImage<float>::load_from_cube_FITS_file(fits_cube_filename,    1);

        CcdImage<float> plane45 =
            CcdImage<float>::load_from_cube_FITS_file(fits_cube_filename,   45);

        CcdImage<float> plane99_from_cube =
            CcdImage<float>::load_from_cube_FITS_file(fits_cube_filename,   99);

        plane1_from_cube.write_to_file(test_suite_output_dir+"cube_plane1.fits", fht);
        plane45.write_to_file(test_suite_output_dir+"cube_plane45.fits", fht);
        plane99_from_cube.write_to_file(test_suite_output_dir+"cube_plane99.fits", fht);


        CcdImage<float> plane1_via_ds9_extraction(lucky_lib_test_resources_dir
                +"cube_float_file_frame1.fits");

        for (PixelIterator i(plane1_from_cube.pix.range()); i!=i.end; ++i) {
            CHECK_CLOSE(plane1_via_ds9_extraction.pix(i), plane1_from_cube.pix(i), 1e-6);
        }

    }

    TEST(init_CCD_grid_for_raw_data) {
        CcdImage<float> test_image;
        test_image.pix = PixelArray2d<float>(5,5,0);
        test_image.initialize_CCD_grid_for_raw_data();
        CcdBoxRegion expected_rgn(0,0,5,5);
        CHECK_EQUAL(test_image.CCD_grid.image_outline_, expected_rgn);
    }

}

