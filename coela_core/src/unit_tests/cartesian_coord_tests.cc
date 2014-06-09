#include <UnitTest++/UnitTest++.h>
#include "../cartesian_coords.h"

#include <iostream>
#include <sstream>

#include <cmath>
using namespace coela;
using namespace std;
SUITE(position_class)
{

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"position\" unit tests running ***" <<endl;
    }

//    TEST(not_a_number){
//        cout << "This test will fail when optimizations are on - this is OK" <<endl;
//        double x = NAN;
//        CHECK(  std::isnan(x));
//
//        double y = 5.0;
//        CHECK(  !std::isnan(y));
//    }

    TEST(constructor) {
        Position<coordinate_types::pixels> null_pixel_Position;

        Position<coordinate_types::pixels> valid_pixel_Position(0,0);

        double x = NAN;
        //Gets disabled when math optimizations are on
        bool isnan_currently_enabled = std::isnan(x);

        if (isnan_currently_enabled) {
            CHECK(null_pixel_Position.is_valid()==false);
        }

        CHECK(valid_pixel_Position.is_valid()==true);
    }

    TEST(pixels_only_function) { //More to check the syntax is all ok, than anything else.
        PixelIndex ind(5,5);

        PixelPosition pixcentre = PixelPosition::centre_of_pixel(ind);
        PixelPosition pixcentre_manual(ind.x-0.5, ind.y-.5);
        CHECK(pixcentre==pixcentre_manual);

    }

    TEST(coord_distance_squared) {
        Position<coordinate_types::pixels> origin(0,0);
        Position<coordinate_types::pixels> five(5,0);

        CHECK_EQUAL(25, coord_distance_squared(origin,five));
    }
}

SUITE(RectangularRegion)
{

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"RectangularRegion\" unit tests running ***" <<endl;
    }


    TEST(constructor) {
        RectangularRegion<coordinate_types::pixels> null_box;

        RectangularRegion<coordinate_types::pixels> valid_box(0,0, 10,10);

        double x = NAN;
        //Gets disabled when math optimizations are on
        bool isnan_currently_enabled = std::isnan(x);

        if (isnan_currently_enabled) {
            CHECK(null_box.is_valid()==false);
        }

        CHECK(valid_box.is_valid()==true);
    }

    TEST(basic_functionality) {
        RectangularRegion<coordinate_types::pixels> box(0,0 , 10,10);

        PixelPosition expected_centre(5,5);

        CHECK(expected_centre == box.centre());

    }

    TEST(Point_region_expansion) {
        PixelPosition centre(5.0,5.0);
        PixelBoxRegion box(centre,centre);

        double pad_width = 5.0;
        box = PixelBoxRegion::padded_copy(box, pad_width);
        CHECK(box.is_valid());
        CHECK_EQUAL(pad_width*2.0, box.x_dim());
        CHECK_EQUAL(pad_width*2.0, box.y_dim());
    }
}
