#include <UnitTest++/UnitTest++.h>
#include "../ds9_region.h"

#include <iostream>
#include <sstream>

#include <cmath>
using namespace coela;
using namespace std;
SUITE(Ds9Region)
{

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"Ds9Region\" unit tests running ***" <<endl;
    }

    TEST(check_linkage_ok) {
        Position<coordinate_types::pixels> x(5,5);
        //Try some basic constructors to make sure all the template syntax is correct
        ds9::Ds9Region x_pt = ds9::Ds9Region::point(x);
    }

    TEST(conversion) {
        Position<coordinate_types::pixels> x1(5,5);
        ds9::Ds9Region ds9_x1 = ds9::Ds9Region::point(x1);
        Position<coordinate_types::pixels> x2(ds9_x1);
        CHECK(x1==x2);
    }

    TEST(serialization_coord_pixels) {
        Position<coordinate_types::pixels> x1(5,5);
        ds9::Ds9Region ds9_x1 = ds9::Ds9Region::point(x1);

        stringstream ss;
        ss<<ds9_x1;

//        cout<<"point string:\n"<<ss.str()<<endl;


        //NB the coord-type does not get stored in the serialization,
        //so we will just init to the right value:
        ds9::Ds9Region ds9_x2(ds9_x1);


        ss>>ds9_x2;

        Position<coordinate_types::pixels> x2(ds9_x2);

        CHECK(x1==x2);
    }



    TEST(file_save_load_pixel_location) {
        Position<coordinate_types::pixels> x1(5,5);
        Position<coordinate_types::pixels> x2(10,5);

        ds9::Ds9Region ds9_x1 = ds9::Ds9Region::point(x1);
        ds9::Ds9Region ds9_x2 = ds9::Ds9Region::point(x2);

        vector<ds9::Ds9Region> v;
        v.push_back(ds9_x1);
        v.push_back(ds9_x2);

        string filename = "test_Ds9Region_pixel.reg";
        ds9::Ds9Region::save_regions_to_file(filename, v);

//        cout<<"Trying to load"<<endl;
        vector<ds9::Ds9Region> loaded = ds9::Ds9Region::load_regions_from_file(filename);

        CHECK_EQUAL(v.size(), loaded.size());
        CHECK(x1 == Position<coordinate_types::pixels>(v.front()));
        CHECK(x2 == Position<coordinate_types::pixels>(v.back()));
    }

    TEST(file_save_load_CCD_boxes) {
        CcdPosition x1(5,5);
        CcdPosition x2(10,5);
        CcdPosition x3(10,10);
        CcdBoxRegion box1(x1,x2);
        CcdBoxRegion box2(x2,x3);

        ds9::Ds9Region ds9_box1 = ds9::Ds9Region::box(box1);
        ds9::Ds9Region ds9_box2 = ds9::Ds9Region::box(box2);

        vector<ds9::Ds9Region> v;
        v.push_back(ds9_box1);
        v.push_back(ds9_box2);

        string filename = "test_Ds9Region_CCD_boxes.reg";
        ds9::Ds9Region::save_regions_to_file(filename, v);


        vector<ds9::Ds9Region> loaded = ds9::Ds9Region::load_regions_from_file(filename);


        CHECK_EQUAL(v.size(), loaded.size());
        CHECK(box1 == CcdBoxRegion(v.front()));
        CHECK(box2 == CcdBoxRegion(v.back()));
    }
}
