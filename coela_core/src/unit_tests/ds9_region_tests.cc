#include <UnitTest++/UnitTest++.h>
#include "../DS9Region.h"

#include <iostream>
#include <sstream>

#include <cmath>
using namespace coela;
using namespace std;
SUITE(DS9Region)
{

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"DS9Region\" unit tests running ***" <<endl;
    }

    TEST(check_linkage_ok) {
        Position<coordinate_types::pixels> x(5,5);
        //Try some basic constructors to make sure all the template syntax is correct
        ds9::DS9Region x_pt = ds9::DS9Region::point(x);
    }

    TEST(conversion) {
        Position<coordinate_types::pixels> x1(5,5);
        ds9::DS9Region ds9_x1 = ds9::DS9Region::point(x1);
        Position<coordinate_types::pixels> x2(ds9_x1);
        CHECK(x1==x2);
    }

    TEST(serialization_coord_pixels) {
        Position<coordinate_types::pixels> x1(5,5);
        ds9::DS9Region ds9_x1 = ds9::DS9Region::point(x1);

        stringstream ss;
        ss<<ds9_x1;

//        cout<<"point string:\n"<<ss.str()<<endl;


        //NB the coord-type does not get stored in the serialization,
        //so we will just init to the right value:
        ds9::DS9Region ds9_x2(ds9_x1);


        ss>>ds9_x2;

        Position<coordinate_types::pixels> x2(ds9_x2);

        CHECK(x1==x2);
    }



    TEST(file_save_load_pixel_location) {
        Position<coordinate_types::pixels> x1(5,5);
        Position<coordinate_types::pixels> x2(10,5);

        ds9::DS9Region ds9_x1 = ds9::DS9Region::point(x1);
        ds9::DS9Region ds9_x2 = ds9::DS9Region::point(x2);

        vector<ds9::DS9Region> v;
        v.push_back(ds9_x1);
        v.push_back(ds9_x2);

        string filename = "test_DS9Region_pixel.reg";
        ds9::DS9Region::save_regions_to_file(filename, v);

//        cout<<"Trying to load"<<endl;
        vector<ds9::DS9Region> loaded = ds9::DS9Region::load_regions_from_file(filename);

        CHECK_EQUAL(v.size(), loaded.size());
        CHECK(x1 == Position<coordinate_types::pixels>(v.front()));
        CHECK(x2 == Position<coordinate_types::pixels>(v.back()));
    }

    TEST(file_save_load_CCD_boxes) {
        CCD_Position x1(5,5);
        CCD_Position x2(10,5);
        CCD_Position x3(10,10);
        CCD_BoxRegion box1(x1,x2);
        CCD_BoxRegion box2(x2,x3);

        ds9::DS9Region ds9_box1 = ds9::DS9Region::box(box1);
        ds9::DS9Region ds9_box2 = ds9::DS9Region::box(box2);

        vector<ds9::DS9Region> v;
        v.push_back(ds9_box1);
        v.push_back(ds9_box2);

        string filename = "test_DS9Region_CCD_boxes.reg";
        ds9::DS9Region::save_regions_to_file(filename, v);


        vector<ds9::DS9Region> loaded = ds9::DS9Region::load_regions_from_file(filename);


        CHECK_EQUAL(v.size(), loaded.size());
        CHECK(box1 == CCD_BoxRegion(v.front()));
        CHECK(box2 == CCD_BoxRegion(v.back()));
    }
}
