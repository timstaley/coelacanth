#include <UnitTest++/UnitTest++.h>
#include <sstream>
#include "../image_grid.h"
#include <iostream>

using namespace coela;
using namespace std;

SUITE(ImageGrid)
{
    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"ImageGrid\" unit tests running ***" <<endl;
    }


    TEST(CCD_image_region_fht_read_write) {
        ImageGrid<coordinate_types::ccd> rgn1(CcdBoxRegion(0,0,10,10),5);

        FitsHeader fht;

        rgn1.write_to_fht(fht);
        ImageGrid<coordinate_types::ccd> rgn2;
        rgn2.load_from_fht(fht);

        CHECK(rgn1.image_outline_==rgn2.image_outline_);
        CHECK(rgn1.pixel_width_==rgn2.pixel_width_);
    }

    TEST(mosaic_image_region_fht_read_write) {

        ImageGrid<coordinate_types::mosaic> rgn1(MosaicBoxRegion(0,0,10,10),5);

        FitsHeader fht;

        rgn1.write_to_fht(fht);
        ImageGrid<coordinate_types::mosaic> rgn2;
        rgn2.load_from_fht(fht);

        CHECK(rgn1.image_outline_==rgn2.image_outline_);
        CHECK(rgn1.pixel_width_==rgn2.pixel_width_);
    }

    TEST(pixel_to_CcdPosition_conversion) {
        ImageGrid<coordinate_types::ccd> trivial(CcdBoxRegion(0,0,10,10), 1.0);
        PixelPosition p1(10,10);

        CcdPosition c1= trivial.corresponding_grid_Position(p1);

        CHECK_EQUAL(p1.x, c1.x);
        CHECK_EQUAL(p1.y, c1.y);

        ImageGrid<coordinate_types::ccd> half_scale(CcdBoxRegion(0,0,10,10), 0.5);

        CcdPosition c2 = half_scale.corresponding_grid_Position(p1);
        CHECK_CLOSE(c2.x, p1.x/2.0, 1e-7);
        CHECK_CLOSE(c2.y, p1.y/2.0, 1e-7);

        ImageGrid<coordinate_types::ccd> half_scale_and_offset(CcdBoxRegion(2,2,10,10), 0.5);
        CcdPosition c3 = half_scale_and_offset.corresponding_grid_Position(p1);
        CHECK_CLOSE(c3.x, p1.x/2.0 +2.0, 1e-7);
        CHECK_CLOSE(c3.y, p1.y/2.0 + 2.0, 1e-7);

    }


}

