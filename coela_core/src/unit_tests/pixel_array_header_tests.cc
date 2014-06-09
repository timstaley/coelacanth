#include <UnitTest++/UnitTest++.h>
#include "../PixelArrayHeader.h"



#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <exception>

using namespace coela;
using namespace std;


extern string lucky_lib_test_resources_dir;

SUITE(PixelArrayHeader)
{

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"PixelArrayHeader\" unit tests running ***" <<endl;
    }


    TEST(header_size_calculation_correct) {

        PixelArray2dHeader tiny_float_hdr(fits_header_conventions::FLOATIMG,1,1);
        PixelArray2dHeader float_hdr(fits_header_conventions::FLOATIMG,10,5);
        CHECK_EQUAL(4, tiny_float_hdr.byte_size_of_uncompressed_data());
        CHECK_EQUAL(50*4, float_hdr.byte_size_of_uncompressed_data());

    }

}
