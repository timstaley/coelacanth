#include <UnitTest++/UnitTest++.h>
#include "../pixel_index.h"

#include <iostream>
#include <sstream>

using namespace coela;
using namespace std;


SUITE(PixelIndex)
{
    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"PixelIndex\" unit tests running ***" <<endl;
    }


    TEST(null_constructor) {
        PixelIndex null_index;
        CHECK(null_index.is_valid()==false);
    }

    TEST(valid_constructor) {
        PixelIndex valid_index(5,-10);
        CHECK(valid_index.is_valid()==true);
    }

}

SUITE(PixelRange)
{
    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"pixel_box\" unit tests running ***" <<endl;
    }


    TEST(null_constructor) {
        PixelRange null_box;
        CHECK(null_box.is_valid()==false);
    }


    TEST(valid_constructor) {
        PixelRange valid_box(5,6,8,10);
        CHECK(valid_box.is_valid()==true);

        PixelRange bad_box(5,6,-10,10);
        CHECK(bad_box.is_valid()==false);
    }

    TEST(PixelIndex_retrieval_for_linear_vector_index) {
        PixelRange image_outline_(1,1,10,10);

        CHECK(PixelIndex(1,1) == image_outline_.get_PixelIndex_for_data_vector_element(0));
        CHECK(PixelIndex(2,1) == image_outline_.get_PixelIndex_for_data_vector_element(1));
        CHECK(PixelIndex(1,2) == image_outline_.get_PixelIndex_for_data_vector_element(
                  image_outline_.x_dim()));
        CHECK(PixelIndex(10,10) == image_outline_.get_PixelIndex_for_data_vector_element(99));

    }

}

