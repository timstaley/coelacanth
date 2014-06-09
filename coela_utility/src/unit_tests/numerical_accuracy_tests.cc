
#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include <iostream>

#include <limits>

//using namespace coela;
using namespace std;

SUITE(numerical_accuracy)
{


    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \""<<UnitTest::CurrentTest::Details()->suiteName
             <<"\" unit tests running ***" <<endl;
    }

    TEST(float_accuracy_limit) {

        //A float is accurate to about 1 in ten million
        //But converting to double, then adding floats will preserve the accuracy of the floats

        float hundred_million=1e8;
        float float_sum = hundred_million;
        float_sum += 1.0f;
        CHECK_EQUAL(float_sum, hundred_million);

        double double_sum = hundred_million;
        double_sum +=1.0f;
        CHECK(double_sum!=hundred_million);


        float ten_million=1e7;
        float_sum = ten_million;
        float_sum+=1.0f;
        CHECK(float_sum != ten_million);


        //If the above example seems contrived, consider the following:
        float_sum = 0.0f;
        double_sum=0.0;
        for (size_t i=0; i!=1e8; ++i) {
            float_sum+=1.0f;
            double_sum+=1.0f;
        }

//        cerr<<"float_sum:"<<float_sum<<endl;
//        cerr<<"double_sum:"<<double_sum<<endl;

        CHECK(float_sum<5e7);
        CHECK_EQUAL(double_sum, 1e8);

        float down_conversion(double_sum);
        CHECK_EQUAL(down_conversion, 1e8);

        //So, consider we sum all the pixels of a 1k squared frame, pixel by pixel for each frame,
        //  for 10,000 frames.
        // --- Not totally implausible
        // --- Definitely need double precision for sums.
        // --- Can down-convert e.g. a bias frame once we have used double precision for the initial creation.

    }
}