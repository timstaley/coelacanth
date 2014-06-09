
#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include "../lucky_math_funcs.h"
#include <iostream>
#include <boost/math/distributions/poisson.hpp>
#include <limits>
#include <boost/math/distributions/hypergeometric.hpp>

using namespace coela;
using namespace std;

SUITE(lucky_math)
{
    using namespace lucky_math;

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"lucky_math\" unit tests running ***" <<endl;
    }

    TEST(integer_factorial) {
        CHECK_EQUAL(1ul, integer_factorial(0));
        CHECK_EQUAL(1ul, integer_factorial(1));
        CHECK_EQUAL(2ul, integer_factorial(2));
        CHECK_EQUAL(6ul, integer_factorial(3));
        CHECK_EQUAL(24ul, integer_factorial(4));

//        cout<<"Unsigned long max: " << numeric_limits<unsigned long>::max() <<endl;
//
//        double max_int_fac =1;
//        while (approx_factorial_double(max_int_fac) < numeric_limits<unsigned long>::max())
//            max_int_fac++;
//
//        cout<<"Implies max input for integer_factorial:"<<max_int_fac<<endl;

    }

    TEST(boost_poisson_dist_function) {
        boost::math::poisson_distribution<> p(0.2);
        CHECK_CLOSE(0.81873 , boost::math::pdf(p, 0), 0.001);
        CHECK_CLOSE(0.16375 , boost::math::pdf(p, 1), 0.001);
        CHECK_CLOSE(0.01637 , boost::math::pdf(p, 2), 0.001);
    }
}