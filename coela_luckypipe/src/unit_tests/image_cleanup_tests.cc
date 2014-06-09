#include <UnitTest++/UnitTest++.h>
#include "../image_cleanup.h"
#include "coela_random/src/random.h"
#include <boost/random/normal_distribution.hpp>
#include <iostream>
#include <sstream>


using namespace coela;
using namespace std;
SUITE(image_cleanup)
{

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"image_cleanup\" unit tests running ***" <<endl;
    }

    TEST(Init_Unuran_Seed) {
        //Need to do this within a function, hence declaration within a TEST macro.
        std::vector < unsigned long > seed;

        seed.push_back(111);
        seed.push_back(222);
        seed.push_back(333);
        seed.push_back(444);
        seed.push_back(555);
        seed.push_back(666);

        //NB changing the seed may cause some tests to fail due to overly tight / crude statistical checks.

        if (! unuran::StreamWrapper::unuran_package_has_been_seeded()) {
            unuran::StreamWrapper::set_unuran_package_seed(seed);
        }
    }

    TEST(bias_pedestal_for_normally_distributed_pixel_image) {
        PixelArray2d<float> img(1000,1000,0.0);

        unuran::StreamWrapper rv;

        double true_bias_pedestal = 750.8;
        double var = 15;

        unuran::GaussianRandomVariate gaussian_dist(true_bias_pedestal, var, rv);

        int n_repeats=1;
        for (int run=0; run!=n_repeats; ++run) {
            for (PixelIterator i(img.range()); i!=i.end; ++i) {
                img(i) = rint(gaussian_dist());
            }

            HistogramContainer14bit hist;
            double est_bias_pedestal =
                image_cleanup::determine_bias_pedestal_from_box_in_raw_image(img,
                        img.range(), hist);

            CHECK_CLOSE(true_bias_pedestal, est_bias_pedestal, 0.5);
        }


    }


}

