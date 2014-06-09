#include <UnitTest++/UnitTest++.h>
#include "../convolution.h"
#include "../image_utils.h"
#include "../pixel_array_routines.h"

#include <iostream>
#include <sstream>

#include <cmath>
using namespace coela;
using namespace std;
SUITE(spatial_convolution)
{

    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"spatial_convolution\" unit tests running ***" <<endl;
    }

    TEST(Convolution_Preserves_Flux) {
        PixelArray2d<double> input(20,20, 0);
        input(10,10)=25;
        input(11,11)=25;

        PixelArray2d<double> kernel(3,3, 1.0/9.0);
        PixelIndex kernel_centre(2,2);

        CHECK_CLOSE(1.0, kernel.sum(), 1e-6);

        PixelRange valid_conv_region_calc =
            spatial_convolution::detail::calculate_valid_region_for_convolution(
                input.range(), kernel.range(), kernel_centre);

        PixelRange valid_manual(2,2, 19,19);

        CHECK(valid_manual==valid_conv_region_calc);

        CHECK_EQUAL(input.region_sum(valid_conv_region_calc), input.sum());

        PixelArray2d<double> conv =
            spatial_convolution::convolve_with_kernel(input,
                    kernel, kernel_centre);

        CHECK(input.range() == conv.range());

        CHECK_CLOSE(1.0, kernel.sum(), 1e-6);

        CHECK_CLOSE(input.sum(), conv.sum(), 1e-6);

//        input.write_to_file("conv_input.fits");
//        conv.write_to_file("conv_output.fits");
    }

    TEST(Convolution_with_symmetric_kernel_Preserves_Centroid) {
        PixelArray2d<double> input(20,20, 0);
        input(10,10)=25;
        input(11,11)=25;

        PixelArray2d<double> kernel(3,3, 1.0/9.0);
        PixelIndex kernel_centre(2,2);

        PixelArray2d<double> conv =
            spatial_convolution::convolve_with_kernel(input,
                    kernel, kernel_centre);

        PixelPosition orig_centroid = pixel_array_routines::centroid(input, input.range());
        PixelPosition conv_centroid = pixel_array_routines::centroid(conv, conv.range());

        CHECK_CLOSE(orig_centroid.x, conv_centroid.x, 1e-6);
        CHECK_CLOSE(orig_centroid.y, conv_centroid.y, 1e-6);
    }

    TEST(convolution_above_threshold) {
        PixelArray2d<double> input(20,20, 0);
        input(5,5)=25;
        input(15, 15)=10;

        PixelArray2d<double> kernel(3,3, 1.0/9.0);
        PixelIndex kernel_centre(2,2);

        PixelArray2d<double> thresh_conv =
            spatial_convolution::convolve_with_kernel_at_pixels_above_threshold(
                input, kernel, kernel_centre, 15.0);

        CHECK(pixel_array_routines::centroid(thresh_conv, thresh_conv.range())
              == PixelPosition(4.5,4.5));
    }

}