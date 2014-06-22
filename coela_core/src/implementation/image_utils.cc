#include "../image_utils.h"
#include "coela_utility/src/microstats.h"

#include "../pixel_array_routines.h"
//#include <deque>
#include <iostream>
#include <stdexcept>
#include <algorithm>
using namespace std;

namespace coela {
namespace image_utils {
//=====================================================================================================================


//-----------------------------------------------------------------------------------------------------------------------------------------
//Indexing from original dataset to resampled dataset requires an offset due to the fencepost problem
// i.e.  N points at original sampling, (N-1)*resample_factor + 1 points in resampled set.
//However when coding it's much easier to map from input pixel index -> output PixelIndex if we can just multiply by the resample factor

//For this reason we add dummy columns for 1...R-1 where R is the interpolation factor. hence (input index x) -> (output index Rx) .
//Hence the *input sample co-ordinate Position* of a resampled pixel is equal to     (output_index/R - 0.5); or  (output_index - R/2) / R
//and if we are working in the co-ordinate numberline of the resampled pixel, then we must add 0.5 to move to the index numberline.


//notes:
//Since the inter-spacing is always 0.5 (,0.25, 0.75) pix we can pre-calculate
//the coefficients for each of the 4 sample points (a-1, a0, a1, a2) per 1 dimensional resampling.
//See wikipedia article 'Bicubic Interpolation' for maths theory behind this
//implementation (originally R. Keys 1981).
//NB! Edges are messy - would need to implement a half-matrix jobby to fix this
//- or linearly interpolate, perhaps

///Interpolates for a mid-way point:
double cubic_resample_mid_point(const vector<double>& v)
{
    return ((-1.0*v[0]  +9.0*v[1] + 9.0*v[2] -1.0*v[3]) /16.0);
}
double cubic_resample_first_quarter_point(const vector<double>& v)
{
    return (-7.0*v[0] + 105.0*v[1] + 35.0*v[2] -5.0*v[3])/128.0 ;
}
double cubic_resample_third_quarter_point(const vector<double>& v)
{
    return (-7.0*v[3] + 105.0*v[2] + 35.0*v[1] -5.0*v[0])/128.0 ;
}

CcdImage<double> initialise_bmp_for_resampling(
    const CcdImage<double>& input_bmp,
    const double pixscale_factor)
{

    /*
     The natural number of pixels in an image of width X, resampled by a factor
     of N is:
          (x-1)*N +1
     This is due to the fencepost problem. (sx samples, sx-1 interpolated sections)

    However to keep indexing simple, the output is generated with x*N pixels
    and the first  (n-1) are left blank.
    We then tweak the CCD reference co-ordinate grid to correctly represent
    this shift.

    in comments below, N == zoom_factor (e.g. 2,4,8)

    NB  the first (N-1) column(s) will always be a zero since they correspond
    to a point left of original grid (dummy cols).
    This is done to keep indexing math simple.
    So the (N)th (e.g. 2nd) pixel corresponds to the sample point originally
    at Cartesian position 0.5 (i.e. sample position represented by first input
    pixel).
    Samples are then spaced at pixscale_factor increments.
    --So the first pixel is centred at (Cartesian position) -0.5 + pixscale_factor
    --So the lower_bound is at ( -0.5 + pixscale_factor/2.0)
    (in the Cartesian co-ordinate system of the frame we are resampling).
    */

    int zoom_factor = 1/pixscale_factor;
    CcdImage<double> output_bmp;
    output_bmp.pix = PixelArray2d<double>(zoom_factor*input_bmp.pix.range().x_dim(),
                                          zoom_factor*input_bmp.pix.range().y_dim(),
                                          0.0);
    //So, calculate the lower-left cartesian co-ordinate of this new pixel array,
    // in the reference frame of the input image:
    PixelPosition lower_bound(-0.5 +pixscale_factor/2.0, -0.5 +pixscale_factor/2.0);
    //Convert to the CCD reference frame:
    CcdPosition output_CCD_region_corner(
        input_bmp.CCD_grid.corresponding_grid_Position(lower_bound));

    //And initialise accordingly
    output_bmp.initialize_CCD_grid_to_specific_offset_and_scale(
        output_CCD_region_corner,
        input_bmp.CCD_grid.pixel_width_*pixscale_factor);

    return output_bmp;
}

//        float_bitmap initialise_bmp_for_resampling(const float_bitmap& input_bmp, const double pixscale_factor, const pixel_box& input_rgn, const int input_padding){
//            PixelIndex lowest_padding_pixel(input_rgn.low.x - input_padding, input_rgn.low.y - input_padding);
//            PixelPosition low_pixel_centre PixelPosition::pixel_centre(lowest_padding_pixel);
//            PixelPosition low_pixel_bount
//            PixelPosition lower_bound()
//        }

CcdImage<double> bicubic_resample_2x(const CcdImage<double>& input_bmp)
{
    const size_t in_sx(input_bmp.pix.range().x_dim()), in_sy(input_bmp.pix.range().y_dim());
    const double pixscale_factor = 0.5;
    size_t resample_factor = 1.0/pixscale_factor;
    CcdImage<double> output_bmp(initialise_bmp_for_resampling(input_bmp, pixscale_factor));

    //Already sampled:
    for (size_t y_in = 1, y_out=resample_factor; y_in<=in_sy;
            ++y_in, y_out+=resample_factor) {
        for (size_t x_in = 1, x_out=resample_factor; x_in<=in_sx;
                ++x_in, x_out+=resample_factor) {
            output_bmp.pix(x_out,y_out) = input_bmp.pix(x_in, y_in);
        }
    }

    //First interpolate along edges where there is room to do standard interpolation in one direction but not the other.
    //top and bottom first
    //Fixme: use deques? Would have to do a comparative test to see if it's really worth it.
    //FIXME high edge minus one co-ords get missed unless we use slow 'ifs' in main loop... to be put here later...(perhaps)
    vector<double> low_samples(4), high_samples(4); //, high_m1_samples(4);
    for (size_t x=2; x<in_sx - 1; ++x) {
        //here the value of x denotes a0
        for (size_t i(0); i<4; ++i) {
            low_samples[i]=input_bmp.pix(x -1 + i , 1);
            high_samples[i]=input_bmp.pix(x -1 + i , in_sy);
        }
        output_bmp.pix(x*2 +1, 2)= cubic_resample_mid_point(low_samples) ;
        output_bmp.pix(x*2 +1, in_sy*2) =cubic_resample_mid_point(high_samples) ;
    }
    //now left and right
    for (size_t y=2; y<in_sy - 1; ++y) {
        for (size_t i(0); i<4; ++i) {
            low_samples[i]=input_bmp.pix(1 , y -1 + i);
            high_samples[i]=input_bmp.pix(in_sx , y -1 + i);
        }
        output_bmp.pix(2, y*2 +1)= cubic_resample_mid_point(low_samples) ;
        output_bmp.pix(in_sx*2, y*2 +1)= cubic_resample_mid_point(high_samples) ;
    }

    //Now interpolate all points away from edges that are on connecting lines between input points (or are input points themselves)
    vector<double> x_samples(4), y_samples(4);
    for (size_t y =2; y<= in_sy -1; ++y) {
        //Odd co-ords are the interpolated points. (row / col 1 is simply zero, for padding)
        for (size_t x =2; x<=in_sx - 1; ++x) {
            //Already sampled:
//                        output_bmp.pix( x*2,y*2)= input_bmp.pix( x, y ) ;

            //Interpolation along grid lines in x and y directions: (NB - this is at the next half point - don't do it in direction of boundary if we're at the final non-edge sample)
            //Fixme: Remove this if statement in favour of an extra loop along the x=in_sx-1 and y=in_sy-1 lines.
            if (x<in_sx -1) {
                for (size_t i(0); i<4; ++i) {
                    x_samples[i]=input_bmp.pix(x -1 + i , y);    //again x denotes a0
                }
                output_bmp.pix(x*2+1 , y*2)= cubic_resample_mid_point(x_samples) ;
            }
            if (y<in_sy -1) {
                for (size_t i(0); i<4; ++i) {
                    y_samples[i]=input_bmp.pix(x   , y -1 + i);
                }
                output_bmp.pix(x*2 , y*2 +1)= cubic_resample_mid_point(y_samples) ;
            }
        } //end input y-scan
    }//end input x-scan.

    //and finally mid-points away from edges (this relies on having interpolated the points on connecting lines, first
    //-i.e. this is the BI-cubic bit. Hence we get values from 'output.'  ) We re-interpolate in the x direction as this improves cache performance
    size_t orig_x,oy; //(out_y)
    for (size_t y =2; y< in_sy -1; ++y) {
        for (size_t x =2; x<in_sx - 1; ++x) {
            orig_x = 2*x;
            oy= 2*y +1;
            x_samples[0]=output_bmp.pix(orig_x -2 , oy);
            x_samples[1]=output_bmp.pix(orig_x  , oy);
            x_samples[2]=output_bmp.pix(orig_x +2 , oy);
            x_samples[3]=output_bmp.pix(orig_x +4 , oy);
            output_bmp.pix(orig_x + 1,oy)= cubic_resample_mid_point(x_samples);
        }
    }

    output_bmp.pix *= pixscale_factor*pixscale_factor; //total flux preserving
    return output_bmp;
}

CcdImage<double> bicubic_resample_4x(const CcdImage<double>& input_bmp)
{
    //for resample factor R, indexing from 1 to N, y(out) = Rx - (R-1) = R(x-1) +1
    const size_t in_sx(input_bmp.pix.range().x_dim()), in_sy(input_bmp.pix.range().y_dim());
    const double pixscale_factor = 0.25;
    CcdImage<double> output_bmp(initialise_bmp_for_resampling(input_bmp, pixscale_factor));

    size_t resample_factor = 1.0/pixscale_factor;
    //Already sampled:
    for (size_t y_in = 1, y_out=resample_factor; y_in<=in_sy;
            ++y_in, y_out+=resample_factor) {
        for (size_t x_in = 1, x_out=resample_factor; x_in<=in_sx;
                ++x_in, x_out+=resample_factor) {
            output_bmp.pix(x_out,y_out) = input_bmp.pix(x_in, y_in);
        }
    }

    //First interpolate along edges of input bitmap; where there is room to do standard interpolation in one direction but not the other. NB if ever fix edges put in pre-sampled edge points here...
    //top and bottom first
    //high or low refers to co-ordinate value in either x or y
    vector<double> low_samples(4), high_samples(4);
    for (size_t x=2; x<in_sx - 1; ++x) {
        for (size_t i(0); i<4; ++i) {
            low_samples[i]=input_bmp.pix(x -1 + i , 1);
            high_samples[i]=input_bmp.pix(x -1 + i , in_sy);
        }
        output_bmp.pix(x*4 +1, 4)= cubic_resample_first_quarter_point(low_samples) ;
        output_bmp.pix(x*4 +2, 4)= cubic_resample_mid_point(low_samples) ;
        output_bmp.pix(x*4 +3, 4)=cubic_resample_third_quarter_point(low_samples) ;
        output_bmp.pix(x*4 +1, in_sy*4)=cubic_resample_first_quarter_point(high_samples) ;
        output_bmp.pix(x*4 +2, in_sy*4)=cubic_resample_mid_point(high_samples) ;
        output_bmp.pix(x*4 +3, in_sy*4)=cubic_resample_third_quarter_point(high_samples) ;
    }
    //now left and right
    for (size_t y=2; y<in_sy - 1; ++y) {
        for (size_t i(0); i<4; ++i) {
            low_samples[i]=input_bmp.pix(1 , y -1 + i);
            high_samples[i]=input_bmp.pix(in_sx , y -1 + i);
        }
        output_bmp.pix(4, y*4 +1)= cubic_resample_first_quarter_point(low_samples) ;
        output_bmp.pix(4, y*4 +2)= cubic_resample_mid_point(low_samples) ;
        output_bmp.pix(4, y*4 +3)= cubic_resample_third_quarter_point(low_samples);
        output_bmp.pix(in_sx*4, y*4 +1)= cubic_resample_first_quarter_point(high_samples) ;
        output_bmp.pix(in_sx*4, y*4 +2)= cubic_resample_mid_point(high_samples) ;
        output_bmp.pix(in_sx*4, y*4 +3)= cubic_resample_third_quarter_point(high_samples) ;
    }

    //Now interpolate all points away from edges where either ox or oy is even. (connecting lines between pre-sampled points)
    vector<double> x_samples(4), y_samples(4);
    for (size_t y =2; y<= in_sy -1; ++y) {
        for (size_t x =2; x<=in_sx - 1; ++x) {
            //Already sampled:
//                    output_bmp.pix( x*4,y*4)= input_bmp.pix( x, y ) ;

            //Interpolation along grid lines in x and y directions: (NB - this is the next interpolation gap - don't try to do it in direction of boundary
            //if we're at the final non-edge sample, or you'll end up referencing a bad index.)
            if (x<in_sx -1) {
                for (size_t i(0); i<4; ++i) {
                    x_samples[i]=input_bmp.pix(x -1 + i , y);
                }
                output_bmp.pix(x*4 +1, y*4)= cubic_resample_first_quarter_point(x_samples) ;
                output_bmp.pix(x*4 +2 , y*4)= cubic_resample_mid_point(x_samples) ;
                output_bmp.pix(x*4 +3 , y*4)= cubic_resample_third_quarter_point(x_samples);
            }
            if (y < in_sy -1) {
                for (size_t i(0); i<4; ++i) {
                    y_samples[i]=input_bmp.pix(x   , y -1 + i);
                }
                output_bmp.pix(x*4 , y*4 +1) =cubic_resample_first_quarter_point(y_samples);
                output_bmp.pix(x*4 , y*4 +2) =cubic_resample_mid_point(y_samples);
                output_bmp.pix(x*4 , y*4 +3) =cubic_resample_third_quarter_point(y_samples);
            }
        } //end input y-scan
    }//end input x-scan.

    //and finally mid-points away from edges - rely on the fact that we have interpolated all points on connecting lines, and increment quarter blocks in y for each input point. This is a very inner loop!
    size_t orig_x,orig_y, oy;
    for (size_t y =2; y< in_sy -1; ++y) {
        for (size_t x =2; x<in_sx - 1; ++x) {
            orig_x = 4*x;
            orig_y= 4*y;

            for (size_t y_quart=1; y_quart!=4; ++y_quart) {
                oy = orig_y + y_quart;
                x_samples[0]=output_bmp.pix(orig_x -4 , oy);
                x_samples[1]=output_bmp.pix(orig_x  , oy);
                x_samples[2]=output_bmp.pix(orig_x +4 , oy);
                x_samples[3]=output_bmp.pix(orig_x +8 , oy);
                output_bmp.pix(orig_x +1 , oy) = cubic_resample_first_quarter_point(x_samples);
                output_bmp.pix(orig_x +2 , oy) = cubic_resample_mid_point(x_samples);
                output_bmp.pix(orig_x +3 , oy) = cubic_resample_third_quarter_point(x_samples);
            }

        }
    }
    output_bmp.pix *= pixscale_factor*pixscale_factor; //flux preserving
    return output_bmp;
}


//These are an attempt to solve the problem of borders getting left blank,
//by using information from outside the selection region to produce the interpolated result
//NB in the same way that the routines above have dummy rows/cols,
//here we actually return pixels that use information from outside the selection region
//but the main thing is that the pixel indices are mapped simply,
//i->j=Ri. (so x = i-0.5 -> j = Ri   so;  j = R(x + 0.5)  so x = j/R - 0.5

template<typename input_datatype>
CcdImage<double> bicubic_resample_non_edge_region(
    const CcdImage<input_datatype>& input_bmp,
    const CcdBoxRegion ccd_roi,    //roi: Region of interest
    const double resample_factor)
{
    PixelRange box =
        (input_bmp.CCD_grid.corresponding_pixel_region(ccd_roi)).bounded_pixels();

    if (box.low.x<2 || box.high.x > ((int)input_bmp.pix.range().x_dim() - 1) || box.low.y<2
            || box.high.y > ((int)input_bmp.pix.range().y_dim() - 1)) {
        throw domain_error("Cannot resample region this close to edge of frame - use whole image routine.");
    }

    if (resample_factor!=2.0 && resample_factor!=4.0) { throw domain_error("Only resample factors implemented are 2.0 and 4.0"); }

    PixelRange padded_box(box.low.x -1, box.low.y-1, box.high.x+1, box.high.y+1);
    CcdImage<double> padded_bmp = CcdImage<double>(CcdImage<input_datatype>::sub_image(
                                      input_bmp, padded_box));
    CcdImage<double> padded_resample;
    if (resample_factor==2.0) { padded_resample = bicubic_resample_2x(padded_bmp); }
    else if (resample_factor==4.0) { padded_resample = bicubic_resample_4x(padded_bmp); }
    else { throw logic_error("Resample factor not supported"); }


    //Now, the CCD_region outline has changed slightly.
    //If we consider the pixel centre as the "sample point" and keep that Position fixed in the CCD_coords
    //then we must shrink the pixel outline around that centre point, to the smaller resampled pixel size.
    //As such, the pixel outline for the whole image shrinks slightly, even though the edge pixels are the same.

    double horizontal_distance_from_corner_of_old_pixel_to_corner_of_shrunk_pixel =
        input_bmp.CCD_grid.pixel_width_/2.0 -
        padded_resample.CCD_grid.pixel_width_/2.0;

    CcdPixelShift diagonal_shrink(
        horizontal_distance_from_corner_of_old_pixel_to_corner_of_shrunk_pixel,
        horizontal_distance_from_corner_of_old_pixel_to_corner_of_shrunk_pixel);

    CcdBoxRegion ccd_roi_with_new_pixel_sizing(ccd_roi.low + diagonal_shrink,
            ccd_roi.high - diagonal_shrink);

    return CcdImage<double>::sub_image(padded_resample, ccd_roi_with_new_pixel_sizing);
}

template
CcdImage<double> bicubic_resample_non_edge_region(const CcdImage<float>& input_bmp,
        const CcdBoxRegion ccd_roi, const double resample_factor);
template
CcdImage<double> bicubic_resample_non_edge_region(const CcdImage<double>& input_bmp,
        const CcdBoxRegion ccd_roi, const double resample_factor);

CcdImage<double> bin_image(const CcdImage<double>& input_image, const int bin_pixel_width)
{
    CcdImage<double> output;
    output.pix = pixel_array_routines::bin_pixel_array(input_image.pix,
                 bin_pixel_width);

    if (input_image.CCD_grid.is_initialized()) {
        output.initialize_CCD_grid_to_specific_offset_and_scale(
            input_image.CCD_grid.image_outline_.low,
            input_image.CCD_grid.pixel_width_ * bin_pixel_width);
    }
    return output;
}
//=====================================================================================================================
}//end of namespace coela::image_utils
}//end of namespace coela




