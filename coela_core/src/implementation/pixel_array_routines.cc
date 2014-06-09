/*
 * File:   pixel_array_routines.cc
 * Author: ts337
 *
 * Created on 26 August 2011, 17:06
 */

#include "../pixel_array_routines.h"
#include "coela_utility/src/microstats.h"
#include <stdexcept>
using std::pair;

namespace coela {
namespace pixel_array_routines {

template<typename T>
PixelPosition centroid(const PixelArray2d<T>& input, const PixelRange& rgn)
{
    PixelShift first_moment(0.0,0.0);
    double sum=0.0;

    for (PixelIterator i(rgn); i!=i.end; ++i) {
        double pixval = input(i);
        first_moment += (PixelPosition::centre_of_pixel(i) - PixelPosition::origin) * pixval;
        sum+=pixval;
    }

    first_moment/=sum;
    return first_moment + PixelPosition::origin;
}

template PixelPosition centroid(const PixelArray2d<double>& input, const PixelRange& rgn);
template PixelPosition centroid(const PixelArray2d<float>& input, const PixelRange& rgn);
//-----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
PixelPosition thresholded_centroid(const PixelArray2d<T>& input, const PixelRange& rgn,
                                   const T pixel_threshold)
{
    PixelShift first_moment(0.0,0.0);
    double sum=0.0;

    for (PixelIterator i(rgn); i!=i.end; ++i) {
        double pixval = input(i);
        if (pixval>pixel_threshold) {
            first_moment += (PixelPosition::centre_of_pixel(i) - PixelPosition::origin) * pixval;
            sum+=pixval;
        }
    }

    first_moment/=sum;
    return first_moment + PixelPosition::origin;
}

template PixelPosition thresholded_centroid(const PixelArray2d<double>& input,
        const PixelRange& rgn,
        const double pixel_threshold);
template PixelPosition thresholded_centroid(const PixelArray2d<float>& input,
        const PixelRange& rgn,
        const float pixel_threshold);
//-----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void set_col_to_value(PixelArray2d<T>& img, const int col_index, const T value)
{
    PixelRange range = img.range();
    PixelRange col_rgn(col_index, range.low.y, col_index, range.high.y);
    for (PixelIterator i(col_rgn); i!=i.end; ++i) {
        img(i)= value;
    }
}

template void set_col_to_value(PixelArray2d<double>& img, const int col_index,
                               const double value);
template void set_col_to_value(PixelArray2d<float>& img, const int col_index,
                               const float value);
//-----------------------------------------------------------------------------------------------------------------------------------------
template<typename T>
vector<T> get_region_pixels_vector(const PixelArray2d<T>& img, const PixelRange& range)
{
    vector<T> vals; vals.reserve(range.n_pix());
    for (PixelIterator it(range); it!=it.end; ++it) {
        vals.push_back(img(it));
    }
    return vals;
}

template vector<double> get_region_pixels_vector(const PixelArray2d<double>&,
        const PixelRange&);
//-----------------------------------------------------------------------------------------------------------------------------------------
template<typename T>
vector<T> get_region_unmasked_pixels_vector(const PixelArray2d<T>& img,
        const PixelArray2d<T>& img_mask,
        const PixelRange& range)
{
    vector<T> vals; vals.reserve(range.n_pix());
    T zero_val(0); //<-- not sure if this is necessary... unit test (to do!)
    for (PixelIterator it(range); it!=it.end; ++it) {
        if (img_mask(it)!=zero_val) { vals.push_back(img(it)); }
    }
    return vals;
}
template vector<double> get_region_unmasked_pixels_vector(const PixelArray2d<double>&,
        const PixelArray2d<double>&, const PixelRange&);
//-----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
bool pixel_is_greater_than_neighbours(const PixelArray2d<T> & img,
                                      const PixelIndex& centre_pixel,
                                      const unsigned int neighbour_depth)
{
    PixelRange scan_box(centre_pixel,centre_pixel);
    PixelIndex diagonal(neighbour_depth,
                        neighbour_depth); //Could change this for larger box size;

    scan_box.low -= diagonal;
    scan_box.high +=diagonal;

    scan_box = PixelRange::overlap(img.range(), scan_box);
    T centre_val = img(centre_pixel);
    for (PixelIterator it(scan_box); it!=it.end; ++it) {
        if (it!=centre_pixel) {
            if (img(it) > centre_val) { return false; }
        }
    }
    return true;
}

template
bool pixel_is_greater_than_neighbours(const PixelArray2d<double> & img,
                                      const PixelIndex& centre_pixel,
                                      const unsigned int neighbour_depth);

//-----------------------------------------------------------------------------------------------------------------------------------------

//For background level and variance estimation
//Sigma clipping eliminates the large variances due to stars, hot pixels, etc.
//Use an equation similar to that used by S-Extractor
template<typename T>
pair<double,double> calculate_clipped_mean_and_standard_deviation(
    const PixelArray2d<T> & img)
{
    vector<T> img_vals = get_region_pixels_vector(img, img.range());

    double pseudo_mode =    2.5*vector_median(img_vals) - 1.5*vector_mean(img_vals);
    double original_sd=vector_std_dev(img_vals);
    vector<T> clipped_vals=clip_by_limit_about_value(img_vals, original_sd*2.5, pseudo_mode);

    double clipped_mean = vector_mean(clipped_vals);
    double clipped_std_dev=vector_std_dev(clipped_vals);

    return pair<double,double>(clipped_mean, clipped_std_dev);
}

template
pair<double,double> calculate_clipped_mean_and_standard_deviation(
    const PixelArray2d<double>&);

//-----------------------------------------------------------------------------------------------------------------------------------------
template<typename T>
std::vector<PixelIndex> find_local_maxima_above_threshold(
    const PixelArray2d<T> & img,
    const PixelArray2d<T> & bad_pixels_mask,
    const PixelRange& search_region,
    const unsigned int neighbour_depth,
    const double& threshold)
{


    PixelRange bordered=img.range();
    PixelIndex one_diag(1,1);
    bordered.low += one_diag;
    bordered.high -= one_diag;

    vector<PixelIndex> local_maxima;

    PixelRange final_search_region = PixelRange::overlap(bordered,
                                     search_region);

    T zero_val(0);
    for (PixelIterator pix(final_search_region); pix!=pix.end; ++pix) {
        if (img(pix)  > threshold) {
            if (bad_pixels_mask(pix)!=zero_val &&
                    pixel_is_greater_than_neighbours(img, pix, neighbour_depth)) {
                local_maxima.push_back(pix);
//                cerr<<"Found maxima at " <<pix<<"; val: " <<img(pix)<<endl;
            }
        }
    }
    return local_maxima;
}

template std::vector<PixelIndex> find_local_maxima_above_threshold(
    const PixelArray2d<double> & img,
    const PixelArray2d<double> & bad_pixels_mask,
    const PixelRange& search_region,
    const unsigned int neighbour_depth,
    const double& thresh);
//---------------------------------------------------------------------------------

PixelArray2d<double> bin_pixel_array(const PixelArray2d<double>& input_image,
                                     const int bin_pixel_width)
{
    if (input_image.range().x_dim() % bin_pixel_width !=0 ||
            input_image.range().y_dim() % bin_pixel_width !=0) {
        throw std::runtime_error("bin_image: Cannot bin this array, dimensions not divisible by bin width");
    }

    PixelArray2d<double> output;
    output = PixelArray2d<double> (input_image.range().x_dim()/bin_pixel_width,
                                   input_image.range().y_dim()/bin_pixel_width,
                                   0.0);



    //TO DO : can (perhaps) optimize this better - would need to profile.
    // Could calculate output row per input row.
    // Then 1d nested loop, incrementing output per bin_size inputs.
    // (That's how it's done in the old code)

    for (PixelIterator out(output.range()); out!=out.end; ++out) {
        PixelRange input_range((out.x-1)*bin_pixel_width +1, (out.y-1)*bin_pixel_width +1,
                               out.x*bin_pixel_width , out.y*bin_pixel_width);
        for (PixelIterator in(input_range); in!=in.end; ++in) {
            output(out) += input_image(in);
        }
    }

    return output;
}

}//end namespace
}//end namespace