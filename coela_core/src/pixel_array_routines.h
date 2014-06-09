/*
 * File:   pixel_array_routines.h
 * Author: ts337
 *
 * Created on 26 August 2011, 17:06
 */

#ifndef COELA_PIXEL_ARRAY_ROUTINES_H
#define COELA_PIXEL_ARRAY_ROUTINES_H

#include "pixel_array2d.h"
#include "cartesian_coords.h"

namespace coela {
namespace pixel_array_routines {

template<typename T>
PixelPosition centroid(const PixelArray2d<T>& , const PixelRange&);

template<typename T>
PixelPosition thresholded_centroid(const PixelArray2d<T>& , const PixelRange&,
                                   const T pixel_threshold);

template<typename T>
void set_col_to_value(PixelArray2d<T>& img, const int col_index, const T value);

template<typename T>
vector<T> get_region_pixels_vector(const PixelArray2d<T>&, const PixelRange&);

template<typename T>
vector<T> get_region_unmasked_pixels_vector(const PixelArray2d<T>& img,
        const PixelArray2d<T>& img_mask,
        const PixelRange&);


template<typename T>
bool pixel_is_greater_than_neighbours(const PixelArray2d<T> & img,
                                      const PixelIndex& centre_pixel,
                                      const unsigned int neighbour_depth=1);

template<typename T>
std::vector<PixelIndex> find_local_maxima_above_threshold(
    const PixelArray2d<T> & image,
    const PixelArray2d<T> & bad_pixels_mask,
    const PixelRange& search_region,
    const unsigned int neighbour_depth,
    const double& threshold);

///Returns pair<clipped_mean, clipped_std_dev>
template<typename T>
std::pair<double,double> calculate_clipped_mean_and_standard_deviation(
    const PixelArray2d<T>&);

///Rapid pixel binning routine. Note, input bitmap dimensions must be divisible by N
PixelArray2d<double> bin_pixel_array(const PixelArray2d<double>& input,
                                     const int bin_pixel_width);

}//end namespace
}//end namespace

#endif  /* PIXEL_ARRAY_ROUTINES_H */

