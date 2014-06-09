#include "../image_cleanup.h"
#include "coela_utility/src/microstats.h"
#include "coela_core/src/image_utils.h"
#include "coela_core/src/pixel_array_routines.h"
#include <stdexcept>
#include <iostream>
using namespace std;
namespace coela {
namespace image_cleanup {
//=====================================================================================================================

void col_10p_debias_preserving_bg_level(PixelArray2d<double>* bmp_ptr, size_t y_min)
{
    PixelArray2d<double>& bmp = *bmp_ptr;
    PixelArray2d<double> bias_frame(bmp.range().x_dim(), bmp.range().y_dim(), 0.0);

//            Original, short code:
    vector<double> col_vals;
    col_vals.reserve(bmp.range().y_dim());

    for (size_t x(1); x<=bmp.range().x_dim(); ++x) {
        col_vals.clear(); //actually unnecessary since we could just write over all the vals....
        for (size_t y(y_min); y<=bmp.range().y_dim(); ++y) {
            col_vals.push_back(bmp(x,y));
        }
        pixel_array_routines::set_col_to_value(bias_frame, x, vector_tenth_percentile(col_vals));
    }

    double min_col_10p = bias_frame.min_val();
    //assumes a positive background. if it's a negative bg then we probably do want to raise the whole frame?
    if (min_col_10p > 0) { bias_frame -= min_col_10p; }
    bmp -= bias_frame;
//    bmp.add_comment("Column 10p debiased (background preserving)");
    return;
}

template<typename T>
void row_10p_debias_preserving_bg_level(PixelArray2d<T>* img_ptr)
{
    PixelArray2d<T>& img=*img_ptr;
    PixelArray2d<T> bias_frame(img);
    vector<T> row_vals; row_vals.resize(img.range().x_dim());

    for (int y_ind=1; y_ind<=(int)img.range().y_dim(); ++y_ind) {

        for (int x_ind=1; x_ind<=(int)img.range().x_dim(); ++x_ind) {
            row_vals[x_ind-1] = img(x_ind,y_ind);
        }

        double row_10p=vector_tenth_percentile(row_vals);
        for (int x_ind=1; x_ind<=(int)img.range().x_dim(); ++x_ind) {
            bias_frame(x_ind,y_ind) = row_10p;
        }
    }
    double min_row_10p = bias_frame.min_val();
    if (min_row_10p > 0) { bias_frame -= min_row_10p; }
    img-= bias_frame;
}
template void row_10p_debias_preserving_bg_level(PixelArray2d<float>* img_ptr);
template void row_10p_debias_preserving_bg_level(PixelArray2d<double>* img_ptr);

template<class HistogramType>
double determine_bias_pedestal_from_box_in_raw_image(const PixelArray2d<float>& bmp,
        const PixelRange& rgn,
        HistogramType& hist
                                                    )
{

    if (hist.pixel_count()) { hist.clear(); }
    hist.set_min_value(-300); //TO DO : fix hard coding.

    for (PixelIterator pix(rgn); pix!=pix.end; ++pix) {
        hist.count_int_value(bmp(pix));  //NB ASSUMES INTEGER VALUED BMP
        //as in raw data (otherwise truncation will cause a value spike at 0)
        //(faster to cast than use rint, apparently)
    }
    return determine_histogram_bias_pedestal_via_thresholded_centroid(hist);
}

template double determine_bias_pedestal_from_box_in_raw_image(const PixelArray2d<float>&
        bmp,
        const PixelRange& rgn,
        HistogramContainer10bit& hist);
template double determine_bias_pedestal_from_box_in_raw_image(const PixelArray2d<float>&
        bmp,
        const PixelRange& rgn,
        HistogramContainer14bit& hist);



template<class HistogramType>
double determine_histogram_bias_pedestal_via_thresholded_centroid(
    const HistogramType& hist,
    const double threshold_factor)
{

    long max_count;
//    if (skip_zero) max_count= hist.get_count(hist.mode_ignoring_zero());
//    else
    max_count = hist.get_count(hist.mode());

    long count_threshold=max_count*threshold_factor;

    double total_counts=0;
    double weighted_sum_val=0;

    for (size_t i=0; i!=hist.size(); ++i) {
        if (hist[i] > count_threshold
//                && (!skip_zero || hist.get_value(i))
           ) {
            total_counts += hist[i];
            weighted_sum_val+= hist[i] * hist.get_value(i);
        }
    }
    return weighted_sum_val / total_counts;
}

template double determine_histogram_bias_pedestal_via_thresholded_centroid(
    const HistogramContainer10bit& hist,
    const double threshold_factor);
template double determine_histogram_bias_pedestal_via_thresholded_centroid(
    const HistogramContainer14bit& hist,
    const double threshold_factor);

CCDImage<float> get_CCD_default_weight_map(const CCD_calibration_info& ccd_inf)
{
    CCDImage<float> default_weight_map;

    default_weight_map.pix = PixelArray2d<float>(ccd_inf.cropped_PixelRange.x_dim(),
                             ccd_inf.cropped_PixelRange.y_dim(), 1.0);

    default_weight_map.initialize_CCD_grid_to_specific_region(ccd_inf.crop_region);

    vector<CCD_BoxRegion> bad_regions = ccd_inf.get_bad_detector_regions();

    for (size_t rgn_num=0; rgn_num!=bad_regions.size(); ++rgn_num) {
        PixelRange bad_box =
            default_weight_map.CCD_grid.corresponding_pixel_region(
                bad_regions[rgn_num]).bounded_pixels();
        bad_box = PixelRange::overlap(bad_box, default_weight_map.pix.range());
        for (PixelIterator pix(bad_box); pix!=pix.end; ++pix) {
            default_weight_map.pix(pix)=0.0;
        }
    }
    return default_weight_map;
}

template<typename T>
void row_median_debias(PixelArray2d<T>* img_ptr)
{
    PixelArray2d<T>& img=*img_ptr;
    PixelArray2d<T> bias_frame(img);
    vector<T> row_vals; row_vals.resize(img.range().x_dim());

    for (int y_ind=1; y_ind<=(int)img.range().y_dim(); ++y_ind) {

        for (int x_ind=1; x_ind<=(int)img.range().x_dim(); ++x_ind) {
            row_vals[x_ind-1] = img(x_ind,y_ind);
        }

        double row_median=vector_median(row_vals);
        for (int x_ind=1; x_ind<=(int)img.range().x_dim(); ++x_ind) {
            bias_frame(x_ind,y_ind) = row_median;
        }

    }
    img-= bias_frame;

}

template void row_median_debias(PixelArray2d<float>* img_ptr);
template void row_median_debias(PixelArray2d<double>* img_ptr);


//float_bitmap& row_10p_debias(float_bitmap& bmp)        {
//   float_bitmap bias_frame((fits_header)bmp);
////            bias_frame.clear();
//    vector<float> row_vals;
//    row_vals.reserve(bmp.x_dim());
//    float row_sky_level;
//    for (size_t y(1); y<=bmp.y_dim(); ++y)
//    {
//        row_vals.clear();
//        for (size_t x(1); x<=bmp.x_dim(); ++x)
//        {
//            row_vals.push_back(bmp(x,y));
//        }
//        row_sky_level=vector_tenth_percentile(row_vals);
//        bias_frame.set_row_to_value(y, (row_sky_level) );
//    }
//
//    bmp-=bias_frame;
//    bmp.add_comment("Row 10p debiased");
//    return bmp;
//}



//=====================================================================================================================
}//end namespace coela::image_cleanup
}//end namespace coela
