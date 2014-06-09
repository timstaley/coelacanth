#include "../psf_characterisation.h"

#include <fstream>
#include <algorithm>
#include "coela_utility/src/microstats.h"
#include "coela_utility/src/lucky_math_funcs.h"

#include <iostream>

using namespace std;

namespace coela {
namespace psf_characterisation {
//
//double estimate_fwhm_in_image_pix(const PixelArray2d<double>& bmp, const CCD_Position& CCD_lock, const double est_peak_value, const double max_radius){
//    PixelPosition frame_pos = bmp.corresponding_frame_Position(CCD_lock);
//    return  estimate_fwhm_in_image_pix(bmp, frame_pos, est_peak_value, max_radius);
//
//}
//
//
double estimate_fwhm_in_image_pix(const PixelArray2d<double>& bmp,
                                  const PixelPosition& origin, const double est_peak_value, const double max_radius)
{
    return 2.0*estimate_radius_of_avg_value(bmp, origin, max_radius,  est_peak_value*0.5);
}


double estimate_radius_of_avg_value(const PixelArray2d<double>& bmp,
                                    const PixelPosition& origin, const double max_radius, const double pixel_level)
{
    using string_utils::ftoa;
    using namespace psf_characterisation::detail;
    vector<pair <double, double> > raw_data = get_radial_data(bmp, origin, max_radius);
    vector<psf_profile_point> profile = bin_average_radial_data(raw_data);
    assert(profile.empty()==false);
//             double max_value = profile[0].mean_value;
    size_t radius_index=0;

    while (profile[radius_index].mean_value > pixel_level && radius_index!=profile.size()) {
        // proceed to one past desired level
        ++radius_index;
    }
    if (radius_index==profile.size()) throw logic_error(
            "estimage_radius_of_avg_value: This profile does not descend past desired value, enlarge the radius or check the debiasing\n"
            "(img coords: " + ftoa(origin.x)+","+ftoa(origin.y)+"), value " + ftoa(bmp(
                        PixelPosition::pixel_containing_point(origin)))+"."
        );
    return lucky_math::linearly_interpolate(
               profile[radius_index-1].mean_value, profile[radius_index-1].avg_radius,
               profile[radius_index].mean_value, profile[radius_index].avg_radius, pixel_level);
}
//
double estimate_radius_to_enclose_flux(const PixelArray2d<double>& bmp,
                                       const PixelPosition& origin,
                                       const double flux_sum_to_enclose,
                                       const double max_radius)
{
    using namespace psf_characterisation::detail;
    vector<pair <double, double> > raw_data = get_radial_data(bmp, origin, max_radius);

    vector<psf_profile_point> profile = bin_average_radial_data(raw_data);
    assert(profile.empty()==false);

    size_t radius_index=0;
    double prev_flux_sum=0.0;

    while ((profile[radius_index].mean_value*profile[radius_index].n_samples+ prev_flux_sum)
            <flux_sum_to_enclose
            && radius_index!=profile.size()) {
        prev_flux_sum+=profile[radius_index].mean_value*profile[radius_index].n_samples;
        ++radius_index; // proceed to one past desired level
    }
    if (radius_index==profile.size()) throw logic_error(
            "estimage_radius_of_enclosed_flux: This profile does not descend past desired flux level, enlarge the radius or check the debiasing\n");

    double next_flux_sum = prev_flux_sum +
                           profile[radius_index].mean_value*profile[radius_index].n_samples;

    return lucky_math::linearly_interpolate(
               prev_flux_sum, profile[radius_index-1].avg_radius, next_flux_sum,
               profile[radius_index].avg_radius, flux_sum_to_enclose);

//     //Sort the pixel values from low radius to high radius
//     sort(raw_data.begin(), raw_data.end(), double_pair_first_member_predicate);
//
//     vector<psf_profile_point> profile = bin_average_radial_data(raw_data);
//     assert (profile.empty()==false);
//
//     double prev_sum_flux=0.0;
//     double sum_flux=0.0;
//
//     size_t index;
//     for (index=0; index!=raw_data.size() && sum_flux < (flux_sum_to_enclose - lucky_math::err_margin) ; ++index){
//         prev_sum_flux = sum_flux;
//         sum_flux += raw_data[index].second;
//     }
//
//     if ( index==1 ) return 0.5;
//     if (index==raw_data.size() && sum_flux < (flux_sum_to_enclose - 1e-4) ) {
//         cerr<<"Wanted "<<flux_sum_to_enclose<<" got: "<< sum_flux<<endl;
//         throw logic_error(
//             "estimate_radius_of_enclosed_flux: This profile does not descend past desired flux level, enlarge the radius or check the debiasing\n");
//     }
//
//     cout<<"Previous point: " <<raw_data[index-2].first<< " "<< prev_sum_flux<<endl;
//     cout<<"Above point: " <<raw_data[index-1].first<< " "<< prev_sum_flux<<endl;
//
//     return lucky_math::linearly_interpolate(
//             prev_sum_flux, raw_data[index-2].first, sum_flux, raw_data[index-1].first, flux_sum_to_enclose   );
}


double estimate_FWHEF_in_image_pix(const PixelArray2d<double>& bg_subbed_bmp,
                                   const PixelPosition& origin, const double total_flux, const double max_radius)
{
    return 2.0*estimate_radius_to_enclose_flux(bg_subbed_bmp, origin, total_flux*0.5 ,
            max_radius);
}

double estimate_encircled_flux_at_pixel_radius(const PixelArray2d<double>& bmp,
        const PixelPosition& origin,
        const double flux_aperture_pixel_radius
                                              )
{
    using namespace psf_characterisation::detail;

    double data_gathering_radius = flux_aperture_pixel_radius + 3;

    PixelBoxRegion flux_region(origin, origin);
    flux_region = PixelBoxRegion::padded_copy(flux_region, data_gathering_radius);

    if (!PixelBoxRegion::pixel_box_outline(bmp.range()).contains_region(flux_region)) {
        throw runtime_error("estimate_encircled_flux_at_pixel_radius -- image input is too small to take pixel sums to given radius");
    }

    vector<pair <double, double> > raw_data = get_radial_data(bmp, origin,
            data_gathering_radius);

    vector<psf_profile_point> profile = bin_average_radial_data(raw_data);
    assert(profile.empty()==false);
//    cerr<<"********"<< profile.back().avg_radius<<"  "<<flux_aperture_pixel_radius<<"********"<<endl;
    assert(profile.back().avg_radius >= flux_aperture_pixel_radius);

    size_t radius_index=0;
    double prev_flux_sum=0.0;
    while (profile[radius_index].avg_radius < flux_aperture_pixel_radius
            && radius_index!=profile.size()) {
        prev_flux_sum+=profile[radius_index].mean_value*profile[radius_index].n_samples;
        ++radius_index; // proceed to one past desired level
    }
    assert(radius_index!=profile.size()) ; //sanity check

    double next_flux_sum = prev_flux_sum +
                           profile[radius_index].mean_value*profile[radius_index].n_samples;

    return lucky_math::linearly_interpolate(
               profile[radius_index-1].avg_radius, prev_flux_sum,
               profile[radius_index].avg_radius, next_flux_sum,
               flux_aperture_pixel_radius);
}

//===================================================================================================================
//===================================================================================================================
namespace detail { ///Implementation substructures follow...
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
bool compare_profile_point_radius(const psf_profile_point& first,
                                  const psf_profile_point& second)
{
    return (first.avg_radius < second.avg_radius);
}

void psf_profile_point::save_profile_information_to_file(const std::string& filename,
        const vector<psf_profile_point>& pts,
//        const double peak_val,
        const std::string output_prefix)
{
    ofstream profile_datfile(filename.c_str());
    //------------------------------------------------
    //Header line
    profile_datfile<<output_prefix+"Radius "
                   <<output_prefix+"mean "
                   <<output_prefix+"median "
                   <<output_prefix+"std_dev "
                   <<output_prefix+"min "
                   <<output_prefix+"max ";
//            if (extended_format){
//        profile_datfile
//                <<output_prefix+"normalised_avg_val "
//                <<output_prefix+"proportional_range_min "
//                <<output_prefix+"proportional_range_max  "
//                <<output_prefix+"sigma_range_min "
//                <<output_prefix+"sigma_range_max "
//                <<output_prefix+"5sigma_detection_ratio ";
//                <<output_prefix+"8sigma_detection_mag_diff";
//        profile_datfile<<"\n#Peak value used: "+string_utils::ftoa(peak_val);
//            }
    profile_datfile<<endl;

    //------------------------------------------------
    //Data lines
    for (size_t i=0; i!=pts.size(); ++i) {
        profile_datfile
                <<pts[i].avg_radius<<" "
                <<pts[i].mean_value<<" "
                <<pts[i].median_value<<" "
                <<pts[i].std_dev<<" "
                <<pts[i].min_val<<" "
                <<pts[i].max_val<<" ";
//                if (extended_format){
//            double prop_range_min = (pts[i].min_val - pts[i].mean_value)/ peak_val;
//            double prop_range_max = (pts[i].max_val - pts[i].mean_value)/ peak_val;
//            double sigma_range_min = (pts[i].min_val - pts[i].mean_value)/ pts[i].std_dev;
//            double sigma_range_max = (pts[i].max_val - pts[i].mean_value)/ pts[i].std_dev;
//            profile_datfile
//                    <<pts[i].mean_value/ peak_val<<" "
//                    <<prop_range_min<<" "
//                    <<prop_range_max<<" "
//                    <<sigma_range_min<<" "
//                    <<sigma_range_max<<" ";
//                }
//        double detection_ratio=pts[i].std_dev*8.0/peak_val;
//        profile_datfile<<detection_ratio<<" ";
//        profile_datfile<<-2.5*log10(detection_ratio)<<" ";

        profile_datfile<<endl;
    }
}

//---------------------------------------------------------------------------------------------------------------

bool double_pair_first_member_predicate(const std::pair <double, double>  & first,
                                        const std::pair <double, double>  & second)
{
    return first.first < second.first;
}

//---------------------------------------------------------------------------------------------------------------
//vector<pair <double, double>  > get_radial_data(const PixelArray2d<double>& bmp, const CCD_Position& CCD_lock, double max_radius){
//    PixelPosition frame_pos = bmp.corresponding_frame_Position(CCD_lock);
//    return get_radial_data(bmp, frame_pos, max_radius);
//}
//---------------------------------------------------------------------------------------------------------------

//std::vector<std::pair <double, double>  >
//get_radial_data(    const PixelArray2d<double>& bmp,
//                    const CCD_Position& origin,
//                    const double max_radius,
//                    const PixelArray2d<double>& data_mask){
//    return get_radial_data(bmp, bmp.corresponding_frame_Position(origin), max_radius,data_mask);
//}


//---------------------------------------------------------------------------------------------------------------
vector<pair <double, double>  > get_radial_data(const PixelArray2d<double>& img,
        const PixelPosition& origin, const double max_radius)
{
    PixelArray2d<double> blank;
    return get_radial_data(img, origin, max_radius, blank);
}
//---------------------------------------------------------------------------------------------------------------

std::vector<std::pair <double, double>  > get_radial_data(
    const PixelArray2d<double>& bmp, const PixelPosition& origin, const double max_radius,
    const PixelArray2d<double>& data_mask)
{
    vector<pair <double,double> > distance_value_pairs;

    PixelBoxRegion original_region(origin, origin);
    original_region = PixelBoxRegion::padded_copy(original_region, max_radius);

    PixelBoxRegion overlap_region =
        PixelBoxRegion::overlapping_region(original_region,
                                           PixelBoxRegion::pixel_box_outline(bmp.range())
                                          );

    PixelRange overlap_box=  overlap_region.shrink_to_pixel_boundaries().bounded_pixels();

    if (data_mask.range().n_pix()) {
        //If a valid data mask is supplied:
        assert(data_mask.range()==bmp.range());
        for (PixelIterator it(overlap_box); it!=it.end; ++it) {
            double distance = coord_distance(PixelPosition::centre_of_pixel(it), origin);
            if (distance <= max_radius && data_mask(it)) {
                //non-zero mask value
                distance_value_pairs.push_back(make_pair(distance, bmp(it)));
            }
        }
    } else { //no masking:
        for (PixelIterator it(overlap_box); it!=it.end; ++it) {
            double distance = coord_distance(PixelPosition::centre_of_pixel(it), origin);
            if (distance <= max_radius) {
                distance_value_pairs.push_back(make_pair(distance, bmp(it)));
            }
        }
    }
    return distance_value_pairs;
}

//---------------------------------------------------------------------------------------------------------------
//
vector<psf_profile_point> bin_average_radial_data(vector< pair<double,double> >& raw_data)
{
    assert(!raw_data.empty());
    const double bin_width=1.0;
    vector <double> radii_bin, values_bin;
    sort(raw_data.begin(),
         raw_data.end()); //Data now sorted by radial distance from Position lock
    const double initial_bin_lower_bound =  -0.5*bin_width;
//        (raw_data[0].first < 0.0 ) ? (   ( floor(raw_data[0].first/bin_width)  -0.5)*bin_width     ) : -0.5*bin_width ; //want bins to be centred about origin.

    vector<psf_profile_point> averaged_profile;
    psf_profile_point pt;
    double bin_upper_bound= initial_bin_lower_bound + bin_width;
//    cout <<"Initial bin bounds: " << initial_bin_lower_bound <<","<<bin_upper_bound<<endl;

    for (size_t i=0; i!=raw_data.size(); ++i) {
//        cout <<"Bin upper bound" << bin_upper_bound<<endl;
//        cout <<"Data point " <<raw_data[i].first <<","<<raw_data[i].second<<endl;
        if (raw_data[i].first<=bin_upper_bound) {
            //this point belongs to this bin, add it and move to next data point
//            cout <<"Adding point to bin" <<endl;
            radii_bin.push_back(raw_data[i].first);
            values_bin.push_back(raw_data[i].second);
        } else if (radii_bin.empty() ==false) {
            //this point belongs to the next bin so process this bin and go back to this datapoint
//            cout <<"taking averages " <<endl;
            pt.avg_radius= vector_mean(radii_bin);
            pt.mean_value = vector_mean(values_bin);
            pt.median_value = vector_median(values_bin);
            pt.std_dev = sqrt(vector_variance(values_bin));
            pt.n_samples = radii_bin.size();
            pt.max_val = *max_element(values_bin.begin(), values_bin.end());
            pt.min_val = *min_element(values_bin.begin(), values_bin.end());
            averaged_profile.push_back(pt);
            radii_bin.clear();
            values_bin.clear();
            bin_upper_bound+=bin_width;
            --i; //Go back to the data point we just rejected so it can go in the next bin.
        } else { //this bin has no points in it at all, increment bin and go back to this datapoint
            bin_upper_bound+=bin_width;
            --i;
        }
    }
    assert(averaged_profile.empty()==false);
    return averaged_profile;
}

//---------------------------------------------------------------------------------------------------------------
psf_profile_point interpolate_profile_to_radius(const vector<psf_profile_point>& pts,
        const double interp_radius)
{
    size_t i=0;
    while (i!=pts.size() && pts[i].avg_radius < interp_radius) { ++i; }

    if (i==pts.size()) throw
        runtime_error("Error in \"interpolate_profile_to_radius\": Profile does not extend this far:"
                      +
                      string_utils::ftoa(interp_radius));

    if (interp_radius==pts[i].avg_radius) { return pts[i]; }

    if (i==0) { return pts[0]; }

    psf_profile_point prev_pt(pts[i-1]), next_pt(pts[i]);
    psf_profile_point interp_pt;
    interp_pt.avg_radius=interp_radius;
    interp_pt.mean_value=
        lucky_math::linearly_interpolate(prev_pt.avg_radius, prev_pt.mean_value,
                                         next_pt.avg_radius, next_pt.mean_value,
                                         interp_radius);
    interp_pt.median_value =
        lucky_math::linearly_interpolate(prev_pt.avg_radius, prev_pt.median_value,
                                         next_pt.avg_radius, next_pt.median_value,
                                         interp_radius);
    interp_pt.std_dev=
        lucky_math::linearly_interpolate(prev_pt.avg_radius, prev_pt.std_dev,
                                         next_pt.avg_radius, next_pt.std_dev,
                                         interp_radius);
    interp_pt.n_samples=1;

    return interp_pt;
}

//---------------------------------------------------------------------------------------------------------------
}//end namespace coela::psf_characterisation::detail
//===================================================================================================================

}//end namespace coela::psf_characterisation
}//end namespace coela