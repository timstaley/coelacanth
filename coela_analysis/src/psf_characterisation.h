#ifndef COELA_PSF_CHARACTERISATION_H
#define COELA_PSF_CHARACTERISATION_H

#include "coela_core/src/pixel_array2d.h"
#include "coela_core/src/cartesian_coords.h"


namespace coela {
namespace psf_characterisation {
//double estimate_fwhm_in_image_pix(const PixelArray2d<double>& bg_subbed_bmp, const CCD_Position& origin, const double est_peak_value, const double max_radius); ///<uses linear interpolation to estimate FWHM, where 'maximum' is value of origin pixel.
double estimate_fwhm_in_image_pix(const PixelArray2d<double>& bg_subbed_bmp,
                                  const PixelPosition& origin,
                                  const double est_peak_value,
                                  const double max_radius);

double estimate_radius_of_avg_value(const PixelArray2d<double>& bmp,
                                    const PixelPosition& origin,
                                    const double max_radius,
                                    const double value_to_determine_radius_of);

double estimate_radius_to_enclose_flux(const PixelArray2d<double>& bmp,
                                       const PixelPosition& origin,
                                       const double flux_sum_to_enclose,
                                       const double max_radius
                                      );

double estimate_FWHEF_in_image_pix(const PixelArray2d<double>& bg_subbed_bmp,
                                   const PixelPosition& origin,
                                   const double total_flux,
                                   const double max_radius);

double estimate_encircled_flux_at_pixel_radius(const PixelArray2d<double>& bmp,
        const PixelPosition& origin,
        const double pixel_radius
                                              );


//===================================================================================================================
//===================================================================================================================
namespace detail { ///Implementation substructures follow...
//---------------------------------------------------------------------------------------------------------------

struct psf_profile_point {
    double avg_radius;
    double mean_value, median_value, std_dev, max_val, min_val;
    size_t n_samples;

    static void save_profile_information_to_file(const std::string& filename,
            const std::vector<psf_profile_point>& pts,
//            const double peak_val=1.0,
            const std::string output_prefix="");
};

bool compare_profile_point_radius(const psf_profile_point& first,
                                  const psf_profile_point& second);
//---------------------------------------------------------------------------------------------------------------

//std::vector<std::pair <double, double>  > get_radial_data(const PixelArray2d<double>& bmp, const CCD_Position& origin, double max_radius);



//---------------------------------------------------------------------------------------------------------------

//std::vector<std::pair <double, double>  >
//get_radial_data(    const PixelArray2d<double>& bmp,
//                    const CCD_Position& origin,
//                    const double max_radius,
//                    const PixelArray2d<double>& data_mask);
//---------------------------------------------------------------------------------------------------------------

bool double_pair_first_member_predicate(const std::pair <double, double>  & first,
                                        const std::pair <double, double>  & second);

///Get a raw vector of pairs of <radius, pixel value > , where radius is from a designated centre (origin)
std::vector<std::pair <double, double>  > get_radial_data(const PixelArray2d<double>& img,
        const PixelPosition& origin,
        const double max_radius);

///Optionally; supply a data mask so that certain pixels are ignored.
std::vector<std::pair <double, double>  > get_radial_data(const PixelArray2d<double>& bmp,
        const PixelPosition& origin,
        const double max_radius,
        const PixelArray2d<double>& data_mask);

//---------------------------------------------------------------------------------------------------------------
std:: vector<psf_profile_point> bin_average_radial_data(
    vector< std::pair<double,double> >& raw_data);
//---------------------------------------------------------------------------------------------------------------
psf_profile_point interpolate_profile_to_radius(const vector<psf_profile_point>& pts,
        const double interp_radius);


//===================================================================================================================
}//end namespace coela::psf_characterisation::detail
//===================================================================================================================
//===================================================================================================================



}//end namespace coela::psf_characterisation
}//end namespace coela

#endif // PSF_CHARACTERISATION_H