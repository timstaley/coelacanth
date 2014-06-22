/*
 * File:   psf_models.cc
 * Author: ts337
 *
 * Created on 13 May 2010, 13:30
 */

#include "../psf_models.h"
#include "coela_utility/src/misc_math.h"
#include "coela_utility/src/string_utils.h"
#include <stdexcept>

namespace coela {
namespace psf_models {

//================================================================================
double calculate_total_flux(const axisymmetric_psf_model_interface& psf_gen,
                            const double step_size,
                            const size_t CCD_pixel_limit,
                            const double precision)
{
    double delta=0;
    double sum=0.0;
    double psf_val=0.0;
    //first calculate for origin patch (0 - r_step/2)
    sum += psf_gen(0.0) * M_PI *(step_size/2)*(step_size/2);

    //now for annuli of increasing radius ( n*r_step-1/2 to n*r_step+1/2
    size_t n_steps=1;
    do {
        double radius = n_steps*step_size;
        psf_val = psf_gen(radius) ;
        delta =2*M_PI*radius*psf_val* step_size;
        sum+=delta;
        n_steps++;
    } while (psf_val > precision && n_steps*step_size<CCD_pixel_limit);

//            cout <<"Truncated at ccd_radius" <<n_steps*step_size<<" at psf val " << psf_val<<endl;

    if (n_steps*step_size>= CCD_pixel_limit) throw
        std::runtime_error("Flux calculation did not converge at "
                           +string_utils::itoa(CCD_pixel_limit)+" CCD pixels radius,  PSF value is still"
                           + string_utils::ftoa(psf_val));
    return sum;
}

double calculate_flux_enclosed_at_radius(
    const axisymmetric_psf_model_interface& psf_gen,
    const double max_radius,
    const double step_size
)
{
    double sum=0.0;
    sum += psf_gen(0.0) * M_PI *(step_size/2)*(step_size/2);
    size_t n_steps=1;
    double radius;
    do {
        radius = n_steps*step_size;
        double psf_val = psf_gen(radius) ;
        sum+= 2*M_PI*radius*psf_val* step_size;
        n_steps++;
    } while (radius < max_radius);
    return sum;
}

double calculate_strehl(const axisymmetric_psf_model_interface& psf_gen,
                        const double airy_total_flux)
{
    double psf_flux_for_unit_peak = calculate_total_flux(psf_gen) / psf_gen(0.0);
    return airy_total_flux / psf_flux_for_unit_peak ;
}

double fwhm_in_CCD_pix(const axisymmetric_psf_model_interface& psf_gen,
                       const double CCD_pixel_width_precision)
{
    double half_max = psf_gen(0.0)/2.0;
    double half_width=0;
    while (psf_gen(half_width)> half_max) { half_width+=1.0; }
    double low_guess= half_width-1.0;
    double high_guess= half_width;
    //current_guess=half_width
    while ((high_guess-low_guess) > CCD_pixel_width_precision) {
        double mid_point = (high_guess + low_guess) / 2.0;

        if (psf_gen(mid_point) > half_max) {
            //need to move to higher radius to achieve lower psf value, raise low guess
            low_guess=mid_point;
        } else {
            high_guess=mid_point;
        }
    }
    half_width=low_guess;
    return half_width*2.0;
}

double find_first_minima_in_CCD_pix(const axisymmetric_psf_model_interface& psf_gen,
                                    const double CCD_pixel_stepsize)
{
    double prev_val = psf_gen(0.0);
    double r=CCD_pixel_stepsize;
    double curr_val=psf_gen(r);
    while (curr_val < prev_val && r<100) {
        prev_val=curr_val;
        r+=CCD_pixel_stepsize;
        curr_val=psf_gen(r);
    }
    if (r>=100) { throw std::runtime_error("No psf minima found"); }
    return r-CCD_pixel_stepsize;
}

//================================================================================
airy_psf_model::airy_psf_model()
    :peak_val(-1.0), rads_per_CCD_pix(-1.0), wavelength(-1.0),
     aperture_diameter(-1.0), obscuration_diameter(-1.0)
{}

airy_psf_model::airy_psf_model(double peak_value, double rads_per_CCD_pixel,
                               double wavelength_,
                               double aperture_diameter_, double central_obscuration_diameter_)
    :peak_val(peak_value), rads_per_CCD_pix(rads_per_CCD_pixel), wavelength(wavelength_),
     aperture_diameter(aperture_diameter_), obscuration_diameter(central_obscuration_diameter_)
{}

double airy_psf_model::operator()(double radius_in_CCD_pix) const
{
    return peak_val*
           misc_math::airy_function(radius_in_CCD_pix*rads_per_CCD_pix,
                                     wavelength, aperture_diameter, obscuration_diameter);
}


//================================================================================

gaussian_psf_model::gaussian_psf_model()
    :peak_val(-1.0),sigma_in_CCD_pix(-1.0)
{}

gaussian_psf_model::gaussian_psf_model(double peak_value,
                                       double sigma_measured_in_CCD_pixels)
    :peak_val(peak_value), sigma_in_CCD_pix(sigma_measured_in_CCD_pixels)
{}

double gaussian_psf_model::operator()(double radius_in_CCD_pix) const
{
    return misc_math::gaussian_1d_function(radius_in_CCD_pix, peak_val, sigma_in_CCD_pix);
}

//================================================================================
empirical_psf_model::empirical_psf_model(
    const std::vector<psf_characterisation::detail::psf_profile_point>& pts)
    :datapoints(pts)
{}

double empirical_psf_model::operator()(double radius_in_CCD_pix) const
{
    return get_interpolated_point_at_radius(radius_in_CCD_pix).median_value;
}

psf_characterisation::detail::psf_profile_point
empirical_psf_model::get_interpolated_point_at_radius(double radius_in_CCD_pix) const
{
    return psf_characterisation::detail::interpolate_profile_to_radius(
               datapoints, radius_in_CCD_pix);
}


//================================================================================

hybrid_radial_switchover_psf_model::hybrid_radial_switchover_psf_model(
    const axisymmetric_psf_model_interface* inner_model_ptr,
    const axisymmetric_psf_model_interface* outer_model_ptr,
    const double switchover_radius,
    const double normalization_factor)
    :inner_model_ptr_(inner_model_ptr), outer_model_ptr_(outer_model_ptr),
     switchover_radius_(switchover_radius),
     normalization_factor_(normalization_factor)
{}

double hybrid_radial_switchover_psf_model::operator()(double radius_in_CCD_pix) const
{
    if (radius_in_CCD_pix<switchover_radius_) {
        return normalization_factor_*(*inner_model_ptr_)(radius_in_CCD_pix) ;
    } else {
        return normalization_factor_*(*outer_model_ptr_)(radius_in_CCD_pix);
    }
}

//================================================================================
} //end namespace coela::psf_models
}//end namespace coela



