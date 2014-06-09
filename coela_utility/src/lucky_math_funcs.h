/*
 * File:   lucky_math_funcs.h
 * Author: ts337
 *
 * Created on 11 April 2009, 19:49
 */

#ifndef COELA_MATH_FUNCS_H
#define COELA_MATH_FUNCS_H
#include <cmath>
#include <utility>

namespace coela {
namespace lucky_math {

const double err_margin=1e-6;

const float gaussian_fwhm_to_sigma_ratio = 2.3548;
const double rads_to_mas  = 1000*180*3600 /M_PI;
const double mas_to_rads =  M_PI  / (180.0*3600.0 *1000.0) ;

const double radians_per_arcsecond = M_PI / (180.0*3600.0);
const double radians_per_milliarcsecond=  M_PI  / (180.0*3600.0 *1000.0);


double linearly_interpolate(double input1, double dependent1,
                            double input2, double dependent2, double input_Position_to_interpolate);
double airy_function(double angle_in_radians, double wavelength,
                     double aperture_diameter, double aperture_obscuration_diameter=0.0);

double airy_fwhm_in_rads(double wavelength, double aperture_diameter,
                         double aperture_obscuration_diameter=0.0);

double gaussian_1d_function(double radius, double peak_value,
                            double sigma); //returns value at distance 'radius' from the peak
double moffat_function(double radius, double sigma, double beta); //Value at origin is 1.0
double moffat_fwhm(double sigma, double beta);
double moffat_sigma(double fwhm, double beta);



double dimensionless_airy_function(double x, double obscuration_ration);

//NB throws if number is greater than ...(TEST ME)... because it will hit the integer limit.
unsigned long integer_factorial(int n);

///Sidesteps the integer limits, but as such is not exact.
// Used for calculating PDFs for an EMCCD (which are already approximations anyway).
double approx_factorial_double(double n);

double integer_power(double x, int n);

///Returns a pair (x,y) representing the x value and y value at the stationary point ( dy/dx =0 )
/// NB here x is the independent variable, and y = f(x).  i.e. this applies to a 1d curve
std::pair<double,double> fit_1d_parabola_to_find_local_maxima(double x0,
        double y_at_x_minus_1, double y_at_x_0, double y_at_x_plus_1);

double bilinear_interpolate(double x, double y,
                            double x0, double y0,
                            double x1, double y1,
                            double f_x0_y0, double f_x1_y0,
                            double f_x0_y1, double f_x1_y1);

///Was used when testing various optimization techniques:
inline int auxCeil(float x) { return (float)(long)(x+0.999999f) ;} //This is the closest we can get to 1.0 without it rounding up by 1.0.


}//end namespace coela::math
}//end namespace coela


#endif  /* _MATH_FUNCS_H */

