#include "../misc_math.h"
#include "../string_utils.h"
#include <boost/math/special_functions/bessel.hpp>
#include <algorithm>
#include <stdexcept>



using namespace std;
namespace coela {
namespace misc_math {



double linearly_interpolate(double x1, double y1, double x2, double y2, double x_between)
{
//    cerr<<"Inputs: "<<endl;
//    cerr<<"x1, y1: "<< x1 <<","<<y1<<endl;
//    cerr<<"x2, y2: "<< x2 <<","<<y2<<endl;
//    cerr<<"xinterp: "<< x_between<<endl;

    assert(std::min(x1,x2) - err_margin <= x_between
           && std::max(x1,x2) + err_margin >= x_between);

    //Flat line case,  pick arbitrary mid point.
    if (x1==x2) { return (y1+y2)/2.0; }

    double grad = (y2 - y1) / (x2 - x1);
    return y1 + grad*(x_between-x1);
}

double dimensionless_airy_function(double x, double ratio)
{
    //ratio = ratio of DIAMETERS ( or equiv, ratio of radii.) (not area)
    double j = 2.0 * (boost::math::cyl_bessel_j(1, x) -  ratio*boost::math::cyl_bessel_j(1,
                      ratio*x)) / (x*(1.0 - ratio*ratio)) ;     // 2[ J(x)/x - f*f*J(fx)/fx ] / (1-f*f)
    return j*j;
}

double airy_function(double angle, double wavelength, double aperture_diameter,
                     double obscuration_diameter)
{
    double f = obscuration_diameter / aperture_diameter;  //fractional obscuration
    double x = M_PI * aperture_diameter * (angle) /
               wavelength; //rescaled distance from centre
//            double x = M_PI * aperture_diameter * sin(angle) / wavelength; //better for large angles? But irrelevant here
    //avoid divide by 0 - return the value "in the limit" as x-->0:
    if (x==0) { x=1e-9; }
    //simple case: I = ( 2 J1[x] / x )^2
//             double j = 2.0 * ( boost::math::cyl_bessel_j(1, x) -  f*boost::math::cyl_bessel_j(1, f*x) ) / (x*(1.0 - f*f )  );  // 2( J(x) - fJ(fx) ) / x(1-f*f)
//            return j*j;
    return dimensionless_airy_function(x, f);
}

double airy_fwhm_in_rads(double wavelength, double aperture_diameter,
                         double obscuration_diameter)
{
    //Monotonic between peak and first minima- use bisection and linear interpolation NB first minima not at 1.22 if obscuration !=0, but will do for a first approximation
    double lower_limit(0.0), upper_limit(1.22*wavelength/aperture_diameter), eps(1e-5),
           delta(1.0);
    double hwhm_angle(upper_limit/2.0);
    while (fabs(delta)>eps) {
        delta = airy_function(hwhm_angle, wavelength, aperture_diameter,
                              obscuration_diameter) -0.5;
        if (delta >0) {
            //too close to origin
            lower_limit = hwhm_angle;
        } else { upper_limit = hwhm_angle; }
        hwhm_angle = (upper_limit + lower_limit) /2.0;
    }
    return hwhm_angle*2.0;
}

double gaussian_1d_function(double radius, double peak_value, double sigma)
{
    if (sigma!=0) { return peak_value*exp(-1.0*radius*radius/(2.0*sigma*sigma)); }
    else if (radius<1e-7) { return peak_value; }
    return 0;
}

double moffat_function(double radius, double sigma, double beta)
{
    assert(sigma!=0);
    double scale_length = radius/sigma;
    return pow((1+scale_length*scale_length), -1.0*beta);
}

double moffat_fwhm(double sigma, double beta)
{
    return 2.0*sigma* sqrt(pow(2.0, 1.0/beta) - 1);
}

double moffat_sigma(double fwhm, double beta)
{
    double radius_of_half_maximum=fwhm/2.0;
    return radius_of_half_maximum / sqrt((pow(2.0, 1.0/beta) - 1.0));
}

unsigned long integer_factorial(int n)
{
    assert(n<20);
    if (n==0) { return 1ul; }
    return n*integer_factorial(n-1);
}


///Sidesteps the integer limits, but as such is not exact. Used for calculating PDFs for an EMCCD (which are already approximations anyway).
double approx_factorial_double(double number)
{
    if (number <= 1.0) { return 1.0; }
    return number * approx_factorial_double(number - 1.0);
}

double integer_power(double x, int n)
{
    if (n==0) { return 1.0; }
    if (n<0) { return 1.0/integer_power(x,-1*n); }

    double result=x;
    while (n!=1) {
        result*=x;
        n--;
    }
    return result;
}

pair<double,double> fit_1d_parabola_to_find_local_maxima(double x0,
        double y_xm1, double y_x0, double y_xp1)
{
    //nb should have been fed a local maximum, so:
//            using string_utils::ftoa;
//            if (!(y_x0 > y_xm1 && y_x0 > y_xp1))
//                throw std::runtime_error("Fit parabola: Not a local max: "+ftoa(y_x0) +" less than "+ftoa(y_xm1)+"or"+ftoa(y_xp1));
//                ;
    assert(y_x0 >= y_xm1 && y_x0 >= y_xp1);
    pair<double,double> xmax_ymax;
//See Lisa Poyneer 2003, "scene-based SHWFS" for a nicely formatted version... but it's just very simple math.
    double b = -0.5*(y_xm1 - y_xp1);
    double two_a =(y_xm1 + y_xp1 -2*y_x0);
    //c = y_x0

//            xmax_ymax.first= x0 + ( 0.5*(y_xm1 - y_xp1)  )/(y_xm1 + y_xp1 -2*y_x0);
    xmax_ymax.first= x0 - b/two_a;
    xmax_ymax.second=y_x0 - b*b/(2*two_a);

    return xmax_ymax;
}

double bilinear_interpolate(double x, double y,
                            double x0, double y0,
                            double x1, double y1,
                            double f_x0_y0, double f_x1_y0, double f_x0_y1, double f_x1_y1)
{
    double x_frac(x-x0 / (x1-x0)), y_frac(y-y0 / (y1-y0));

    return f_x0_y0 *(1 - x_frac)*(1 - y_frac)
           + f_x1_y0 *(x_frac)*(1 - y_frac)
           + f_x0_y1 *(1-x_frac)*(y_frac)
           + f_x1_y1 *(x_frac)*(y_frac)    ;
}

}
}
