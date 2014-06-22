/*
 * File:   fitting_noise_models.h
 * Author: ts337
 *
 * Created on 13 May 2011, 18:31
 */

#ifndef COELA_FITTING_NOISE_MODELS_H
#define COELA_FITTING_NOISE_MODELS_H
#include "coela_core/src/pixel_array2d.h"

namespace coela {
namespace psf_fitting {
//==================================================================================
//Interface class allowing for fitting with different noise models
class variance_estimation_base {
public:
    //return the estimated variance for this pixel:
    virtual double operator()(const PixelIndex & input_data_pixel_Position) const=0;
};
//==================================================================================
class least_squares_weighting: public variance_estimation_base {
public:
    double operator()(const PixelIndex&) const {return 1.0;} //Weights all pixels equally.
};

/////Essentially weights pixels as photon variance + gaussian BG variance. (see details below for stochastic weighting...)
//class photon_counted_weighting: public variance_estimation_base{
//public:
//    photon_counted_weighting(const float_bitmap& input_image, int n_frames_used, double background_variance);
//    photon_counted_weighting(const float_bitmap& input_image, const float_bitmap& weight_map, double background_variance);
//    float operator()(size_t input_x_index, size_t input_y_index) const {return pixel_noise_estimates(input_x_index,input_y_index);}
//private:
//    float_bitmap pixel_noise_estimates;
//};
//
/////Assumes the input image has been normalised to avg photon counts, and should be fed the background gaussian variance in the averaged image
/////Essentially weights pixels as 2 * photon variance + gaussian BG variance.
////Standard constructor assumes uniform drizzle weighting used to create an avg image (i.e. divided by a factor N, so var divided by N^2.).
////For non-uniform drizzled images, construct with the weighting bitmap
////Formula for noise variance is VAR = N*(avg photon rate)*2/N^2 + BG Var. = 2*avg_photon_rate/N + BG_Var.
//class stochastic_EMCCD_weighting: public variance_estimation_base{
//public:
//    stochastic_EMCCD_weighting(const float_bitmap& input_image, int n_frames_used, double background_variance);
//    stochastic_EMCCD_weighting(const float_bitmap& input_image, const float_bitmap& weight_map, double background_variance);
//    float operator()(size_t input_x_index, size_t input_y_index) const {return pixel_noise_estimates(input_x_index,input_y_index);}
//private:
//    float_bitmap pixel_noise_estimates; //Initialised in constructor
//
//};

} //end namespace coela::PsfFitting
}//end namespace coela

#endif  /* FITTING_NOISE_MODELS_H */

