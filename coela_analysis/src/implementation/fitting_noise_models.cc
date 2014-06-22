/*
 * File:   fitting_noise_models.cc
 * Author: ts337
 *
 * Created on 13 May 2011, 18:31
 */

#include "../fitting_noise_models.h"

namespace coela {
namespace psf_fitting {


//photon_counted_weighting::photon_counted_weighting(const float_bitmap& input_image, int n_frames_used, double background_variance){
//    pixel_noise_estimates=input_image;
//    for (size_t i=0; i!=input_image.n_elements();++i){
//        pixel_noise_estimates.element(i)=  max(input_image.element(i),0.f)/n_frames_used + background_variance;
//    }
//}
//
//photon_counted_weighting::photon_counted_weighting(const float_bitmap& input_image, const float_bitmap& weight_map, double background_variance){
//    pixel_noise_estimates=input_image;
//    for (size_t i=0; i!=input_image.n_elements();++i){
//        pixel_noise_estimates.element(i)=  max(input_image.element(i),0.f)/max(weight_map.element(i),1.f) + background_variance;
//    }
//}
//
//stochastic_EMCCD_weighting::stochastic_EMCCD_weighting(const float_bitmap& input_image,
//        int n_frames_used, double background_variance){
//    pixel_noise_estimates=input_image;
//    for (size_t i=0; i!=input_image.n_elements();++i){
//        pixel_noise_estimates.element(i)=  max(input_image.element(i),0.f)*2.0/n_frames_used + background_variance;
//    }
//}
//
//stochastic_EMCCD_weighting::stochastic_EMCCD_weighting(const float_bitmap& input_image, const float_bitmap& weight_map, double background_variance){
//    pixel_noise_estimates=input_image;
//    for (size_t i=0; i!=input_image.n_elements();++i){
//        pixel_noise_estimates.element(i)=  max(input_image.element(i),0.f)*2.0/max(weight_map.element(i),1.f) + background_variance;
//        //nb if the weight map is 0 for this pixel, then the pixel should be masked during the fitting procedure
//        //so the noise weighting in that case is moot.
//    }
//}


} //end namespace coela::PsfFitting
}//end namespace coela