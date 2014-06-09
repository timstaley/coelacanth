/*
 * File:   image_noise_sim.h
 * Author: ts337
 *
 * Created on 22 March 2011, 18:18
 */

#ifndef COELA_IMAGE_NOISE_SIM_H
#define COELA_IMAGE_NOISE_SIM_H

#include "coela_random/src/random.h"
#include "coela_random/src/emccd_model.h"
#include "coela_core/src/ccd_image.h"
namespace coela {
namespace image_noise_simulation {

/// Picks out a pixel at random (e.g. for use simulating uniform illumination)
class uniform_random_pixel_variate {
public:
    uniform_random_pixel_variate(const PixelRange& image_outline_,
                                 unuran::StreamWrapper & rns);
    PixelIndex operator()();


private:
    unuran::UniformRandomVariate urv;
    const PixelRange range;

};


class intensity_map_random_pixel_variate {
public:
    intensity_map_random_pixel_variate(
        const PixelArray2d<double>& intensity_map,
        unuran::StreamWrapper& rns,
        const size_t expected_number_of_uses=0); //dummy variable, may be implemented later

    PixelIndex operator()();

private:
    UNUR_GEN *gen;
    const PixelRange range;
};



PixelArray2d<int> generate_photon_arrival_map(
    const PixelArray2d<double>& source_intensity_map,
    const double mean_total_number_of_source_photons,
    const double mean_bg_photons_per_pixel,
    unuran::StreamWrapper& rns);


PixelArray2d<int> EMCCD_simulated_image(const PixelArray2d<int> & photon_arrival_map,
                                        EmccdModel& model_variate ,
                                        unuran::StreamWrapper& rns);

}//end namespace coela::image_noise_simulation
}//end namespace coela

#endif  /* IMAGE_NOISE_SIM_H */

