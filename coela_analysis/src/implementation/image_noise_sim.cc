/*
 * File:   image_noise_sim.cc
 * Author: ts337
 *
 * Created on 22 March 2011, 18:49
 */


#include "../image_noise_sim.h"
#include <boost/random/poisson_distribution.hpp>
#include<iostream>

namespace coela {
namespace image_noise_simulation {


uniform_random_pixel_variate::uniform_random_pixel_variate(
    const PixelRange& image_outline_,
    unuran::StreamWrapper& rns)
    : urv(0.0, image_outline_.n_pix(), rns), range(image_outline_)
{}

PixelIndex uniform_random_pixel_variate::operator()()
{
    int data_vec_index = int(urv()) ;
    return range.get_PixelIndex_for_data_vector_element(data_vec_index);
}


intensity_map_random_pixel_variate::intensity_map_random_pixel_variate(
    const PixelArray2d<double>& intensity_map,
    unuran::StreamWrapper& rns,
    const size_t /*expected_number_of_uses*/)
    :gen(NULL), range(intensity_map.range())
{
    UNUR_DISTR *distr;    /* distribution object                   */
    UNUR_PAR   *par;      /* parameter object                      */

    distr = unur_distr_discr_new();

    const vector<double>& prob_vec(intensity_map.get_raw_data());

    unur_distr_discr_set_pv(distr, &prob_vec[0] , prob_vec.size());
//    par = unur_dgt_new(distr);
//    unur_dgt_set_guidefactor(par, std::min(3.0/sqrt(prob_vec.size()), 1.0) );

    par = unur_dau_new(distr);
    unur_set_urng(par, rns.get_stream_pointer());
    gen = unur_init(par);
    unur_distr_free(distr);

    if (gen==NULL) { throw std::runtime_error("intensity_map_random_pixel_variate::intensity_map_random_pixel_variate - could not construct generator."); }

}

PixelIndex intensity_map_random_pixel_variate::operator()()
{
    int data_vec_index = unur_sample_discr(gen);
    return range.get_PixelIndex_for_data_vector_element(data_vec_index);
}


PixelArray2d<int> generate_photon_arrival_map(
    const PixelArray2d<double>& source_intensity_map,
    const double mean_total_number_of_source_photons,
    const double mean_bg_photons_per_pixel,
    unuran::StreamWrapper& rns)
{
    
    unuran::PoissonRandomVariate source_total_prv(
        mean_total_number_of_source_photons, rns);

    int n_source_photons = source_total_prv();
//    return photon_arrival_map_fixed_n_low_flux_optimized(intensity_map, n_photons, rns);

    PixelArray2d<int> source_photon_arrival_map(
        source_intensity_map.range().x_dim(), source_intensity_map.range().y_dim(), 0);

    intensity_map_random_pixel_variate pixel_picker(source_intensity_map, rns);

    for (int phot_index=0; phot_index!=n_source_photons; ++phot_index) {
        source_photon_arrival_map(pixel_picker())+=1;
    }

    if (mean_bg_photons_per_pixel!=0.0) {
        double mean_bg_sum_flux =
            source_photon_arrival_map.range().n_pix() * mean_bg_photons_per_pixel;
        
        unuran::PoissonRandomVariate bg_total_prv(mean_bg_sum_flux, rns);
        int n_bg_photons = bg_total_prv();
        uniform_random_pixel_variate bg_pixel_picker(source_photon_arrival_map.range(),
                rns);
        for (int bg_phot_gen=0; bg_phot_gen!=n_bg_photons; ++bg_phot_gen) {
            source_photon_arrival_map(bg_pixel_picker()) +=1;
        }
    }
    return source_photon_arrival_map;

}





PixelArray2d<int> EMCCD_simulated_image(
    const PixelArray2d<int> & photon_arrival_map,
    EmccdModel& model_variate,
    unuran::StreamWrapper& rns)
{
    PixelArray2d<int> sim_image(photon_arrival_map.range().x_dim(),
                                photon_arrival_map.range().y_dim(), 0);

    //Simulate photon events
    for (PixelIterator i(photon_arrival_map.range()); i!=i.end; ++i) {
        if (photon_arrival_map(i)) {
            sim_image(i) +=
                model_variate.stochastically_multiply_photons(photon_arrival_map(i));
        }
    }

    //simulate CIC events
    if (model_variate.params.serial_CIC_rate!=0) {
        int n_CIC_events;

        double expected_num_CIC_events = model_variate.params.serial_CIC_rate *
                                         sim_image.range().n_pix() ;

        unuran::PoissonRandomVariate n_CIC_events_variate(
            expected_num_CIC_events,
            rns);
        n_CIC_events = n_CIC_events_variate();
        uniform_random_pixel_variate pixel_picker(sim_image.range(), rns);

        for (int i=0; i!=n_CIC_events; ++i) {
//            int pixel_number = int( pixel_picker() ) + 1;
//            PixelIndex pix = sim_image.range().get_nth_PixelIndex(pixel_number);
            PixelIndex pix = pixel_picker();
            sim_image(pix) += model_variate.CICIR_variate();
        }
    }
    //Simulate readout
    for (PixelIterator i(sim_image.range()); i!=i.end; ++i) {
        sim_image(i)+=model_variate.readout_noise_variate();
    }

    return sim_image;
}

}//end namespace coela::image_noise_simulation
}//end namespace coela
