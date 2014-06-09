/*
 * File:   psf_generators.h
 * Author: tim
 *
 * Created on 06 January 2010, 21:25
 */

#ifndef COELA_PSF_GENERATORS_H
#define COELA_PSF_GENERATORS_H



#include "coela_core/src/ccd_image.h"
#include "psf_models.h"

//To do - merge all generators into a single class.
//Models should inherit an abc interface to allow psf generation

namespace coela {
namespace psf_models {

//==================================================================================
struct reference_psf {
    CCDImage<double> psf_image;
    PixelIndex central_pixel;
    PixelPosition exact_centre;
    CCDImage<double> mask;

    void apply_mask() {psf_image.pix *= mask.pix;}
//    void set_mask_as_aperture(const PixelPosition& aperture_centre, const double radius_in_frame_pixels);
    vector<double> get_unmasked_ref_psf_pixel_values() const;
};

//==================================================================================
///"oversampling_factor" refers to how many times the construction image is shrunk to produce the output (i.e. degree of numerical integration)
/// -(thereby performing crude integration of the continuous func. to get more accurate pixel values)
//FIXME to use gen_psf_with_params....

reference_psf generate_psf(const psf_model_interface&,
                           const PixelRange& output_frame_shape,
                           const CCD_Position& output_CCD_region_low_corner,
                           const CCD_Position& centre_point,
                           const double desired_pixel_scale_relative_to_CCD_pix,
                           const double construction_oversampling_factor,
                           const bool radius_limited=false,
                           const double max_CCD_radius_required=0.0
                          );

//==================================================================================

reference_psf generate_centred_psf(const axisymmetric_psf_model_interface& psf_func,
                                   const double outer_radius_in_CCD_pix,
                                   const double desired_pixel_scale_relative_to_CCD_pix,
                                   const int oversampling_factor
                                  );
//==================================================================================
//==================================================================================




//        class moffat_psf_generator: public axisymmetric_psf_base, public psf_models::moffat_psf_model{
//        public:
//            moffat_psf_generator(){}
//            moffat_psf_generator(const psf_models::moffat_psf_model& mm):moffat_psf_model(mm){}
//            moffat_psf_generator(double peak_value, double sigma_in_CCD_pixels, double beta_):moffat_psf_model(peak_value,sigma_in_CCD_pixels,beta_){}
//            double operator()(double radius_in_CCD_pix) const;
//        };
//
//        class gauss_moffat_psf_generator: public axisymmetric_psf_base{
//        public:
//            gauss_moffat_psf_generator():combined_peak_val(-1.0),gauss_weight(-1.0){}
//            gauss_moffat_psf_generator(const gauss_psf_generator& g_component, const moffat_psf_generator& moff_component);
//            double operator()(double radius_in_CCD_pix) const;
//            double get_peak_value()const {return combined_peak_val;}
//            void set_peak_value(double combined_peak_value);
//            double get_gauss_component_weight()const {return gauss_weight;}
//            void set_gauss_component_weight(double gauss_peak_weight);
//
//            gauss_psf_generator gauss;
//            moffat_psf_generator moffat;
//        private:
//            double combined_peak_val, gauss_weight;
//            void derive_and_set_component_peaks();
//
//        private:
//            friend class boost::serialization::access;
//            template<class Archive>
//            void serialize(Archive & ar, const unsigned int /* file_version */){
//                ar  & BOOST_SERIALIZATION_NVP(combined_peak_val)
//                    & BOOST_SERIALIZATION_NVP(gauss_weight)
//                    & BOOST_SERIALIZATION_NVP(gauss)
//                    & BOOST_SERIALIZATION_NVP(moffat);
//            }
//        };


}
}//end namespaces

#endif  /* _PSF_GENERATORS_H */

