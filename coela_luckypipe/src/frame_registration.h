/*
 * File:   frame_registration.h
 * Author: ts337
 *
 * Created on 18 February 2011, 15:15
 */

#ifndef COELA_FRAME_REGISTRATION_H
#define COELA_FRAME_REGISTRATION_H

#include "coela_core/src/ccd_image.h"
#include "coela_analysis/src/psf_generation.h"

namespace coela {

namespace frame_registration {

//=========================================================================================
///"Guide star lock" - A struct to hold a co-ordinate Position, and a quality rating for a guide star
template<class Position_class>
struct GuideStarLock {
public:
    //-------------------------------------------------
    //data
    Position_class Position;
    double signal;
    //-------------------------------------------------
    GuideStarLock() {}
    GuideStarLock(Position_class p): Position(p),signal(0.0) {}
    GuideStarLock(Position_class p, double s): Position(p), signal(s) {}
};

template<class Position_type>
std::ostream& operator<<(std::ostream& os, const GuideStarLock<Position_type>& gs) {os << gs.signal <<" " << gs.Position; return os;}

template<class Position_type>
std::istream& operator>>(std::istream& is, GuideStarLock<Position_type>& gl) {is>>gl.signal; is>> gl.Position; return is;}

//==============================================================================================================

//TO DO: ad hoc, needs thought, don't have time.


psf_models::ReferencePsf convert_to_log_values(const psf_models::ReferencePsf&);
//==============================================================================================================


///This is a shorthand for "resample, convolve with kernel (perhaps only at threshold), find max pixel and work back to co-ordinates of original frame"
template<typename input_datatype>
GuideStarLock<CcdPosition> find_best_psf_match(
    const CcdImage<input_datatype>& input_image,
    const CcdBoxRegion& gs_region,
    const psf_models::ReferencePsf& ref_psf,
    const double input_resample_factor,
    const double convolution_threshold_factor,
    const bool fit_parabola=false);

//
//GuideStarLock calculate_Position(const convolution_kernel& kernel, const image<double>& resampled_gs_rgn, const image<double>& convolved_gs_rgn, bool fit_parabola=true);
//
//GuideStarLock convolve_and_calculate_Position(const convolution_kernel& kernel, const image<double>& resampled_gs_rgn, const double convolution_threshold_factor=0.0, bool fit_parabola=true);


//==============================================================================================================
} //end namespace coela::frame_registration
}


#endif  /* FRAME_REGISTRATION_H */

