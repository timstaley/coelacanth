/*
 * File:   photometry_recipes.h
 * Author: ts337
 *
 * Created on 25 May 2011, 13:21
 */

#ifndef COELA_PHOTOMETRY_RECIPES_H
#define COELA_PHOTOMETRY_RECIPES_H

#include "regions.h"
#include "coela_core/src/pixel_array2d.h"

namespace coela {
namespace photometry {

double aperture_flux(const PixelArray2d<double>& img,
                     const PixelPosition& centre, const double pixel_radius,
                     size_t* n_pixels_in_aperture_return_value=NULL);


}
}


#endif  /* PHOTOMETRY_RECIPES_H */

