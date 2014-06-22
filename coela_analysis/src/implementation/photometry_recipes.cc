#include "../photometry_recipes.h"


namespace coela {
namespace photometry {

double aperture_flux(const PixelArray2d<double>& img,
                     const PixelPosition& centre, const double pixel_radius,
                     size_t* n_pixels_in_aperture_ret_ptr)
{

    regions::CircularAperture<coordinate_types::pixels>
    photometry_ap(PixelPosition::centre_of_pixel(
                      PixelPosition::pixel_containing_point(centre)),
                  pixel_radius);

    PixelRange photometry_box =
        photometry_ap.get_covering_region().expand_to_pixel_boundaries().bounded_pixels();

    size_t n_pix_summed=0;
    double flux_sum=0.0;
    for (PixelIterator it(photometry_box); it!=it.end; ++it) {
        if (photometry_ap.contains_point(PixelPosition::centre_of_pixel(it))) {
            flux_sum+= img(it);
            n_pix_summed++;
        }
    }
    if (n_pixels_in_aperture_ret_ptr) { *n_pixels_in_aperture_ret_ptr=n_pix_summed; }
    return flux_sum;
}


}
}
