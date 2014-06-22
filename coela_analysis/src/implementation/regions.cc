/*
 * File:   regions.cc
 * Author: ts337
 *
 * Created on 16 May 2011, 17:24
 */

#include "../regions.h"


namespace coela {
namespace regions {

//============================================================================
template<class coord_type>
CircularAperture<coord_type>::CircularAperture(
    const Position<coord_type>& centre_init, const double radius_in_these_coords)
    :centre(centre_init), radius(radius_in_these_coords)
{}


template<class coord_type>
bool CircularAperture<coord_type>::contains_point(const Position<coord_type>& testpoint)
const
{
    if ((testpoint - centre).length() < radius) {
        return true;
    }
    return false;
}
template<class coord_type>
RectangularRegion<coord_type> CircularAperture<coord_type>::get_covering_region() const
{
    Position<coord_type> low_corner(centre), high_corner(centre);
    TwoVector<coord_type> square_half_diagonal(radius,radius);
    low_corner-=square_half_diagonal;
    high_corner+=square_half_diagonal;

    return RectangularRegion<coord_type>(low_corner,high_corner);
}

template<class coord_type>
double CircularAperture<coord_type>::analytic_area() const
{
    return M_PI*radius*radius;
}

template<>
CcdImage<double> CircularAperture<coordinate_types::pixels>::
generate_mask_for_covered_portion_of_image(
    const CcdImage<double>& input_image) const
{

    assert(input_image.CCD_grid.is_initialized());
    PixelBoxRegion img_region = PixelBoxRegion::pixel_box_outline(input_image.pix.range());

    PixelBoxRegion overlap_region =
        PixelBoxRegion::overlapping_region(this->get_covering_region(), img_region);

    overlap_region.expand_to_pixel_boundaries();

    CcdImage<double> mask;
    mask.pix = PixelArray2d<double>(overlap_region.x_dim(), overlap_region.y_dim(), 0.0);

    mask.initialize_CCD_grid_to_specific_region(
        input_image.CCD_grid.corresponding_grid_region(overlap_region));

    PixelShift mask_offset = overlap_region.low - PixelPosition::origin;

    for (PixelIterator it(mask.pix.range()); it!=it.end; ++it) {
        PixelPosition input_pixel_centre = PixelPosition::centre_of_pixel(it) + mask_offset;
        if (this->contains_point(input_pixel_centre)) {
            mask.pix(it)=1.0;
        }
    }
    return mask;

}

template<>
CcdImage<double> CircularAperture<coordinate_types::ccd>::
generate_mask_for_covered_portion_of_image(
    const CcdImage<double>& input_image) const
{
    CcdBoxRegion CCD_overlap_region =
        CcdBoxRegion::overlapping_region(this->get_covering_region(),
                                          input_image.CCD_grid.image_outline_);

    PixelBoxRegion image_pix_overlap_region =
        input_image.CCD_grid.corresponding_pixel_region(CCD_overlap_region);
    image_pix_overlap_region.expand_to_pixel_boundaries();

    CcdImage<double> mask;
    mask.pix = PixelArray2d<double>(image_pix_overlap_region.x_dim(),
                                    image_pix_overlap_region.y_dim(), 0.0);
    mask.initialize_CCD_grid_to_specific_region(
        input_image.CCD_grid.corresponding_grid_region(image_pix_overlap_region));

    PixelShift mask_offset = image_pix_overlap_region.low - PixelPosition::origin;

    for (PixelIterator it(mask.pix.range()); it!=it.end; ++it) {
        CcdPosition input_pixel_centre =
            input_image.CCD_grid.corresponding_grid_Position(
                PixelPosition::centre_of_pixel(it) + mask_offset);
        if (this->contains_point(input_pixel_centre)) {
            mask.pix(it)=1.0;
        }
    }
    return mask;
}

template class CircularAperture<coordinate_types::pixels>;
template class CircularAperture<coordinate_types::ccd>;

//============================================================================


}//end namespace coela::regions
}//end namespace coela
