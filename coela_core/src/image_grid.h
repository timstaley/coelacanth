/*
 * File:   image_grid.h
 * Author: ts337
 *
 * Created on 14 February 2011, 16:51
 */

#ifndef COELA_IMAGE_REGION_H
#define COELA_IMAGE_REGION_H

#include "cartesian_coords.h"
#include "fits_header.h"

namespace coela {
//===============================================================
/**
A struct to represent a region of a co-ordinate grid.
The idea being that an 'Image' is a set of pixels associated with a particular
region of a co-ordinate space. So by tying together a PixelArray and an ImageGrid
we get a meaningful construct (an 'Image').

The pixel_scale of the image in this co-ordinate grid is recorded -
this allows conversion between pixel co-ordinates and reference system
coordinates, thus allowing e.g. a Gaussian fit to be located in the external
co-ordinate system.
(Square pixels assumed currently.)
*/

template<class coord_type>
struct ImageGrid {
public:
    //-----------------------------------------------------
    //data
    RectangularRegion<coord_type> image_outline_;
    double pixel_width_; ///<width of an image pixel in units of the grid coord system
    //-----------------------------------------------------

    ImageGrid():pixel_width_(0.0)/*, grid_initialized_(false)*/{}
    ImageGrid(const RectangularRegion<coord_type> region, double pixel_scale);
    ImageGrid(const PixelRange& range, const Position<coord_type> origin_posn,
              double pixel_scale);

    bool is_initialized() const;

    void write_to_fht(FitsHeader&) const;
    void load_from_fht(const FitsHeader&);


    PixelPosition corresponding_pixel_Position(const Position<coord_type>&) const;
    PixelShift corresponding_pixel_shift(const TwoVector<coord_type>&) const;
    PixelBoxRegion corresponding_pixel_region(const RectangularRegion<coord_type>&)const;

    Position<coord_type> corresponding_grid_Position(const PixelPosition&) const;
    RectangularRegion<coord_type> corresponding_grid_region(const PixelBoxRegion&) const;

//    bool operator==(const ImageGrid& rhs) const{
//          if (rhs.image_outline_ == image_outline_ && rhs.pixel_width_ ==pixel_width_ ) return true; else return false;}

    //To do -- implement own state tracking (convert to class)
//private:
//    bool grid_initialized_;
};

//================================================================
}//end namespace coela

#endif  /* IMAGE_REGION_H */

