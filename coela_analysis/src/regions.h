/*
 * File:   regions.h
 * Author: ts337
 *
 * Created on 16 May 2011, 17:24
 */

#ifndef COELA_REGIONS_H
#define COELA_REGIONS_H

#include "coela_core/src/cartesian_coords.h"
#include "coela_core/src/ccd_image.h"
namespace coela {
namespace regions {

//============================================================================
template<class coord_type>
class region_interface {
public:
    virtual ~region_interface() {}
    virtual bool contains_point(const Position<coord_type>& testpoint) const=0;
    virtual RectangularRegion<coord_type> get_covering_region() const=0;
    virtual double analytic_area() const=0; //NB in units of co-ordinate system!
    virtual CCDImage<double> generate_mask_for_covered_portion_of_image(
        const CCDImage<double>& image_to_mask) const=0;

//        ds9::region DS9Region(); //implemented for some but not yet all
};
//============================================================================

//Handy generic subroutines:
//CCD_Position get_region_centroid(const image<double>& image, const region_base<CCD_Position>& rgn, const double threshold=0.0);
//double get_region_weight(const float_bitmap& image, const region_base<frame_Position>& frm_rgn);

//Internal generic subroutine:
//image<double> generate_mask_for_image(const ImageGrid<coord_type>& image_oultin)
//============================================================================

template <class coord_type>
class circular_aperture: public region_interface<coord_type> {
public:
    //constructor
    circular_aperture(const Position<coord_type>& centre_init,
                      const double radius_in_these_coords);

//    circular_aperture(const ds9::region& rgn):region_base<Coord_Type>(rgn.coord_type){
//        if (rgn.region_type!="circle") throw std::runtime_error("Cannot convert this region to a circular aperture");
//        assert(rgn.parameters.size()==3);
//        centre = Coord_Type(rgn.parameters[0], rgn.parameters[1]);
//        pixel_radius = rgn.parameters[2];
//        initialise_cover_region();
//    }

    bool contains_point(const Position<coord_type>& testpoint) const;
    RectangularRegion<coord_type> get_covering_region() const;
    double analytic_area() const;
    CCDImage<double> generate_mask_for_covered_portion_of_image(const CCDImage<double>&
            image_to_mask) const;

//    ds9::region DS9Region() const;
    Position<coord_type> centre;
    double radius;


};




}//end namespace coela::regions
}//end namespace coela

#endif  /* REGIONS_H */

