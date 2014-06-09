/*
 * File:   cartesian_coords.cc
 * Author: ts337
 *
 * Created on 10 February 2011, 12:02
 */



#include "../cartesian_coords.h"

#include <cmath>
#include <map>

namespace coela {

namespace coordinate_types {
std::vector<std::string> all_ds9_strings()
{
    using std::vector; using std::string;
    vector<string> v;
    v.push_back(pixels::ds9_string());
    v.push_back(CCD::ds9_string());
    v.push_back(mosaic::ds9_string());
    v.push_back(wcs::ds9_string());
    return v;
}
}

const double position_precision_limit=1e-7;
//=======================================================================================================
//
// Position
//
template<>
Position<coordinate_types::pixels> Position<coordinate_types::pixels>::centre_of_pixel(
    const PixelIndex& pi)
{
    Position<coordinate_types::pixels> p(pi.x-0.5, pi.y-0.5);
    return p;
}

template<>
PixelIndex Position<coordinate_types::pixels>::pixel_centred_at(
    const Position<coordinate_types::pixels>& p)
{
    assert(fabs(fmod(p.x, 1.0)-0.5) < position_precision_limit);
    assert(fabs(fmod(p.y, 1.0)-0.5) < position_precision_limit);
    return PixelIndex(p.x+0.5, p.y+0.5);
}
template<>
PixelIndex Position<coordinate_types::pixels>::pixel_containing_point(
    const Position<coordinate_types::pixels>& p)
{
    return PixelIndex(ceil(p.x), ceil(p.y));
}

//=======================================================================================================
//
//Rectangular region
//
template<class ref_frame>
Position<ref_frame> RectangularRegion<ref_frame>::centre() const
{
    return Position<ref_frame>::mid_point(high,low) ;
}

template<class ref_frame>
bool RectangularRegion<ref_frame>::contains_point(const Position<ref_frame>& pt) const
{
    if (pt.x >= low.x && pt.y >= low.y
            &&
            pt.x <= high.x && pt.y <= high.y) { return true; }
    else { return false; }
}

template<class ref_frame>
bool RectangularRegion<ref_frame>::contains_region(const RectangularRegion<ref_frame>&
        testbox) const
{
    if (this->contains_point(testbox.low)
            && this->contains_point(testbox.high)) { return true; }
    else { return false; }
}

template<class ref_frame>
RectangularRegion<ref_frame>&
RectangularRegion<ref_frame>::shrink_to_pixel_boundaries()
{
    assert(is_valid());
    *this =  RectangularRegion<ref_frame>(ceil(low.x), ceil(low.y),
                                          floor(high.x), floor(high.y));
    return *this;
}

template<class ref_frame>
RectangularRegion<ref_frame>&
RectangularRegion<ref_frame>::expand_to_pixel_boundaries()
{
    assert(is_valid());
    *this= RectangularRegion<ref_frame>(floor(low.x), floor(low.y),
                                        ceil(high.x), ceil(high.y));
    return *this;
}

template<class ref_frame>
RectangularRegion<ref_frame>&
RectangularRegion<ref_frame>::enlarge_to_cover_point(const Position<ref_frame>&
        point_to_cover)
{
    *this = RectangularRegion<ref_frame>(
                std::min(low.x, point_to_cover.x),std::min(low.y, point_to_cover.y),
                std::max(high.x, point_to_cover.x),std::max(high.y, point_to_cover.y));
    return *this;
}

template<class ref_frame>
RectangularRegion<ref_frame>&
RectangularRegion<ref_frame>::enlarge_to_cover_region(
    const RectangularRegion<ref_frame>& rgn_to_cover)
{
    enlarge_to_cover_point(rgn_to_cover.low);
    enlarge_to_cover_point(rgn_to_cover.high);
    return *this;
}
template<class ref_frame>
RectangularRegion<ref_frame>
RectangularRegion<ref_frame>::padded_copy(const RectangularRegion<ref_frame>& orig,
        double pad_width)
{
    return RectangularRegion<ref_frame>(orig.low.x -pad_width, orig.low.y-pad_width,
                                        orig.high.x +pad_width, orig.high.y+pad_width);
}

template<class ref_frame>
RectangularRegion<ref_frame>
RectangularRegion<ref_frame>::overlapping_region(const RectangularRegion<ref_frame>& rgn1,
        const RectangularRegion<ref_frame>& rgn2)
{
    using std::max; using std::min;
    return RectangularRegion<ref_frame>(
               max(rgn1.low.x, rgn2.low.x) ,   max(rgn1.low.y, rgn2.low.y) ,
               min(rgn1.high.x, rgn2.high.x) ,   min(rgn1.high.y, rgn2.high.y));
}

template<class ref_frame>
bool RectangularRegion<ref_frame>::is_valid() const
{
    if (high.x>=low.x && high.y>=low.y && high.is_valid() && low.is_valid()) { return true; }
    return false;
}

//----------------------------------------------------------------------------------------
//Only defined for coord-type == pixels:
template<>
PixelRange PixelBoxRegion::bounded_pixels() const
{
    //Check that this PixelBoxRegion outline actually does lie along valid pixel boundaries:
    assert(fmod(low.x, 1.0)==0.0);
    assert(fmod(low.y, 1.0)==0.0);
    assert(fmod(high.x, 1.0)==0.0);
    assert(fmod(high.y, 1.0)==0.0);
    return PixelRange(low.x+1.0, low.y+1.0,
                      high.x, high.y);
}

template<>
PixelBoxRegion PixelBoxRegion::pixel_box_outline(const PixelRange& box)
{
    return PixelBoxRegion(box.low.x-1.0, box.low.y-1.0, box.high.x, box.high.y);
}
//

//=======================================================================================================

template class RectangularRegion<coordinate_types::pixels>;
template class RectangularRegion<coordinate_types::CCD>;
template class RectangularRegion<coordinate_types::mosaic>;


}//end namespace coela