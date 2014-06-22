#include "../image_grid.h"
#include "coela_utility/src/string_utils.h"

using namespace std;
namespace coela {
using namespace string_utils;

const double ds9_display_offset = 0.5;

//================================================================================
template<class coord_type>
ImageGrid<coord_type>::ImageGrid(const RectangularRegion<coord_type> region_init,
                                 double pixel_scale):
    image_outline_(region_init), pixel_width_(pixel_scale)
//    , grid_initialized_(true)
{}

template<class coord_type>
ImageGrid<coord_type>::ImageGrid(const PixelRange& range,
                                 const Position<coord_type> origin_posn,
                                 const double pixel_scale):
    pixel_width_(pixel_scale)
{
    PixelShift pixarray_diagonal(range.x_dim(), range.y_dim());
    TwoVector<coord_type> grid_diagonal(pixarray_diagonal.x * pixel_width_,
                                        pixarray_diagonal.y*pixel_width_);
    image_outline_=RectangularRegion<coord_type>(origin_posn,
                   origin_posn+grid_diagonal);
//    , grid_initialized_(true)
}

template<class coord_type>
bool ImageGrid<coord_type>::is_initialized() const
{
    if (pixel_width_==0.0) {
        return false;
    }
    return true;
}

template<>
void ImageGrid<coordinate_types::pixels>::write_to_fht(FitsHeader& fht) const
{
    if (!is_initialized()) {
        return;
    }
    //Set keywords so that image coords are matched to DS9 readout.
    fht.set_keyword("LTV1",  string_utils::ftoa(image_outline_.low.x +
                    ds9_display_offset*pixel_width_));
    fht.set_keyword("LTV2", string_utils::ftoa(image_outline_.low.y +
                    ds9_display_offset*pixel_width_));
    fht.set_keyword("LTM1_1", string_utils::ftoa(1.0));
    fht.set_keyword("LTM2_2", string_utils::ftoa(1.0));
}



///NB The idea of the ImageGrid<pixels> is simply to set the DS9 coords displaying correctly.
///The actual pixel region is defined entirely by the PixelArrayHeader and so we do not duplicate the information here.
template<>
void ImageGrid<coordinate_types::pixels>::load_from_fht(const FitsHeader&)
{
    pixel_width_ = 1.0;
}


template<>
void ImageGrid<coordinate_types::ccd>::write_to_fht(FitsHeader& fht) const
{
    if (!is_initialized()) {
        return;
    }
    fht.set_keyword("CCD_RGN", "T");
    //Set keywords so that image coords are matched to DS9 readout.
    fht.set_keyword("ATM1_1", string_utils::ftoa(pixel_width_));
    fht.set_keyword("ATM2_2", string_utils::ftoa(pixel_width_));
    fht.set_keyword("ATV1",  string_utils::ftoa(image_outline_.low.x));
    fht.set_keyword("ATV2", string_utils::ftoa(image_outline_.low.y));

    //Set keywords for easy loading:
    fht.set_keyword("CRGN_XLO",string_utils::ftoa(image_outline_.low.x));
    fht.set_keyword("CRGN_XHI", string_utils::ftoa(image_outline_.high.x));
    fht.set_keyword("CRGN_YLO", string_utils::ftoa(image_outline_.low.y));
    fht.set_keyword("CRGN_YHI", string_utils::ftoa(image_outline_.high.y));
}

template<>
void ImageGrid<coordinate_types::ccd>::load_from_fht(const FitsHeader& fht)
{
    if (fht.key_exists("CCD_RGN")) {
        assert(fht.get_key_value("ATM1_1")==fht.get_key_value("ATM2_2"));
        pixel_width_ = atof(fht.get_key_value("ATM1_1"));

        image_outline_.low.x =  atof(fht.get_key_value("CRGN_XLO"));
        image_outline_.high.x =  atof(fht.get_key_value("CRGN_XHI"));
        image_outline_.low.y =  atof(fht.get_key_value("CRGN_YLO"));
        image_outline_.high.y =  atof(fht.get_key_value("CRGN_YHI"));
    }
}



template<>
void ImageGrid<coordinate_types::mosaic>::write_to_fht(FitsHeader& fht) const
{
    if (!is_initialized()) {
        return;
    }
    //Set keywords so that image coords are matched to DS9 readout.
    fht.set_keyword("DTM1_1", string_utils::ftoa(pixel_width_));
    fht.set_keyword("DTM2_2", string_utils::ftoa(pixel_width_));
    fht.set_keyword("DTV1", string_utils::ftoa(image_outline_.low.x));
    fht.set_keyword("DTV2", string_utils::ftoa(image_outline_.low.y));

    fht.set_keyword("DRGN_XLO", string_utils::ftoa(image_outline_.low.x));
    fht.set_keyword("DRGN_XHI", string_utils::ftoa(image_outline_.high.x));
    fht.set_keyword("DRGN_YLO", string_utils::ftoa(image_outline_.low.y));
    fht.set_keyword("DRGN_YHI", string_utils::ftoa(image_outline_.high.y));

}

template<>
void ImageGrid<coordinate_types::mosaic>::load_from_fht(const FitsHeader& fht)
{
    assert(fht.get_key_value("DTM1_1")==fht.get_key_value("DTM2_2"));
    pixel_width_ = atof(fht.get_key_value("DTM1_1"));

    image_outline_.low.x  = atof(fht.get_key_value("DRGN_XLO"));
    image_outline_.low.y  = atof(fht.get_key_value("DRGN_YLO"));
    image_outline_.high.x = atof(fht.get_key_value("DRGN_XHI"));
    image_outline_.high.y = atof(fht.get_key_value("DRGN_YHI"));
}

template<class coord_type>
PixelPosition ImageGrid<coord_type>::corresponding_pixel_Position(
    const Position<coord_type>& ref_posn) const
{
    return PixelPosition((ref_posn.x - image_outline_.low.x)/pixel_width_,
                         (ref_posn.y - image_outline_.low.y)/pixel_width_);
}

template<class coord_type>
PixelShift ImageGrid<coord_type>::corresponding_pixel_shift(const TwoVector<coord_type>&
        ref_shift) const
{
    return PixelShift(ref_shift.x / pixel_width_,
                      ref_shift.y/ pixel_width_);
}

template<class coord_type>
PixelBoxRegion ImageGrid<coord_type>::corresponding_pixel_region(
    const RectangularRegion<coord_type>& ref_rgn)const
{
    return PixelBoxRegion(corresponding_pixel_Position(ref_rgn.low),
                          corresponding_pixel_Position(ref_rgn.high));
}

template<class coord_type>
Position<coord_type> ImageGrid<coord_type>::corresponding_grid_Position(
    const PixelPosition& p) const
{
    return Position<coord_type>(p.x*pixel_width_ + image_outline_.low.x,
                                p.y*pixel_width_ + image_outline_.low.y);
}
template<class coord_type>
RectangularRegion<coord_type> ImageGrid<coord_type>::corresponding_grid_region(
    const PixelBoxRegion& r) const
{
    return RectangularRegion<coord_type>(corresponding_grid_Position(r.low),
                                         corresponding_grid_Position(r.high));
}



template class ImageGrid<coordinate_types::pixels>;
template class ImageGrid<coordinate_types::ccd>;
template class ImageGrid<coordinate_types::mosaic>;


//================================================================================

}//end namespace coela


