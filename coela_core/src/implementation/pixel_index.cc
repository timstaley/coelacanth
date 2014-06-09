#include "../pixel_index.h"
#include <algorithm>
#include <cassert>

using std::max; using std::min;
namespace coela {

PixelRange PixelRange::pad(const PixelRange& orig, int npix)
{
    return  PixelRange(orig.low.x - npix, orig.low.y-npix, orig.high.x+npix,
                       orig.high.y+npix);
}

PixelRange PixelRange::overlap(const PixelRange& box1, const PixelRange& box2)
{
    //check there is some overlap, if not return uninitialized pixel_box;
    if (box1.low.x > box2.high.x || box1.high.x < box2.low.x
            || box1.low.y > box2.high.y || box1.high.y < box2.low.y) {
        return PixelRange();
    }

    return PixelRange(max(box1.low.x, box2.low.x),  max(box1.low.y, box2.low.y),
                      min(box1.high.x, box2.high.x), min(box1.high.y, box2.high.y)
                     );
}

PixelRange PixelRange::stretch_to_pixel(const PixelRange& box, const PixelIndex& pix)
{
    return PixelRange(
               min(box.low.x, pix.x), min(box.low.y, pix.y),
               max(box.high.x, pix.x), max(box.high.y, pix.y)
           );
}

bool PixelRange::contains_range(const PixelRange& testbox) const
{
    assert(is_valid());
    assert(testbox.is_valid());
    if (testbox.high.x<=high.x && testbox.high.y<=high.y  && testbox.low.x >=low.x
            && testbox.low.y >=low.y) { return true; }
    return false;
}

bool PixelRange::contains_pixel(const PixelIndex& p) const
{
    assert(is_valid());
    assert(p.is_valid());
    if (p.x >= low.x  && p.y >=low.y &&
            p.x <= high.x && p.y <=high.y) { return true; }
    return false;
}

PixelIndex PixelRange::get_PixelIndex_for_data_vector_element(
    const size_t datavec_index) const
{
    //Equations incorrect unless mapping to a window covering whole data array.
    assert(low.x ==1 && low.y ==1);
    //(Of course, if the high pixel is off then we're also up the creek, but this is unknowable without access to the pixel_array.
//    assert(datavec_index>=0); //Does nothing since we're taking a size_t
    assert(datavec_index<n_pix());
    return PixelIndex((datavec_index)% x_dim() +1, (datavec_index) / x_dim() +1);
}

PixelIndex operator-(const PixelIndex& lhs,const PixelIndex& rhs)
{
    PixelIndex temp=lhs;
    temp-=rhs;
    return temp;
}
PixelIndex operator+(const PixelIndex& lhs,const PixelIndex& rhs)
{
    PixelIndex temp=lhs;
    temp+=rhs;
    return temp;
}

//pixel_box pixel_box::bounded_by(const frame_box_region& rgn){
//    assert(rgn.low.x >=0.0);
//    assert(rgn.low.y >=0.0);
//    assert( fmod(rgn.low.x, 1.0)==0.0);
//    assert( fmod(rgn.low.y, 1.0)==0.0);
//    assert( fmod(rgn.high.x, 1.0)==0.0);
//    assert( fmod(rgn.high.y, 1.0)==0.0);
//    PixelIndex low(rgn.low.x+1.0, rgn.low.y+1.0), high(rgn.high.x, rgn.high.y);
//    return pixel_box(low,high);
//}




}//end namespace coela
