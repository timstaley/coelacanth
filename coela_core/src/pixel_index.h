/*
 * File:   PixelIndex.h
 * Author: ts337
 *
 * Created on 10 February 2011, 12:08
 */

#ifndef COELA_PIXEL_INDEX_H
#define COELA_PIXEL_INDEX_H

#include <istream>
#include <ostream>

#ifndef INT_MAX
#include <limits>
#endif

namespace coela {

//=========================================================================================================
//PixelIndex:
/**
 A simple struct for hanging an x and a y index value together.

 Has automatic stream readout formatting and equality testing.
 NB *signed* integers are used to represent the index values, since we may want
 to use a negative pixel index to indicate an offset.
 This is generally not a problem since even a large image will rarely have an
 x or y dimension large enough to require unsigned values,
 even though the product x*y might well exceed that limit.
 In extreme cases the int could of course be bumped up to a long int with no problem.
 */

struct PixelIndex {

//---------------------------------
//Data
    int x, y;
//---------------------------------
    //default values should hopefully throw errors if used accidentally
    PixelIndex() :
            x(std::numeric_limits<int>::max()), y(
                    std::numeric_limits<int>::max()) {
    }
    PixelIndex(int x_in, int y_in) :
            x(x_in), y(y_in) {
    }

//    //Named constructors for various derivations from co-ordinate Positions:
//    static PixelIndex pixel_containing_point(const frame_Position&);
//    static PixelIndex pixel_centred_at(const PixelPosition&);

    bool is_valid() const {
        if (x == std::numeric_limits<int>::max()
                || y == std::numeric_limits<int>::max()) {
            return false;
        } else {
            return true;
        }
    }

    bool operator==(const PixelIndex& rhs) const {
        if (rhs.x == x && rhs.y == y) {
            return true;
        } else {
            return false;
        }
    }
    bool operator!=(const PixelIndex& rhs) const {
        if (rhs.x != x || rhs.y != y) {
            return true;
        } else {
            return false;
        }
    }

    PixelIndex& operator -=(const PixelIndex& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    PixelIndex& operator +=(const PixelIndex& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

};
inline std::ostream& operator<<(std::ostream& os, const PixelIndex& pi) {
    os << "(" << pi.x << "," << pi.y << ")";
    return os;
}
inline std::istream& operator>>(std::istream& is, PixelIndex& pi) {
    char c;
    is >> c >> pi.x >> c >> pi.y >> c;
    return is;
}
PixelIndex operator-(const PixelIndex& lhs, const PixelIndex& rhs);
PixelIndex operator+(const PixelIndex& lhs, const PixelIndex& rhs);

//=========================================================================================================
///Defines a selection of pixels from [ (lowx,lowy),(highx,highy) ] inclusive.
struct PixelRange {
    //---------------------------------
    //Data
    PixelIndex low, high; //define the corners of the box
    //---------------------------------

    PixelRange() {
    }
    PixelRange(const PixelIndex& low_in, const PixelIndex& high_in) :
            low(low_in), high(high_in) {
    }
    PixelRange(int xlow, int ylow, int xhigh, int yhigh) :
            low(xlow, ylow), high(xhigh, yhigh) {
    }

    ///Get an enlarged / shrunk copy of a box - i.e. a copy padded by an integer width in all directions
    //+ve val to grow the region, -ve val to shrink it
    static PixelRange pad(const PixelRange&, int pixel_pad_width);

    ///Get overlapping region - if no overlap an uninitialized pixel_box() is returned
    static PixelRange overlap(const PixelRange& box1, const PixelRange& rgn2);

    ///Get copy enlarged to cover the given pixel
    static PixelRange stretch_to_pixel(const PixelRange& original_box,
            const PixelIndex&);

    size_t x_dim() const {
        return high.x - low.x + 1;
    }
    size_t y_dim() const {
        return high.y - low.y + 1;
    }
    size_t n_pix() const {
        return x_dim() * y_dim();
    }

    bool is_valid() const {
        if (high.x >= low.x && high.y >= low.y && high.is_valid()
                && low.is_valid()) {
            return true;
        } else {
            return false;
        }
    }
    bool contains_range(const PixelRange& testbox) const;
    bool contains_pixel(const PixelIndex&) const;

    //Where data index ranges from 0 to N-1 (for N pixels in the range)
    //TO DO : Move this to become a member of PixelArray2d.
    PixelIndex get_PixelIndex_for_data_vector_element(
            const size_t data_vec_index) const;

    bool operator==(const PixelRange& rhs) const {
        if (rhs.low == low && rhs.high == high) {
            return true;
        } else {
            return false;
        }
    }
    bool operator!=(const PixelRange& rhs) const {
        if (*this == rhs) {
            return false;
        }
        return true;
    }
//    static pixel_box bounded_by(const frame_box_region& rgn);///<Conversion operator, checks that the box boundaries are at integer values and returns the corresponding pixel range

};
inline std::ostream& operator<<(std::ostream& os, const PixelRange& pb) {
    os << "[" << pb.low << "," << pb.high << "]";
    return os;
}
inline std::istream& operator>>(std::istream& is, PixelRange& pb) {
    char c;
    is >> c >> pb.low >> c >> pb.high >> c;
    return is;
}

//=========================================================================================================

}//end namespace coela

#endif  /* PIXEL_INDEX_H */

