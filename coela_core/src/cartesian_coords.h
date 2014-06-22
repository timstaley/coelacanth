/*
 * File:   cartesian_coords.h
 * Author: ts337
 *
 * Created on 10 February 2011, 12:02
 */

#ifndef COELA_CARTESIAN_COORDS_H
#define COELA_CARTESIAN_COORDS_H

#include<string>
#include<vector>
#include<limits>
#include <ostream>
#include <istream>
#include "pixel_index.h"

#include <cassert>
#include <cmath>

using std::string;


namespace coela {

//Fwd declaration:
template<class ref_frame>
struct TwoVector;

//NB See end of file for shorthand typedefs.
//=======================================================================================================
//=======================================================================================================
namespace coordinate_types {
//Dummy or "tag" classes to denote separate frames of reference for use in templatized cartesian coordinate classes:

///This refers to the subpixel co-ordinates of whatever pixel array we are currently indexing into.
///Cartesian convention, so the lower-left corner of image is (0,0).
struct pixels {
    static string ds9_string() {return "physical";}
};

///This refers to the subpixel co-ordinates of the original CCD array from which an image is derived.
///Of course, if there are multiple CCDs this will need combining with the CCD id, or... see below.
struct ccd {
    static  string ds9_string() {return "amplifier";}
};

///This refers to an arbitrary co-ordinate system spanning a mosaic of multiple CCDs
struct mosaic {
    static  string ds9_string() {return "detector";}
};

//Refers to the world co-ordinate system (RA,Dec)
struct wcs {
    static  string ds9_string() {return "fk5";}
};

std::vector<std::string> all_ds9_strings();

}//end namespace coela::coordinate_types
//=======================================================================================================
//=======================================================================================================
///Origin - a helper class to optimize conversion from vectors to Positions and vice versa.
/// Explanation:
/// pos(x,y) = vec(x,y) + pos (0,0) , i.e. a vector and a Position have the same components if the vector is measured relative to the origin.
/// But we still wish to differentiate the concepts, as this results in clearer statement of intentions in the code.
/// We use origin instead of Position(0,0) since this allows the compiler to use optimized conversions. (And arguably results in clearer code.)
struct origin_helper_class {
};

//=======================================================================================================
///A Position in a 2d plane.
template<class ref_frame>
struct Position {
public:

    //------------------------------------------------------------------------------------------
    //Data members:
    double x,y;
    //------------------------------------------------------------------------------------------

    ///Default init should flag errors early if used accidentally ( see also is_valid() )
    Position():x(NAN),y(NAN) {}
    Position(double x_, double y_):x(x_),y(y_) {}

    ///NB Only performs correctly if the -ffast-math compilation optimizations are off:
    bool is_valid() const { return !(std::isnan(x) || std::isnan(y));}


    //------------------------------------------------------------------------------------------
    //Functions specific to certain types:

    ///We enforce that this function is only valid for the pixels system,
    ///to remind the user that pixel indices always refer to the *current* array
    ///and not the original co-ordinate Position...
    static Position<coordinate_types::pixels> centre_of_pixel(const PixelIndex&);
    static PixelIndex pixel_centred_at(const Position<coordinate_types::pixels>&);
    static PixelIndex pixel_containing_point(const Position<coordinate_types::pixels>&);


    //------------------------------------------------------------------------------------------
    static Position<ref_frame> mid_point(const Position& lhs,const Position& rhs) {
        return Position((lhs.x+rhs.x) /2.0, (lhs.y+rhs.y)/2.0);
    }

    ///This is used for generating .reg files that can be loaded into the DS9 fits viewer.
    static string ds9_reference_frame_string() {return ref_frame::ds9_string();}

    bool operator==(const Position& rhs) const {return (x==rhs.x && y==rhs.y);}

    Position& operator +=(const TwoVector<ref_frame>& rhs) {x+=rhs.x; y+=rhs.y; return *this;}
    Position& operator -=(const TwoVector<ref_frame>& rhs) {x-=rhs.x; y-=rhs.y; return *this;}
    //
    static origin_helper_class origin;


};

template<class T>
origin_helper_class Position<T>::origin;  ///Instantiate the static class member


template<class ref_frame>
std::ostream& operator<<(std::ostream& os, const Position<ref_frame>& posn)
{
    os <<"(" <<posn.x<<","<<posn.y<<")"; return os;
}

template<class ref_frame>
std::istream& operator>>(std::istream& is, Position<ref_frame>& posn)
{
    char c; is>>c>>posn.x>>c>>posn.y>>c; return is;
}


//=======================================================================================================

template<class T>
double coord_distance_squared(const Position<T>& p1, const Position<T>& p2)
{
    return (p1-p2).length_squared();
}

template<class T>
double coord_distance(const Position<T>& p1, const Position<T>& p2)
{
    return (p1-p2).length();
}

//=======================================================================================================
///A translation vector in 2d.
template<class ref_frame>
struct TwoVector {
public:
    //------------------------------------------------------------------------------------------
    //Data members:
    double x,y;
    //------------------------------------------------------------------------------------------
    ///default init should flag errors early if used accidentally ( see also is_valid() )
    TwoVector():x(NAN),y(NAN) {}
//    Position(){}
    TwoVector(double x_, double y_):x(x_),y(y_) {}

    double length_squared() const {return x*x + y*y;}
    double length() const {return sqrt(x*x + y*y);}

    bool is_valid() const { return !(std::isnan(x) || std::isnan(y));} ///<NB Only performs correctly if the -ffast-math compilation optimizations are off.
    bool operator==(const TwoVector& rhs) const {return (x==rhs.x && y==rhs.y);}

    TwoVector& operator /=(const double d) {x/=d; y/=d; return *this;}
    TwoVector& operator *=(const double factor) {x*=factor; y*=factor; return *this;}

    TwoVector& operator +=(const TwoVector& rhs) {x+=rhs.x; y+=rhs.y; return *this;}
    TwoVector& operator -=(const TwoVector& rhs) {x-=rhs.x; y-=rhs.y; return *this;}
};

template<class ref_frame>
std::ostream& operator<<(std::ostream& os, const TwoVector<ref_frame>& v) { os <<"(" <<v.x<<","<<v.y<<")"; return os;}

template<class ref_frame>
std::istream& operator>>(std::istream& is, TwoVector<ref_frame>& v) { char c; is>>c>>v.x>>c>>v.y>>c; return is;}

/// vec= vec + vec
template<class T>
TwoVector<T> operator+(const TwoVector<T>& lhs,const TwoVector<T>& rhs) { return TwoVector<T>(lhs.x+rhs.x, lhs.y+rhs.y);}

/// vec= vec - vec
template<class T>
TwoVector<T> operator-(const TwoVector<T>& lhs,const TwoVector<T>& rhs) { return TwoVector<T>(lhs.x-rhs.x, lhs.y-rhs.y);}

/// vec = vec / scalar
template<class T>
TwoVector<T> operator/(const TwoVector<T>& lhs,const double rhs) {TwoVector<T> p(lhs); p/=rhs; return p;}

/// vec = vec * scalar
template<class T>
TwoVector<T> operator*(const TwoVector<T>& lhs,const double rhs) {TwoVector<T> p(lhs); p*=rhs; return p;}

//=======================================================================================================
///Addition and subtraction operators between Positions and vectors:


/// ( conversion ) vec = posn - origin
template<class T>
TwoVector<T> operator-(const Position<T>& lhs,const origin_helper_class&) { return TwoVector<T>(lhs.x, lhs.y);}

/// ( conversion ) posn = vec + origin
template<class T>
Position<T> operator+(const TwoVector<T>& lhs,const origin_helper_class&) { return Position<T>(lhs.x, lhs.y);}

/// vec = posn - posn
template<class T>
TwoVector<T> operator-(const Position<T>& lhs,const Position<T>& rhs) { return TwoVector<T>(lhs.x-rhs.x, lhs.y - rhs.y);}

/// posn = posn + vec
template<class T>
Position<T> operator+(const Position<T>& lhs,const TwoVector<T>& rhs) { return Position<T>(lhs.x + rhs.x, lhs.y + rhs.y); }

/// posn = posn - vec
template<class T>
Position<T> operator-(const Position<T>& lhs,const TwoVector<T>& rhs) {  return Position<T>(lhs.x - rhs.x, lhs.y - rhs.y); }


//=======================================================================================================
///A region in a 2d plane.
template<class ref_frame>
struct RectangularRegion {
public:
    //---------------------------------------------------------------------
    //Data:
    Position<ref_frame> low, high; //define the corners of the box
    //---------------------------------------------------------------------


    RectangularRegion() {}

    RectangularRegion(double xlow, double ylow, double xhigh, double yhigh):
        low(xlow, ylow), high(xhigh, yhigh) { assert(is_valid()); }

    RectangularRegion(const Position<ref_frame>& low_corner,
                      const Position<ref_frame>& high_corner):
        low(low_corner), high(high_corner) {assert(is_valid());}

    static RectangularRegion padded_copy(const RectangularRegion<ref_frame>&,
                                         double pad_width);  //+ve val to grow the region, -ve val to shrink it

    ///NB unsafe; assumes that the regions DO have an overlap
    static RectangularRegion overlapping_region(const RectangularRegion<ref_frame>& rgn1,
            const RectangularRegion<ref_frame>& rgn2);

    static RectangularRegion<coordinate_types::pixels> pixel_box_outline(const PixelRange&);

    //Informative functions:
    double x_dim() const { return  high.x - low.x; } //Now the edges are infinitely thin (rather than 1 pixel thick)
    double y_dim() const { return  high.y - low.y; } //and therefore are bounds rather than indices (so no +1 here)
    double area() const {return x_dim()*y_dim(); }

    Position<ref_frame> centre() const;

    bool contains_point(const Position<ref_frame>& pt) const;
    bool contains_region(const RectangularRegion<ref_frame>& testbox) const;

    ///Only defined for coord-type==pixels:
    PixelRange bounded_pixels() const;

    //----------------------------------------------------------------------------------------------
    //Modifier functions:
    RectangularRegion& shrink_to_pixel_boundaries() ;
    RectangularRegion& expand_to_pixel_boundaries() ;

    RectangularRegion<ref_frame>& enlarge_to_cover_point(const Position<ref_frame>&) ;
    RectangularRegion<ref_frame>& enlarge_to_cover_region(const RectangularRegion<ref_frame>&)
    ;

    //----------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------------

    //Handy bits and pieces
    static string ds9_reference_frame_string() {return ref_frame::ds9_string();}
    bool is_valid()
    const; ///<NB Only performs correctly if the -ffast-math compilation optimizations are off.
    bool operator==(const RectangularRegion& rhs) const {
        if (rhs.low ==low && rhs.high==high) { return true; } else { return false; }
    }
    bool operator!=(const RectangularRegion& rhs) const {
        if (*this==rhs) { return false; } return true;
    }
};

template<class ref_frame>
std::ostream& operator<<(std::ostream& os, const RectangularRegion<ref_frame>& rgn)
{
    os<<'['<<rgn.low <<" , "<< rgn.high<<']'; return os;
}

template<class ref_frame>
std::istream& operator>>(std::istream& is, RectangularRegion<ref_frame>& rgn)
{
    char c; is>>c>>rgn.low >>c>> rgn.high >>c; return is;
}


//=======================================================================================================

//typedefs:

typedef Position<coordinate_types::pixels> PixelPosition ;
typedef Position<coordinate_types::ccd> CcdPosition ;
typedef Position<coordinate_types::mosaic> MosaicPosition ;
typedef Position<coordinate_types::wcs> WcsPosition ;

typedef TwoVector<coordinate_types::pixels> PixelShift ;
typedef TwoVector<coordinate_types::ccd> CcdPixelShift ;
typedef TwoVector<coordinate_types::mosaic> MosaicPixelShift ;
typedef Position<coordinate_types::wcs> WcsShift;

typedef RectangularRegion<coordinate_types::pixels> PixelBoxRegion;
typedef RectangularRegion<coordinate_types::ccd> CcdBoxRegion;
typedef RectangularRegion<coordinate_types::mosaic> MosaicBoxRegion;


}//end namespace coela

#endif  /* CARTESIAN_COORDS_H */

