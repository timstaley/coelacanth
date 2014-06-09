/*
 * File:   pixel_array2d.h
 * Author: ts337
 *
 * Created on 09 February 2011, 10:32
 */

#ifndef COELA_PIXEL_ARRAY_H
#define COELA_PIXEL_ARRAY_H


#include "PixelArrayHeader.h"
#include "array_compression_info.h"


#include "coela_utility/src/string_utils.h"
#include <algorithm>
//#include <functional>
#include <numeric>
#include <stdexcept>

namespace coela {
//=========================================================================================================
///Pixel: Just tacks a value onto a pixel index. Useful for returning e.g. max pixel from a pixel_array

template<typename T>
struct Pixel: public PixelIndex {
    Pixel() {}
    Pixel(PixelIndex c, T v=0.0):PixelIndex(c), value(v) {}
    T value;
//    static bool descending_value_predicate(const pixel& first, const pixel& second);
};

//=========================================================================================================
/// PixelIterator : Incrementable pixel index covering a pixel_array, in row major ordering.
class PixelIterator: public PixelIndex {
public:
    PixelIterator(const PixelRange& rgn);
    PixelIterator& operator++(); //prefix
    PixelIndex operator++(int); //postfix
    const PixelIndex end;

private:
    PixelRange range;
};
//=======================================================================================================

template<typename T>
class PixelArray2d {
    //-----------------------------------------------------------------
    //friend declaration allows private access between template types (necessary for e.g. converting type)
    template<typename U> friend class PixelArray2d;

public:
    //-----------------------------------------------------------------
    //Public data members:
    //-----------------------------------------------------------------
    PixelArray2dHeader hdr;

    //-----------------------------------------------------------------
    //Constructors
    //-----------------------------------------------------------------
    PixelArray2d() {}   ///< Uninitialized:
    ///Construct from parameters
    PixelArray2d(const size_t x_dim_init, const size_t y_dim_init, T init_value);

    ///Promotion from any other data type to double data type.
    ///(this is unit tested only for float->double currently, but should work fine, uses std::copy on the pixel values).
    template<typename U>
    explicit PixelArray2d(const PixelArray2d<U>& rhs);

    ///subarray constructor
    static PixelArray2d sub_array(const PixelArray2d& , const PixelRange&);

    ///Since it is used in some compression schemes, we provide functionality
    ///to enable loading of an image stored in what is effectively a tarball -
    ///hence the "byte_offset" parameter.
    static PixelArray2d load_from_file(const string& filename,
                                       const size_t hdr_begin_byte_offset=0);

    ///Load from buffered file
    static PixelArray2d load_from_buffer(const FileBuffer&,
                                         const size_t hdr_begin_byte_offset=0);

    void write_to_file(const string& filename,
                       const ArrayCompressionInfo compression_inf=
                           ArrayCompressionInfo::no_compression()) const;

    //-----------------------------------------------------------------
    //Methods
    //-----------------------------------------------------------------
    fits_header_conventions::fits_imagetype bitpix() const {
        check_header_type_matches_array_type();
        return hdr.bitpix();
    }

    //pixel access
    PixelRange range() const {return hdr.range();}
    size_t n_pixels() const {return hdr.n_pixels();}


    std::vector<T>& get_raw_data_unsafe() {return data_;} ///< Use at own risk!
    const std::vector<T>& get_raw_data() const {return data_;}

    void assign(T value); ///<sets all pixels to this value.

    PixelIndex min_PixelIndex() const; ///<get index of pixel in img with min value
    PixelIndex max_PixelIndex() const; ///<get index of pixel in img with max value
    ///get index of pixel in specified pixel range with max value
    PixelIndex max_PixelIndex(const PixelRange&) const;

    T min_val() const { return *std::min_element(data_.begin(), data_.end()); }
    T max_val() const { return *std::max_element(data_.begin(), data_.end()); }

    T sum() const { return std::accumulate(data_.begin(), data_.end(), T(0));} //NB must init with "T(0)" not just "0" otherwise always call <int> version.
    //... it's stuff like this that make you glad of unit testing!

    T region_sum(const PixelRange&) const;

    ///Pixel access operators:
    ///NB x,y, vary from [1,x/ydim] inclusive in the style of DS9 pixel indexing.
    /// ( one-based indexing ).
    ///Essential a pixel is indexed by the cartesian co-ordinates of its upper corner
    const T& operator()(const int x, const int y) const;  //Read-only Access pixel(x,y)
    T & operator()(const int x,const  int y); //RW access pixel(x,y)

    const T& operator()(const PixelIndex& pi) const { return (*this)(pi.x,pi.y);} //Read-only Access pixel(x,y)
    T& operator()(const PixelIndex& pi) { return (*this)(pi.x,pi.y);}  //rw

    void operator*=(const PixelArray2d& mask);
    void operator/=(const PixelArray2d& divisor);
    void operator-=(const PixelArray2d& bias_frame);


    void operator+=(const PixelArray2d& summation_frame);

    template<typename U>
    void operator+=(const PixelArray2d<U>& summation_frame);

    void operator*=(const double factor)
    {for (typename std::vector<T>::iterator i=data_.begin(); i!=data_.end(); ++i) { (*i)*=factor; }}
    void operator/=(const double divisor)
    {for (typename std::vector<T>::iterator i=data_.begin(); i!=data_.end(); ++i) { (*i)/=divisor; }}
    void operator-=(const double subtraction)
    {for (typename std::vector<T>::iterator i=data_.begin(); i!=data_.end(); ++i) { (*i)-=subtraction; }}
    void operator+=(const double addition)
    {for (typename std::vector<T>::iterator i=data_.begin(); i!=data_.end(); ++i) { (*i)+=addition; }}

    bool operator==(const PixelArray2d& rhs) const;

    //Operators that return a new object, protected so we don't initialize a derived class badly in client code.
    PixelArray2d<double> operator-(const double subtractor);

    void check_header_type_matches_array_type() const; //If not, throw an exception

//-----------------------------------------------------------------------------------------------------------
public:
    ///The serialization functions for the header and array data are decoupled,
    ///so that extra stuff can be put into / read from the header
    ///after the essential "PixelArrayHeader" params.



    ///The serialization functions for the header and array data are decoupled
    /// so that extra stuff can be put into / read from the header after the essential "PixelArrayHeader" params.

    ///NB extra parameter allows you to specify compression level to encode in header table.

    void write_header_to_fht(FitsHeader& fht,
                             const ArrayCompressionInfo = ArrayCompressionInfo::no_compression()) const;

    std::istream& seek_and_load_data_from_stream(std::istream& is,
            const size_t data_begin_byte_offset,
            const ArrayCompressionInfo&);

    std::istream& load_data_from_positioned_stream(
        std::istream& is,
        const ArrayCompressionInfo& compression_inf);


    std::ostream& write_data_to_stream(std::ostream& os,
                                       const ArrayCompressionInfo&) const;



    void load_data_from_buffer(const FileBuffer&,
                               const size_t data_begin_byte_offset,
                               const ArrayCompressionInfo&);



//-----------------------------------------------------------------------------
private:
    //-----------------------------------------------------------
    //Data members:
    std::vector<T> data_;
    //-----------------------------------------------------------
    ///Construct and allocate array according to a given header, but do not initialize the array:
    PixelArray2d(const PixelArray2dHeader& hdr);

    void resize_data_vec_according_to_header_requirements();///< Always deletes data pointer first:

    PixelIndex PixelIndex_of_data_element(size_t i) const;

    //as otherwise we may be allocating the wrong kind of array for say, data loaded from a file.

    static fits_header_conventions::fits_imagetype array_FITS_imagetype();

    static size_t calculate_minimum_data_size_in_bytes(
        size_t n_pixels,const ArrayCompressionInfo&);

    //---------------------------------------------------------------------------------------------------------
    //load / write subroutines for the various formats:

    void load_array_from_regular_FITS_buffer(const FileBuffer& buf,
            const size_t data_begin_offset);

    std::istream& load_array_from_regular_FITS_stream(std::istream& is);
    std::ostream& write_array_to_regular_FITS_stream(std::ostream& os) const;

    void load_array_from_lcz_bits_buffer(
        const char * data_begin_ptr, const char * buffer_end_ptr,
        const ArrayCompressionInfo& compression_info
    );

    std::istream& load_array_from_lcz_bits_stream(std::istream& is,
            const ArrayCompressionInfo& compression_info);
//    std::ostream& write_array_to_lcz_bits_stream(std::ostream& os) const; //to do!

    void load_array_from_lcc_byte_buffer(
        const char * data_begin_ptr, const char * buffer_end_ptr);

    std::istream& load_array_from_lcc_byte_stream(std::istream& is);
};

//=========================================================================================================



}//end namespace coela

#include "implementation/pixel_array2d.icc"

#endif  /* PIXEL_ARRAY_H */

