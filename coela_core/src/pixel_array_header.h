/*
 * File:   pixel_array_header.h
 * Author: ts337
 *
 * Created on 08 February 2011, 17:05
 */
#ifndef COELA_PIXEL_ARRAY_HEADER_H
#define COELA_PIXEL_ARRAY_HEADER_H

//#include "../level1/coords.h"

#include "fits_header_conventions.h"
#include "fits_header.h"
#include "pixel_index.h"


namespace coela {

//=========================================================================================================
//PixelArrayHeader:
//A very lightweight class to represent the shape of a pixel array.
//NB not templatized, as the data type is unknown until initialized from the fits header! (And it would be unnecessary in any case)
//However, we must store the data type so we can calculate how much space to allocate for it later on.
//So the data type is encoded via the fits_imagetype enumeration, as it is assumed this will be used with FITS files

//So, functionally it:
//  * Provides interaction with the FITS_header_table to pick out these essential parameters.
//  * Provide the information to know what kind of pixel_array to instantiate.
//  * Supplies pixel_regions for pixel iterators

class PixelArray2dHeader {
public:
    PixelArray2dHeader():x_dim_(0),y_dim_(0),n_pix(0),
        bitpix_((fits_header_conventions::fits_imagetype)0) {}
    PixelArray2dHeader(fits_header_conventions::fits_imagetype img_type, size_t x,
                       size_t y);  ///< Create from bits
    PixelArray2dHeader(const FitsHeader& fht) {load_from_fht(fht);}

//    PixelArrayHeader(const FitsHeader& fht){load_array_dims_and_type_from_fht(fht);}
    static PixelArray2dHeader subarray_header(const PixelArray2dHeader& original_array,
            const PixelRange& subarray_box);

    PixelRange range() const {return PixelRange(1,1,x_dim(),y_dim());}

    void write_to_fht(FitsHeader&,
                      const fits_header_conventions::fits_imagetype bitpix_to_set) const;
    void load_from_fht(const FitsHeader&);

    //Accessors
    size_t x_dim() const {return x_dim_;}
    size_t y_dim() const {return y_dim_;}
    size_t n_pixels() const {return n_pix;}
    size_t byte_size_of_uncompressed_data() const;

    ///Get the current bitpix, usually as loaded from a FitsHeader
    ///NB we reset the bitpix when writing to relect the requested write method.
    ///(see "write_to_fht" parameters)
    fits_header_conventions::fits_imagetype bitpix() const { return bitpix_; }

    bool operator==(const PixelArray2dHeader& rhs) const;

private:
    //-------------------------------------------------------------------------------------------------
    //Data
    //Basic FITS stuff
//   int naxes; //Must be 2, currently
    size_t x_dim_, y_dim_; //naxis1, naxis2
    size_t n_pix; //number of elements, sx*sy
    fits_header_conventions::fits_imagetype bitpix_;
    //-------------------------------------------------------------------------------------------------
};

//=========================================================================================================

//Helper class to load basic header info from a pixel cube; added "Z" dimension is used as an index
class PixelCubeHeader {
public:
    PixelCubeHeader() {}
    PixelCubeHeader(const FitsHeader& fht) {load_from_fht(fht);}

    void load_from_fht(const FitsHeader&);
    PixelArray2dHeader get_xy_array2d_header() const { return PixelArray2dHeader(bitpix_,x_dim_,y_dim_);}
    size_t z_dim() const {return z_dim_;}


private:
    //-------------------------------------------------------------------------------------------------
    //Data
    size_t x_dim_, y_dim_, z_dim_;
    fits_header_conventions::fits_imagetype bitpix_;
    //-------------------------------------------------------------------------------------------------
};


//=========================================================================================================
//struct multi_extension_header; //To do.
//=========================================================================================================
}//end namespace coela
#endif  /* PIXEL_ARRAY_HEADER_H */

