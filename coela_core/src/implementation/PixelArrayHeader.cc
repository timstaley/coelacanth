/*
 * File:   PixelArrayHeader.cc
 * Author: ts337
 *
 * Created on 08 February 2011, 17:05
 */

#include "../PixelArrayHeader.h"
#include "coela_utility/src/string_utils.h"
#include <stdexcept>

namespace coela {
using namespace string_utils;

//=========================================================================================================
//PixelArrayHeader:

PixelArray2dHeader::PixelArray2dHeader(
    fits_header_conventions::fits_imagetype img_type,
    size_t x, size_t y):
//naxes(2),
    x_dim_(x), y_dim_(y), bitpix_(img_type)
{
    n_pix=x_dim()*y_dim();
}

PixelArray2dHeader PixelArray2dHeader::subarray_header(const PixelArray2dHeader&
        original_array, const PixelRange& subarray_box)
{
    if (!original_array.range().contains_range(subarray_box)) {
        throw std::runtime_error("subarray_header: pixel_box is not contained in original array");
    }
    return PixelArray2dHeader(original_array.bitpix(), subarray_box.x_dim(),
                              subarray_box.y_dim());
}


void PixelArray2dHeader::write_to_fht(FitsHeader& fht,
                                      const fits_header_conventions::fits_imagetype bitpix_to_set) const
{
    using string_utils::itoa;
    FitsHeader prefix_fht;
    prefix_fht.set_keyword("SIMPLE","T");
    prefix_fht.set_keyword("BITPIX",itoa((int)bitpix_to_set));
    prefix_fht.set_keyword("NAXIS",itoa(2)); //naxes
    prefix_fht.set_keyword("NAXIS1",itoa(x_dim_));
    prefix_fht.set_keyword("NAXIS2",itoa(y_dim_));
    prefix_fht.set_keyword("BSCALE","1.0");
    prefix_fht.set_keyword("BZERO","0.0");

    fht.merge_with_prefix_table(prefix_fht);
    fht.remove_key("NAXIS3"); //in case we loaded from a fits cube.
//    if (bitpix==FLOATIMG) {
//        fht.remove_key("BZERO");
//        fht.remove_key("BSCALE");
//    }
}

void PixelArray2dHeader::load_from_fht(const FitsHeader& fht)
{
    if (!fht.key_exists("SIMPLE")) {
        throw std::runtime_error("fits_header::load_essential_FITS_keys "
                                 "- table Does not appear to have a valid header (No SIMPLE keyword): ");
    }
    bitpix_=(fits_header_conventions::fits_imagetype) string_utils::atoi(
                fht.get_key_value("BITPIX"));
    int naxes= atoi(fht.get_key_value("NAXIS"));
    if (naxes!=2) { throw std::runtime_error("File is not 2D"); }
    else {
        x_dim_= atoi(fht.get_key_value("NAXIS1"));
        y_dim_= atoi(fht.get_key_value("NAXIS2"));
    }
    n_pix=x_dim_*y_dim_;
}

size_t PixelArray2dHeader::byte_size_of_uncompressed_data() const
{
    return n_pix * fits_header_conventions::byte_size_of_fits_imagetype_pixel(bitpix());
}

bool PixelArray2dHeader::operator==(const PixelArray2dHeader& rhs) const
{
    return //n_axes()==rhs.n_axes() &&
        n_pixels()==rhs.n_pixels() &&
        x_dim()==rhs.x_dim() &&
        y_dim()==rhs.y_dim() &&
        bitpix() ==rhs.bitpix();
}

//=========================================================================================================
//PixelCubeHeader

void PixelCubeHeader::load_from_fht(const FitsHeader& fht)
{
    if (!fht.key_exists("SIMPLE")) {
        throw std::runtime_error("PixelCubeHeader::load_from_fht "
                                 "- table Does not appear to have a valid header (No SIMPLE keyword): ");
    }
    bitpix_=(fits_header_conventions::fits_imagetype) string_utils::atoi(
                fht.get_key_value("BITPIX"));
    int naxes= atoi(fht.get_key_value("NAXIS"));
    if (naxes!=3) { throw std::runtime_error("File is not a pixel cube"); }
    else {
        x_dim_= atoi(fht.get_key_value("NAXIS1"));
        y_dim_= atoi(fht.get_key_value("NAXIS2"));
        z_dim_= atoi(fht.get_key_value("NAXIS3"));
    }
    return;
}

//=========================================================================================================
}//end namespace coela
