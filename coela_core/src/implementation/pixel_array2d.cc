/*
 * File:   pixel_array.cc
 * Author: ts337
 *
 * Created on 09 February 2011, 10:32
 */

#include "../pixel_array2d.h"
#include "coela_utility/src/string_utils.h"
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <cstring>
#include <stdint.h>
using std::runtime_error;
namespace coela {

//=========================================================================================================
//Pixel iterator
PixelIterator::PixelIterator(const PixelRange& rgn):PixelIndex(rgn.low), end(rgn.low.x,
            rgn.high.y+1), range(rgn)
{
    //NB end must be what we get when we increment 1 past rgn.high - see increment operator below:
}

//prefix increment
PixelIterator& PixelIterator::operator++()
{
    if (x!=range.high.x) {
        ++x;
    } else {
        x=range.low.x;
        ++y;
    }
    return *this;
}

//suffix increment
PixelIndex PixelIterator::operator++(int)
{
    PixelIndex pre_val(*this);
    ++(*this); //call the prefix oper;
    return pre_val;
}

//=========================================================================================================


//=========================================================================================================
//Pixel array
//--------------------------------------------------------------------------------------------
//These must be place before other partial specializations, apparently.
//Returns the fits bitpix flag for each data type:
template<>
fits_header_conventions::fits_imagetype PixelArray2d<float>::array_FITS_imagetype()
{
    return fits_header_conventions::FLOATIMG;
}

template<>
fits_header_conventions::fits_imagetype PixelArray2d<double>::array_FITS_imagetype()
{
    return fits_header_conventions::DOUBLEIMG;
}

template<>
fits_header_conventions::fits_imagetype PixelArray2d<uint32_t>::array_FITS_imagetype()
{
    return fits_header_conventions::LONGIMG;
}
template<>
fits_header_conventions::fits_imagetype PixelArray2d<int32_t>::array_FITS_imagetype()
{
    return fits_header_conventions::LONGIMG;
}

template<>
fits_header_conventions::fits_imagetype PixelArray2d<uint16_t>::array_FITS_imagetype()
{
    return fits_header_conventions::SHORTIMG;
}
template<>
fits_header_conventions::fits_imagetype PixelArray2d<int16_t>::array_FITS_imagetype()
{
    return fits_header_conventions::SHORTIMG;
}

//--------------------------------------------------------------------------------------------
//constructors:
template<typename T>
PixelArray2d<T>::PixelArray2d(const size_t x_dim_init, const size_t y_dim_init,
                              T init_value):
    hdr(PixelArray2dHeader(PixelArray2d<T>::array_FITS_imagetype(), x_dim_init, y_dim_init))
{
    resize_data_vec_according_to_header_requirements();
    assign(init_value);
}


//template<typename T>
//template<typename U>
//PixelArray2d<T>::PixelArray2d(const PixelArray2d<U>&)
//{
//    throw std::logic_error("PixelArray2d<T>::PixelArray2d(const PixelArray2d<U>& rhs):\n"
//                            "Promotion from any type to non floating point types is not implemented.");
//}

template< typename T >
template < typename U >
PixelArray2d<T>::PixelArray2d(const PixelArray2d<U>& rhs):
    hdr(PixelArray2d<T>::array_FITS_imagetype() , rhs.range().x_dim(),
        rhs.range().y_dim())
{
    resize_data_vec_according_to_header_requirements();
    std::copy(rhs.data_.begin(), rhs.data_.end(), data_.begin());
}

//template< >
//template <typename U>
//PixelArray2d<float>::PixelArray2d(const PixelArray2d<U>& rhs):
//hdr(PixelArray2d<float>::array_FITS_imagetype() , rhs.range().x_dim(),
//        rhs.range().y_dim()){
//    resize_data_vec_according_to_header_requirements();
//    std::copy(rhs.data_.begin(), rhs.data_.end(), data_.begin());
//}

template PixelArray2d<double>::PixelArray2d(const PixelArray2d<int >& rhs);
template PixelArray2d<double>::PixelArray2d(const PixelArray2d<unsigned >& rhs);
template PixelArray2d<double>::PixelArray2d(const PixelArray2d<unsigned short>& rhs);
template PixelArray2d<double>::PixelArray2d(const PixelArray2d<float>& rhs);

template PixelArray2d<float>::PixelArray2d(const PixelArray2d<int >& rhs);
template PixelArray2d<float>::PixelArray2d(const PixelArray2d<double>& rhs);

template<typename T>
PixelArray2d<T>::PixelArray2d(const PixelArray2dHeader& hdr_init):
    hdr(hdr_init)
{
    resize_data_vec_according_to_header_requirements();
}



template<typename T>
PixelArray2d<T> PixelArray2d<T>::sub_array(const PixelArray2d& orig,
        const PixelRange& sub_box)
{
    PixelArray2d<T> subarray(PixelArray2dHeader::subarray_header(orig.hdr, sub_box));

    //Regular version
//    for (PixelIterator pixin(sub_box), pixout(subarray.outline()); pixin!=pixin.end; ++pixin, ++pixout){
//        subarray(pixout)=orig(pixin);
//    }

    //Optimization (un-benchmarked) (TO DO - benchmark!)
    //copy row at a time (and hope std::copy is well optimized )
    size_t orig_width = orig.range().x_dim();
    size_t subarray_width = sub_box.x_dim();

    for (int y_in=sub_box.low.y, y_out=1; y_in<=sub_box.high.y; ++y_in, ++y_out) {
        //increment the row
        int in_row_start = (y_in-1)*orig_width+sub_box.low.x-1;
        int in_row_end = (y_in-1)*orig_width+sub_box.high.x;
        int out_row_start = (y_out-1)*subarray_width;
        std::copy(orig.data_.begin() + in_row_start, orig.data_.begin() + in_row_end,
                  subarray.data_.begin()+ out_row_start);
    }

    return subarray;
}

template<typename T>
PixelArray2d<T> PixelArray2d<T>::load_from_file(const string& filename,
        const size_t hdr_begin_byte_offset)
{
    std::ifstream infile(filename.c_str(), std::ios::binary);
    if (!infile.is_open()) { throw std::runtime_error(filename+" won't open"); }
    PixelArray2d<T> pix;
    FitsHeader fht(infile, hdr_begin_byte_offset);
    pix.hdr = PixelArray2dHeader(fht);
    pix.seek_and_load_data_from_stream(infile,
                                       hdr_begin_byte_offset+fht.header_file_length_in_bytes(),
                                       ArrayCompressionInfo(fht));
    infile.close();
    return pix;
}


template<class T>
PixelArray2d<T> PixelArray2d<T>::load_from_buffer(const FileBuffer& buf,
        const size_t hdr_begin_byte_offset)
{
    FitsHeader fht(buf, hdr_begin_byte_offset);
    PixelArray2d<T> pix;
    pix.hdr = PixelArray2dHeader(fht);
    pix.load_data_from_buffer(buf,
                              hdr_begin_byte_offset + fht.header_file_length_in_bytes(),
                              ArrayCompressionInfo(fht));
    return pix;
}


template<typename T>
void PixelArray2d<T>::write_to_file(const string& filename,
                                    const FitsHeader additional_header_info,
                                    const ArrayCompressionInfo compression_inf) const
{
    std::ofstream outfile(filename.c_str(), std::ios::binary);
    if (!outfile.is_open()) throw std::runtime_error(
            "pixel_array<T>::write_to_file - error writing to file "+ filename);

    FitsHeader fht(additional_header_info);
    PixelArray2d<T>::write_header_to_fht(fht, compression_inf);
    outfile<<fht;

    PixelArray2d<T>::write_data_to_stream(outfile, compression_inf);
    outfile.close();
    return;
}


template<typename T>
std::istream& PixelArray2d<T>::seek_and_load_data_from_stream(std::istream& is,
        const size_t data_begin_byte_offset,
        const ArrayCompressionInfo& compression_inf)
{
    assert(hdr.n_pixels()); //Check header initialized

    //Doesn't matter what type if loading from compressed file, as it will
    //just get cast to whatever data type as we process it.
    //However, for a regular fits the pixel byte lengths must be correct
    if (compression_inf.is_lcc_compressed==false) {
        check_header_type_matches_array_type();
    }


    is.seekg(0, std::ios::end);
    //NB seek stream position set to end - this allows us to check that
    //the file is long enough that it probably contains the data the header suggests it will.

    size_t expected_minimum_data_length =
        calculate_minimum_data_size_in_bytes(hdr.n_pixels(), compression_inf);

    if (is.tellg() < std::streampos(
                data_begin_byte_offset + expected_minimum_data_length)) {
        throw runtime_error("PixelArray2d<T>::load_from_stream - "
                            "stream is too short, presumed data missing");
    }
    is.seekg(data_begin_byte_offset, std::ios::beg);

    return load_data_from_positioned_stream(is, compression_inf);
}
//
//    resize_data_vec_according_to_header_requirements();
//
//    if (compression_inf.is_lcc_compressed==false){
//        check_header_type_matches_array_type();
//        load_array_from_regular_FITS_stream(is);
//    } //else, compression must be on...
//    else if (compression_inf.lcc_compression_type==1){
//        load_array_from_lcz_bits_stream(is, compression_inf);
//    }
//    else if (compression_inf.lcc_compression_type==0){
//        load_array_from_lcc_byte_stream(is);
//    }
//    else throw runtime_error("pixel_array<T>::load_array_from_stream - Un-implemented compression type for array input");
//    return is;
//}

template<typename T>
std::istream& PixelArray2d<T>::load_data_from_positioned_stream(std::istream& is,
        const ArrayCompressionInfo& compression_inf)
{
    resize_data_vec_according_to_header_requirements();

    if (compression_inf.is_lcc_compressed==false) {
        check_header_type_matches_array_type();
        load_array_from_regular_FITS_stream(is);
    } //else, compression must be on...
    else if (compression_inf.lcc_compression_type==1) {
        load_array_from_lcz_bits_stream(is, compression_inf);
    } else if (compression_inf.lcc_compression_type==0) {
        load_array_from_lcc_byte_stream(is);
    } else throw runtime_error("pixel_array<T>::load_array_from_stream - "
                                   "Un-implemented compression type for array input");
    return is;

}


template<typename T>
void PixelArray2d<T>::load_data_from_buffer(const FileBuffer& buf,
        const size_t data_begin_offset,
        const ArrayCompressionInfo& compression_inf)
{
    size_t expected_minimum_data_length= calculate_minimum_data_size_in_bytes(hdr.n_pixels(),
                                         compression_inf);

    if (buf.data_vec().size() < data_begin_offset + expected_minimum_data_length) {
        throw runtime_error("PixelArray2d<T>::load_from_buffer - "
                            "buf is too short, presumed data missing");
    }

    resize_data_vec_according_to_header_requirements();

    if (compression_inf.is_lcc_compressed==false) {
        check_header_type_matches_array_type();
        load_array_from_regular_FITS_buffer(buf, data_begin_offset);
    } //else, compression must be on...
    else if (compression_inf.lcc_compression_type==1) {
        load_array_from_lcz_bits_buffer(&buf.data_vec()[0]+ data_begin_offset    ,
                                        &buf.data_vec()[ buf.data_vec().size() ]  ,
                                        compression_inf);
    } else if (compression_inf.lcc_compression_type==0) {
        load_array_from_lcc_byte_buffer(&buf.data_vec()[0]+ data_begin_offset    ,
                                        &buf.data_vec()[ buf.data_vec().size() ]);
    } else throw runtime_error("pixel_array<T>::load_data_from_buffer - "
                                   "Un-implemented compression type for array input");
}

template<typename T>
size_t PixelArray2d<T>::calculate_minimum_data_size_in_bytes(size_t n_pixels,
        const ArrayCompressionInfo& compression_inf)
{
    if (compression_inf.is_lcc_compressed==false) {
        //regular FITS case
        return n_pixels*sizeof(T);
    } else if (compression_inf.lcc_compression_type==1) {
        //compressed case, bitstream
        return n_pixels*compression_inf.lc_compressed_pixel_bits/8;
    } else if (compression_inf.lcc_compression_type==0) {
        //compressed case, lcc byte stream
        return n_pixels; //min 1 byte per pixel
    } else { throw runtime_error("PixelArray2d<T>::load_from_stream - unknown compression type"); }
}


template<typename T>
void PixelArray2d<T>::write_header_to_fht(FitsHeader& fht,
        const ArrayCompressionInfo compression_inf) const
{
    hdr.write_to_fht(fht, array_FITS_imagetype());
    compression_inf.write_to_fht(fht);
}





template<typename T>
std::ostream& PixelArray2d<T>::write_data_to_stream(std::ostream& os,
        const ArrayCompressionInfo& compression_inf) const
{
    if (compression_inf.is_lcc_compressed==false) {
        write_array_to_regular_FITS_stream(os);
    } else throw runtime_error("pixel_array<T>::write_array_to_stream -"
                                   " Un-implemented compression type for array output");
    return os;
}

template<typename T>
void PixelArray2d<T>::assign(T value)
{
    data_.assign(hdr.n_pixels(), value);
}

template<typename T>
bool PixelArray2d<T>::operator==(const PixelArray2d& rhs) const
{
    return hdr==rhs.hdr &&
           std::equal(data_.begin(), data_.end(), rhs.data_.begin());
}

template<typename T>
PixelArray2d<double> PixelArray2d<T>::operator-(const double subtractor)
{
    PixelArray2d<double> temp(*this);
    temp -=subtractor;
    return temp;
}

template<typename T>
PixelIndex PixelArray2d<T>::min_PixelIndex() const
{
    assert(!data_.empty());
    T min = data_.front();
    size_t min_index=0;
    for (size_t i=0; i!=data_.size(); ++i) {
        if (data_[i]<min) {
            min = data_[i];
            min_index = i;
        }
    }
    return PixelIndex_of_data_element(min_index);
}

template<typename T>
PixelIndex PixelArray2d<T>::max_PixelIndex() const
{
    assert(!data_.empty());
    T max = data_.front();
    size_t max_index=0;
    for (size_t i=0; i!=data_.size(); ++i) {
        if (data_[i] > max) {
            max = data_[i];
            max_index = i;
        }
    }
    return PixelIndex_of_data_element(max_index);
}

template<typename T>
PixelIndex PixelArray2d<T>::max_PixelIndex(const PixelRange& rgn_box) const
{
    assert(range().contains_range(rgn_box));
    PixelIndex max_index = rgn_box.low;
    T max = (*this)(rgn_box.low);

    for (PixelIterator it(rgn_box); it!=it.end; ++it) {
        if ((*this)(it) > max) {
            max = (*this)(it);
            max_index = it;
        }
    }
    return max_index;
}

template<typename T>
T PixelArray2d<T>::region_sum(const PixelRange& rgn) const
{
    T sum(0);
    for (PixelIterator i(rgn); i!=i.end; ++i) { sum+=(*this)(i); }
    return sum;
}

template<typename T>
void PixelArray2d<T>::operator*=(const PixelArray2d& mask)
{
    assert(mask.range()==range());
    for (size_t i=0; i!=data_.size(); ++i) {
        data_[i]*=mask.data_[i];
    }
}

template<typename T>
void PixelArray2d<T>::operator/=(const PixelArray2d& divisor)
{
    assert(divisor.range()==range());
    for (size_t i=0; i!=data_.size(); ++i) {
        data_[i]/=divisor.data_[i];
    }
}


template<typename T>
void PixelArray2d<T>::operator-=(const PixelArray2d& subtraction_frame)
{

    assert(subtraction_frame.range()==range());
    for (size_t i=0; i!=data_.size(); ++i) {
        data_[i]-=subtraction_frame.data_[i];
    }
}

template<typename T>
void PixelArray2d<T>::operator+=(const PixelArray2d& summation_frame)
{
    assert(summation_frame.range()==range());
    for (size_t i=0; i!=data_.size(); ++i) {
        data_[i]+=summation_frame.data_[i];
    }
}

template<typename T>
template<typename U>
void PixelArray2d<T>::operator+=(const PixelArray2d<U>& summation_frame)
{
    assert(summation_frame.range()==range());
    for (size_t i=0; i!=data_.size(); ++i) {
        data_[i]+=summation_frame.data_[i];
    }
}

template void PixelArray2d<float>::operator +=(const PixelArray2d<int >& rhs);
template void PixelArray2d<float>::operator +=(const PixelArray2d<unsigned short >& rhs);
template void PixelArray2d<float>::operator +=(const PixelArray2d<unsigned  >& rhs);

template void PixelArray2d<double>::operator +=(const PixelArray2d<float >& rhs);
template void PixelArray2d<double>::operator +=(const PixelArray2d<unsigned short >& rhs);
template void PixelArray2d<double>::operator +=(const PixelArray2d<unsigned  >& rhs);
template void PixelArray2d<double>::operator +=(const PixelArray2d<int >& rhs);


template<typename T>
void PixelArray2d<T>::resize_data_vec_according_to_header_requirements()
{
    if (data_.size()!= hdr.n_pixels()) {
        data_.resize(hdr.n_pixels());
    }
}
template<typename T>
PixelIndex PixelArray2d<T>::PixelIndex_of_data_element(size_t i) const
{
    return PixelIndex(i% hdr.x_dim() +1,  i / hdr.x_dim() +1);
}

template<typename T>
void PixelArray2d<T>::check_header_type_matches_array_type() const
{
    if (array_FITS_imagetype() != hdr.bitpix())
        throw std::runtime_error("PixelArray2d<T>::check_header_type_matches_array_type() \n"
                                 "array / header data type mismatch, header: " + string_utils::itoa(
                                     hdr.bitpix()) +" array: " + string_utils::itoa(array_FITS_imagetype()));
}

template<typename T>
std::istream& PixelArray2d<T>::load_array_from_regular_FITS_stream(std::istream& is)
{
    is.read((char*)(&data_[0]),hdr.n_pixels() * sizeof(T));
    //Swap byte (and pixel) ordering.
    std::reverse((char*)&data_[0], (char*)(&data_[0] + data_.size()));
    std::reverse(data_.begin(), data_.end());    //swap pixel ordering back.
//May be quicker to byte swap each individual pixel, just once - but not sure how best to code that. This works well enough.
    return is;
}

template<typename T>
void PixelArray2d<T>::load_array_from_regular_FITS_buffer(const FileBuffer& buf,
        const size_t data_begin_offset)
{
    std::reverse_copy(buf.data_vec().begin()+data_begin_offset,
                      buf.data_vec().begin() +data_begin_offset + hdr.n_pixels()*sizeof(T),
                      (char*)&data_[0]
                     );//Copy, swapping byte (and pixel) ordering.

    std::reverse(data_.begin(), data_.end());    //swap pixel ordering back.
//May be quicker to byte swap each individual pixel, just once - but not sure how best to code that. This works well enough.
}

template<typename T>
std::ostream& PixelArray2d<T>::write_array_to_regular_FITS_stream(std::ostream& os) const
{
    T* swap_buffer;
    try {        swap_buffer = new T[hdr.n_pixels()];    }
    catch (std::bad_alloc b) {
        throw std::runtime_error("pixel_array<T>::write_to_regular_stream - Could not allocate enough memory for output swap buffer");
    }

//    memcpy(swap_buffer, data_begin,header.n_elements()*sizeof(T)); //copy array
//    std::reverse((char*)swap_buffer, (char*)swap_buffer+header.n_elements()*sizeof(T));  //Reverse byte ordering of array

    //copy array in reversed byte order, into the swap buffer;
    std::reverse_copy((char*)&data_[0], (char*)(&data_[0] + data_.size()),
                      (char*)swap_buffer);


    std::reverse(swap_buffer, swap_buffer+hdr.n_pixels());   //Switch back pixel ordering.
    os.write((char*)swap_buffer,hdr.n_pixels()*sizeof(T));
    delete[] swap_buffer;
    return os;
}

union sixty_four_bits_union {
    struct        {
        uint32_t LowPart;
        uint32_t HighPart;
    };
    uint64_t AllBits;
};


template<typename T>
void PixelArray2d<T>::load_array_from_lcz_bits_buffer(
    const char * data_begin_ptr, const char * buffer_end_ptr,
    const ArrayCompressionInfo& compression_info
)
{
    using string_utils::atoi;
    int CompressedBits   = compression_info.lc_compressed_pixel_bits ;
    int UncompressedBits = compression_info.lc_uncompressed_pixel_bits;
//        unsigned short LevelOffset = (unsigned short)atoi(get_key_value("LC_OFFST"));
    float LevelOffset = compression_info.lc_level_offset;
    typedef uint32_t DWORD; //perhaps more readable, could just find and replace though.

    //NB CompressedBits==0 means all pixels are stored as 14 bits.

    //All 1s, bit shifted right until only the compressed bits are 1 and the rest is zero
    DWORD CompressedPixelMask   = (DWORD) 0xffffffff >> (32 - CompressedBits);
    // similar but with more 1s.
    DWORD UncompressedPixelMask = (DWORD) 0xffffffff >> (32 - UncompressedBits);


    //destination - "float * data" has already been allocated for n_elements

    //set up bitwise operation pointers:
    DWORD *pSrc = (DWORD*)(data_begin_ptr);
    DWORD *pEnd = (DWORD*)(buffer_end_ptr);
    int WordCount = 0;
    sixty_four_bits_union shift_register;

    shift_register.LowPart = pSrc[WordCount++];
    shift_register.HighPart =pSrc[WordCount++];

    //tracks how many bits the current 64 bits in shift_register have been shifted.
    DWORD BitShift = 0;
    //When this gets to 32 we must load in another DWORD

    for (size_t i = 0; i <hdr.n_pixels() && pSrc< pEnd; i++) {
        DWORD BitCount;  //Counts how many bits we have just read and are about to shift

        if (CompressedBits) {
            BitCount = CompressedBits;
            // Fetch compressed datum
            DWORD Pixel = shift_register.LowPart & CompressedPixelMask;

            if (Pixel)  {
                // compressed pixel value (6 zero bits signifies uncompressed)
//                    data[i] = LevelOffset + (unsigned short) Pixel;  // compressed pixel value
                data_[i] = LevelOffset + Pixel;  // compressed pixel value
            } else { // uncompressed flag followed by uncompressed pixel value
                //About to shift away the flag bits, so check we have room...
                // If too many bits to shift in one go (register would overflow), split register shift into two
                if (BitShift + BitCount > 32) {
                    //we are moving the last few bits of the high part into the low part
                    BitShift = 32 - BitShift;

                    shift_register.AllBits >>= BitShift;
//                        ShiftRegister.QuadPart = Int64ShrlMod32( ShiftRegister.QuadPart, BitShift );  // first part of shift

                    BitCount -= BitShift; //so bitcount can be reduced by amount just shifted
                    BitShift  = 0; //and boundaries are now aligned

                    shift_register.HighPart = pSrc[ WordCount++ ];  // fetch next 32-bit word
                }
                shift_register.AllBits >>= BitCount;  //Complete the 6bit shift.
//                    ShiftRegister.QuadPart = Int64ShrlMod32( ShiftRegister.QuadPart, BitCount );  // shift past the 6-bit flag

                BitShift += BitCount; //Update the shift tracker
                BitCount  = UncompressedBits;

//                    data[ i ] = (unsigned short) ( shift_register.LowPart & UncompressedPixelMask );  // uncompressed pixel value
                // uncompressed pixel value
                data_[ i ] = (shift_register.LowPart & UncompressedPixelMask);
            }
        } else {
            BitCount = UncompressedBits;  // all pixels are uncompressed (well, 14 bits)
//                    data[ i ] = (unsigned short) ( shift_register.LowPart & UncompressedPixelMask );  // uncompressed pixel value
            // uncompressed pixel value
            data_[ i ] = (shift_register.LowPart & UncompressedPixelMask);
        }
        //Now shift again
        // If too many bits to shift in one go (register would overflow), split register shift into two
        if (BitShift + BitCount > 32) {
            //about to shift away what we have just read
            BitShift = 32 - BitShift;
            shift_register.AllBits>>= BitShift;
//                    ShiftRegister.QuadPart = Int64ShrlMod32( ShiftRegister.QuadPart, BitShift );  // first part of shift
            BitCount -= BitShift;
            BitShift  = 0;

            shift_register.HighPart = pSrc[ WordCount++ ];  // fetch next 32-bit word
        }

        shift_register.AllBits>>= BitCount;
//          ShiftRegister.QuadPart = Int64ShrlMod32( ShiftRegister.QuadPart, BitCount );  // complete the shift
        BitShift += BitCount;
    }
}

template<typename T>
std::istream& PixelArray2d<T>::load_array_from_lcz_bits_stream(
    std::istream& is, const ArrayCompressionInfo& compression_info)
{
    //read 24 bits per pixel- input will always be 20bits per pixel or less
    char* input_buffer=new char[hdr.n_pixels()*3];
    //so if we read in this much, we're covered

    //this reads in 24 bits * n_pixels - more than enough
    is.read((char*)input_buffer,hdr.n_pixels() * 3);

    load_array_from_lcz_bits_buffer(input_buffer, input_buffer+hdr.n_pixels()*3,
                                    compression_info);

    delete[] input_buffer;
    return is;
}

template<typename T>
void PixelArray2d<T>::load_array_from_lcc_byte_buffer(
    const char * data_begin_ptr, const char * buffer_end_ptr)
{
    signed char * buffer =(signed char*)(void*)data_begin_ptr;
    int16_t *two_bytes;
    const size_t n_elts = hdr.n_pixels();
    size_t out(0), n_bytes(buffer_end_ptr-data_begin_ptr);
    for (size_t in(0); in<n_bytes && out < n_elts; ++in, ++out) {
        if (buffer[in]!= -128) { data_[out]= buffer[in]; }
        else {
            two_bytes = (int16_t*)&buffer[in+1]; //move to next byte, since current char is a prefix.
            data_[out] = *two_bytes;
            in+=2;
        }
    }
    if (out!=n_elts) {
        throw runtime_error("PixelArray2d<T>::load_array_from_lcc_byte_buffer - Could not load enough pixels from file");
    }
}
template<typename T>
std::istream& PixelArray2d<T>::load_array_from_lcc_byte_stream(std::istream& is)
{
    using std::ptrdiff_t;
    const size_t n_elts = hdr.n_pixels();
    //If we wanted to be extra safe we could initialise this but it will slow up.
    char* buffer=new char[n_elts*3];
    ptrdiff_t stream_start_pos = is.tellg();
    is.read((char*)buffer, n_elts * 3); //need the cast to char for is.read
    ptrdiff_t stream_end_pos = is.tellg();
    //will either read n*3, or less if we hit EoF
    size_t bytes_read = stream_end_pos - stream_start_pos;
    load_array_from_lcc_byte_buffer(buffer, buffer+bytes_read);
    delete[] buffer;
    return is;
}
//template<typename T>
//std::ostream& pixel_array<T>::write_array_to_lcz_bits_stream(std::ostream& os) const{
//    throw logic_error()
//}




//Supported types:
//NB! THIS MUST GO AT THE END OF THIS FILE!
//Otherwise the compiler tries to instantiate the templatized classes before it has all the relevant information. D'oh!
template class PixelArray2d<float>;
template class PixelArray2d<double>;

template class PixelArray2d<uint16_t>;
template class PixelArray2d<uint32_t>;
template class PixelArray2d<int16_t>;
template class PixelArray2d<int32_t>;



//=========================================================================================================
}//end namespace coela
