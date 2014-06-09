/*
 * File:   array_compression_info.h
 * Author: ts337
 *
 * Created on 10 February 2011, 15:44
 */

#ifndef COELA_ARRAY_COMPRESSION_INFO_H
#define COELA_ARRAY_COMPRESSION_INFO_H

#include "FitsHeader.h"

namespace coela {
//=============================================================================

/**

 Struct to encode all information stored in a fits header about how the image
 data is compressed.
 NB currently this is passed around *by value* so that we can use a default
 "no_compression" option without creating a dangling pointer
 (Which would happen if we passed it by reference, but then tried to pass a
 temporarily constructed "no_compression" object).

 */
//THEREFORE MEMORY FOOTPRINT SHOULD BE KEPT MINIMAL
struct ArrayCompressionInfo {
    ///No compression is the default.
//    ArrayCompressionInfo(){*this = ArrayCompressionInfo::no_compression();}

    static ArrayCompressionInfo no_compression();
    ArrayCompressionInfo(const FitsHeader& fht) {load_from_fht(fht);}

    void write_to_fht(FitsHeader& fht) const;
    void load_from_fht(const FitsHeader& fht);



    bool is_lcc_compressed;

    ///0 for NLaw byte compression, 1 for FSuess bit compression.
    int lcc_compression_type;

    int lc_compressed_pixel_bits, lc_uncompressed_pixel_bits, lc_level_offset;

private:
    //prevent creation of uninitialized instances.
    ArrayCompressionInfo() {}

};
//======================================================================================
}//end namespace
#endif  /* ARRAY_COMPRESSION_INFO_H */

