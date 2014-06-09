/*
 * File:   array_compression_info.cc
 * Author: ts337
 *
 * Created on 10 February 2011, 15:44
 */

#include "../array_compression_info.h"
#include "coela_utility/src/string_utils.h"

namespace coela {
using namespace string_utils;
//======================================================================================

ArrayCompressionInfo ArrayCompressionInfo::no_compression()
{
    ArrayCompressionInfo inf;
    inf.is_lcc_compressed=false;

    inf.lcc_compression_type=-1;
    inf.lc_compressed_pixel_bits=-1;
    inf.lc_uncompressed_pixel_bits=-1;
    inf.lc_level_offset=-1;

    return inf;
}

void ArrayCompressionInfo::write_to_fht(FitsHeader& fht) const
{
    if (is_lcc_compressed) {
        fht.set_keyword("LC_FORMT", string_utils::itoa(lcc_compression_type));

        if (lcc_compression_type==0) {
            fht.set_keyword("LCL_ZERO", string_utils::itoa(lc_level_offset));
        }
        if (lcc_compression_type==1) {
            fht.set_keyword("LC_OFFST", string_utils::itoa(lc_level_offset));
            fht.set_keyword("LC_CBITS", string_utils::itoa(lc_compressed_pixel_bits));
            fht.set_keyword("LC_UBITS", string_utils::itoa(lc_uncompressed_pixel_bits));
        }

    } else {
        fht.remove_key("LCL_ZERO");
        fht.remove_key("LC_FORMT");
        fht.remove_key("LC_OFFST");
    }

}
void ArrayCompressionInfo::load_from_fht(const FitsHeader& fht)
{

    if (fht.key_exists("LC_FORMT")) {
        is_lcc_compressed=true;
        lcc_compression_type=atof(fht.get_key_value("LC_FORMT"));
        if (lcc_compression_type==1) {
            lc_level_offset = atoi(fht.get_key_value("LC_OFFST"));
            lc_compressed_pixel_bits = atoi(fht.get_key_value("LC_CBITS"));
            lc_uncompressed_pixel_bits = atoi(fht.get_key_value("LC_UBITS"));
        } else if (lcc_compression_type==0) {
            lc_level_offset= atoi(fht.get_key_value("LCL_ZERO"));
        }

    } else if (fht.key_exists("LCL_ZERO")==true) {
        //legacy for NLaw lcc files before introduction of LC_FORMT keyword
        is_lcc_compressed=true;
        lcc_compression_type=0;
        lc_level_offset=atoi(fht.get_key_value("LCL_ZERO"));
    } else { is_lcc_compressed=false; }
}

//======================================================================================
}//end namespace