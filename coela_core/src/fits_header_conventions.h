/*
 * File:   fits_header_conventions.h
 * Author: ts337
 *
 * Created on 09 February 2011, 12:29
 */

#ifndef COELA_FITS_HEADER_CONVENTIONS_H
#define COELA_FITS_HEADER_CONVENTIONS_H
#include <string>

namespace coela {
//===================================================================================
///Fits conventions

namespace fits_header_conventions {
const unsigned int key_max=8, line_max=80,
                   rows_per_card=36, value_col_max=30,
                   value_max=20, kv_comment_max = 46 , bytes_per_card=2880;

const std::string standard_first_ten_chars ="SIMPLE  = ";
const std::string comment_prefix="COMMENT ";

enum fits_imagetype { FLOATIMG=-32, USHORTIMG=20, DOUBLEIMG=-64, BYTEIMG=8,
                      SHORTIMG=16, LONGIMG=32
                    };

size_t byte_size_of_fits_imagetype_pixel(const fits_imagetype T);


};

//=======================================================================================


//=======================================================================================
//FITS conventions:
//const std::string fits_header_conventions::comment_prefix="COMMENT ";

//std::ostream& operator<<(ostream& os, const fits_imagetype& i) //Displays enumerations rather than integer value of 'fits_imagetype'
//{
//    switch(i){
//        case FLOATIMG:
//            os << "FLOATIMG";
//            break;
//
//        case USHORTIMG:
//             os << "USHORTIMG";
//            break;
//
//       case DOUBLEIMG:
//             os << "DOUBLEIMG";
//            break;
//
//        case BYTEIMG:
//             os << "BYTEIMG";
//            break;
//
//        case SHORTIMG:
//             os << "SHORTIMG";
//            break;
//
//       case LONGIMG:
//             os << "LONGIMG";
//            break;
//    }
//
//     return os;
//}

//=======================================================================================

}//end namespace coela

#endif  /* FITS_HEADER_CONVENTIONS_H */

