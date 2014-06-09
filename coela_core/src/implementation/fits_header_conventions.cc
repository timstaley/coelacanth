#include "../fits_header_conventions.h"

#include <stdexcept>


namespace coela {
namespace fits_header_conventions {
size_t byte_size_of_fits_imagetype_pixel(const fits_imagetype T)
{
    switch (T) {
    case FLOATIMG:
        return 4;
    case DOUBLEIMG:
        return 8;
    case SHORTIMG:
        return 2;
    default:
        throw std::runtime_error("byte_size_of_fits_imagetype_pixel: uknown imagetype");
    }
}



}//end namespace coela::fits_header_conventions
}//end namespace coela