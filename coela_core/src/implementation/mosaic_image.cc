#include "../mosaic_image.h"
#include <stdint.h>
namespace coela {

const double error_margin = 1e-6;

template<typename T>
MosaicImage<T>::MosaicImage() {}

template<typename T>
MosaicImage<T>::MosaicImage(const CcdImage<T>& rhs):
    CcdImage<T>(rhs)
{}

template<typename T>
void MosaicImage<T>::initialize_mosaic_grid_to_specific_region(
    const MosaicBoxRegion& mosaic_region, bool no_throw)
{
    if (mosaic_grid.is_initialized() && !no_throw)
        throw std::runtime_error("Image<T>::initialize_mosaic_grid_to_specific_layout() -"
                                 " mosaic_grid already initialized");

    double deduced_pixel_x_scale = mosaic_region.x_dim() /
                                   (double)this->pix.range().x_dim();
    assert(
        fabs(deduced_pixel_x_scale  - mosaic_region.y_dim() /
             (double)this->pix.range().y_dim())
        < error_margin
    ); //assert (square pixels)
    mosaic_grid = ImageGrid<coordinate_types::mosaic>(mosaic_region,
                  deduced_pixel_x_scale);
}

template<typename T>
void MosaicImage<T>::write_to_file(const std::string& filename,
                                   FitsHeader fht,
                                   const ArrayCompressionInfo aci) const
{
    mosaic_grid.write_to_fht(fht);
    CcdImage<T>::write_to_file(filename, fht, aci);
}


//Supported types:
//NB THIS MUST GO AT THE END OF THE FILE
//Otherwise the compiler tries to instantiate the templatized classes before
//it has all the relevant information.
template class MosaicImage<float>;
template class MosaicImage<double>;

template class MosaicImage<int>;
template class MosaicImage<uint16_t>;
template class MosaicImage<uint32_t>;



} //end namespace coela
