
#include <fstream>
#include <stdint.h>
#include "../ccd_image.h"
namespace coela {

const double error_margin = 1e-6;

//=========================================================================================================

template<typename T>
CcdImage<T>::CcdImage() {}

template<typename T>
CcdImage<T>::CcdImage(const std::string& filename, const size_t header_offset)
{
    std::ifstream infile(filename.c_str(), std::ios::binary);
    if (!infile.is_open()) { throw std::runtime_error(filename+" won't open"); }
    FitsHeader fht(infile, header_offset);
    load_from_stream(fht, infile, header_offset);
    infile.close();
}

template<typename T>
template<typename U>
CcdImage<T>::CcdImage(const CcdImage<U>& rhs):
    pix(rhs.pix),
    CCD_grid(rhs.CCD_grid) {}

//template<>
//template<typename U>
//Image<float>::Image(const Image<U>& rhs):
//        pix(rhs.pix),
//        CCD_gridinitialized_(rhs.CCD_gridinitialized_),
//        CCD_grid(rhs.CCD_grid),
//        mosaic_grid_initialized_(rhs.mosaic_grid_initialized_),
//        mosaic_grid_(rhs.mosaic_grid_)
//{}

template CcdImage<double>::CcdImage(const CcdImage<int >& rhs);
template CcdImage<double>::CcdImage(const CcdImage<unsigned >& rhs);
template CcdImage<double>::CcdImage(const CcdImage<float>& rhs);
template CcdImage<float>::CcdImage(const CcdImage<double>& rhs);
template CcdImage<float>::CcdImage(const CcdImage<int >& rhs);

template <>
CcdImage<double> CcdImage<double>::load_from_unknown_filetype(
    const std::string& filename,
    const size_t header_begin_offset_in_bytes)
{
    FitsHeader fht(filename, header_begin_offset_in_bytes);
    if (fht.FITS_imagetype()==fits_header_conventions::DOUBLEIMG) {
        return CcdImage<double>(filename, header_begin_offset_in_bytes);
    }
    if (fht.FITS_imagetype()==fits_header_conventions::FLOATIMG) {
        return CcdImage<double>(CcdImage<float>(filename,header_begin_offset_in_bytes));
    }

    throw std::logic_error("Image<double>::load_image_from_unknown_filetype - "
                           "Need to implement more loading types");
}

template<typename T>
CcdImage<T> CcdImage<T>::load_from_cube_FITS_file(const std::string& filename,
        const size_t z_index)
{
    CcdImage img;
    std::ifstream infile(filename.c_str(), std::ios::binary);
    if (!infile.is_open()) { throw std::runtime_error(filename+" won't open"); }
    size_t header_offset =0;
    FitsHeader fht(infile, header_offset);
    img.load_from_cube_FITS_stream(fht, infile, header_offset, z_index);

    infile.close();
    return img;
}

template<typename T>
CcdImage<T> CcdImage<T>::load_image_from_buffered_data(
    const FileBuffer& buf,
    const FitsHeader& preloaded_fht,
    const size_t data_begin_byte_offset)
{

    CcdImage img;
    img.CCD_grid.load_from_fht(preloaded_fht);
    img.pix.hdr = PixelArray2dHeader(preloaded_fht);
    img.pix.load_data_from_buffer(buf, data_begin_byte_offset,
                                  ArrayCompressionInfo(preloaded_fht));
    return img;
}

//template<typename T>
//Image<T> Image<T>::load_image_from_buffered_cube_slice(
//        const FileBuffer& buf,
//        const FitsHeader& preloaded_fht,
//        const PixelArray2dHeader& preloaded_array_header,
//        const size_t data_begin_byte_offset ){
//
//    Image img;
//    img.load_grid_info_from_fht(preloaded_fht);
//
//    img.PixelArray2d<T>::set_header_to(preloaded_array_header);
//
//    img.PixelArray2d<T>::load_data_from_buffer(buf,
//            data_begin_byte_offset,
//            ArrayCompressionInfo(preloaded_fht) );
//
//    return img;
//}

template<typename T>
CcdImage<T> CcdImage<T>::sub_image(const CcdImage& full_img, const CcdBoxRegion& img_rgn)
{
    return CcdImage::sub_image(full_img,
                               full_img.CCD_grid.corresponding_pixel_region(img_rgn));
}

template<typename T>
CcdImage<T> CcdImage<T>::sub_image(const CcdImage<T>& full_img,
                                   const PixelBoxRegion& img_rgn)
{
    return CcdImage::sub_image(full_img, img_rgn.bounded_pixels());
}

template<typename T>
CcdImage<T> CcdImage<T>::sub_image(const CcdImage<T>& full_img, const PixelRange& img_rgn)
{
    CcdImage sub_img;
    sub_img.pix = (PixelArray2d<T>::sub_array(full_img.pix, img_rgn));

    if (full_img.CCD_grid.is_initialized()) {
        sub_img.CCD_grid =
            ImageGrid<coordinate_types::ccd> (
                full_img.CCD_grid.corresponding_grid_region(PixelBoxRegion::pixel_box_outline(img_rgn)),
                full_img.CCD_grid.pixel_width_
            );
    }

//
//    if (full_img.mosaic_grid_initialized()) {
//        sub_img.mosaic_grid_ =
//            ImageGrid<coordinate_types::mosaic> (
//                full_img.mosaic_grid().corresponding_grid_region(PixelBoxRegion::pixel_box_outline(img_rgn)),
//                full_img.mosaic_grid().pixel_width_
//            );
//    }


    return sub_img;
}

template<typename T>
void CcdImage<T>::write_to_file(const std::string& filename,
                                FitsHeader fht,
                                const ArrayCompressionInfo aci) const
{
    std::ofstream outfile(filename.c_str(), std::ios::binary);
    if (!outfile.is_open()) {
        throw std::runtime_error(
            "Image<T>::write_to_file - error writing to file "+ filename);
    }

    write_header_to_fht(fht, aci);
    outfile<<fht;
    write_data_to_stream(outfile,aci);
    outfile.close();
}


template<typename T>
void CcdImage<T>::write_header_to_fht(FitsHeader& fht,
                                      const ArrayCompressionInfo aci) const
{
    pix.write_header_to_fht(fht, aci);

    //Create a pixel_region matching the array, output to the header so that
    //DS9 shows "physical" co-ords matching our "pixel" cartesian coords
    ImageGrid<coordinate_types::pixels> pixel_grid(
        PixelBoxRegion::pixel_box_outline(pix.range()), 1.0);
    pixel_grid.write_to_fht(fht);
    CCD_grid.write_to_fht(fht);
}

template<typename T>
void CcdImage<T>::write_data_to_stream(std::ostream& os,
                                       const ArrayCompressionInfo& aci) const
{
    pix.write_data_to_stream(os, aci);
}

template<typename T>
void CcdImage<T>::load_from_stream(const FitsHeader& preloaded_fht,
                                   std::istream& is, const size_t header_offset)
{
    pix.hdr = PixelArray2dHeader(preloaded_fht);
    pix.seek_and_load_data_from_stream(is,
                                       header_offset+preloaded_fht.header_file_length_in_bytes(),
                                       ArrayCompressionInfo(preloaded_fht));
    CCD_grid.load_from_fht(preloaded_fht);
}

template<typename T>
void CcdImage<T>::load_from_cube_FITS_stream(const FitsHeader& preloaded_fht,
        std::istream& is,
        const size_t hdr_begin_byte_offset,
        const size_t z_index
                                            )
{

    CCD_grid.load_from_fht(preloaded_fht);

    PixelCubeHeader cube_header(preloaded_fht);
    if (z_index < 1 || z_index > cube_header.z_dim()) {
        throw std::runtime_error("Invalid Z-index for this cube file.");
    }

    pix.hdr =cube_header.get_xy_array2d_header();

    size_t data_step_per_plane = pix.n_pixels() * sizeof(T);
    pix.seek_and_load_data_from_stream(is,
                                       hdr_begin_byte_offset + preloaded_fht.header_file_length_in_bytes() +
                                       data_step_per_plane*(z_index-1),
                                       ArrayCompressionInfo(preloaded_fht));

}

//template<typename T>
//void CcdImage<T>::load_grid_info_from_fht(const FitsHeader& fht)
//{
//    if (fht.key_exists("CCD_RGN")) {
//        CCD_gridinitialized_=true;
//        CCD_grid.load_from_fht(fht);
//    } else(CCD_gridinitialized_=false);
//
//    if (fht.key_exists("MOS_RGN")) {
//        mosaic_grid_initialized_=true;
//        mosaic_grid_.load_from_fht(fht);
//    } else(mosaic_grid_initialized_=false);
//
//}

template<typename T>
void CcdImage<T>::initialize_CCD_grid_for_raw_data()
{
    if (CCD_grid.is_initialized())
        throw std::runtime_error("Image<T>::initialize_CCD_grid_for_raw_data() -"
                                 " CCD_grid already initialized");

    CCD_grid=ImageGrid<coordinate_types::ccd>(pix.range(), CcdPosition(0.,0.),
             1.0);

}

template<typename T>
void CcdImage<T>::initialize_CCD_grid_to_specific_region(
    const CcdBoxRegion& ccd_region)
{
    if (CCD_grid.is_initialized())
        throw std::runtime_error("Image<T>::initialize_CCD_grid_to_specific_layout() -"
                                 " CCD_grid already initialized");

    double deduced_pixel_x_scale = ccd_region.x_dim() / (double)pix.range().x_dim();
    //assert (square pixels)
    assert(deduced_pixel_x_scale == ccd_region.y_dim() / (double)pix.range().y_dim());
    CCD_grid = ImageGrid<coordinate_types::ccd>(ccd_region, deduced_pixel_x_scale);
}



template<typename T>
void CcdImage<T>::initialize_CCD_grid_to_specific_offset_and_scale(
    const CcdPosition& ccd_low_corner, const double pixel_width_in_CCD_coords)
{
    if (CCD_grid.is_initialized())
        throw std::runtime_error("Image<T>::initialize_CCD_grid_to_specific_scale() -"
                                 " CCD_grid already initialized");
    CCD_grid = ImageGrid<coordinate_types::ccd>(pix.range(),ccd_low_corner,
               pixel_width_in_CCD_coords);
}

//Supported types:
//NB THIS MUST GO AT THE END OF THE FILE
//Otherwise the compiler tries to instantiate the templatized classes before
//it has all the relevant information.
template class CcdImage<float>;
template class CcdImage<double>;

template class CcdImage<int>;
template class CcdImage<uint16_t>;
template class CcdImage<uint32_t>;

//=========================================================================================================
}//end namespace coela

