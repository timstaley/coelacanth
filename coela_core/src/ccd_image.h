/*
 * File:   ccd_image.h
 * Author: ts337
 *
 * Created on 15 February 2011, 11:28
 */

#ifndef COELA_CCD_IMAGE_H
#define COELA_CCD_IMAGE_H

#include "pixel_array2d.h"
#include "ImageGrid.h"

namespace coela {
//=========================================================================================================
// CCDImage:
/** Ties together a pixel_array with associated spatial regions.
Provides:
    - Convenience routines to load both pixel data and 'CCDGrid' data from a FITS file.
    - Co-ordinate to pixel conversion and vice versa
    - Sub-image creation with appropriate spatial region assignment.
*/
template<typename data_type>
struct CCDImage {
public:
    //------------------------------------------------------------------------
    //public Data:
    //------------------------------------------------------------------------
    ///Array member is made public
    PixelArray2d<data_type> pix;
    ImageGrid<coordinate_types::CCD> CCD_grid;
    //-------------------------------------------------------------------------
    //Constructors & input / output
    //-------------------------------------------------------------------------
    CCDImage();///< Null init.

    //initialize from different data type
    //Explicit, since it may involve loss of precision.
    template<typename U>
    explicit CCDImage(const CCDImage<U>& lower_precision_data_type_image);

    ///Creating sub-image via CCD / frame / pixel box specs.
    /// This manages the PixelArray and also updates co-ord grids.
    static CCDImage sub_image(const CCDImage&, const CCD_BoxRegion& img_rgn);
    static CCDImage sub_image(const CCDImage&, const PixelBoxRegion& img_rgn);
    static CCDImage sub_image(const CCDImage&, const PixelRange& img_rgn);

    ///File load constructor / write to file (just wrappers to stream funcs)
    CCDImage(const std::string& filename,
             const size_t header_begin_offset_in_bytes=0);

    void write_to_file(const std::string& filename,
                       FitsHeader additional_header_info = FitsHeader(),
                       const ArrayCompressionInfo=ArrayCompressionInfo::no_compression()
                      ) const;

    static CCDImage load_from_unknown_filetype(const std::string& filename,
            const size_t header_begin_offset_in_bytes=0);

    static CCDImage load_image_from_buffered_data(
        const FileBuffer&,
        const FitsHeader& preloaded_header,
        const size_t data_begin_byte_offset);


    ///NB!!! z_index = 1,2,3.... N   (not 0,1,...,N-1) in accordance with FITS convention
    static CCDImage load_from_cube_FITS_file(const std::string& filename,
            const size_t z_index);

//    static image load_image_from_buffered_cube_slice(
//            const FileBuffer&,
//            const FitsHeader& preloaded_header,
//            const PixelArray2dHeader& preloaded_array_header,
//            const size_t data_begin_byte_offset
//            );

    //File load from buffer:

    //Decoupled header and array write to stream:
    void write_header_to_fht(FitsHeader& fht,
                             const ArrayCompressionInfo = ArrayCompressionInfo::no_compression()) const;
    void write_data_to_stream(std::ostream&, const ArrayCompressionInfo&) const;

    //Load from stream + preloaded FitsHeader
    void load_from_stream(const FitsHeader& preloaded_fht, std::istream& is,
                          const size_t hdr_begin_byte_offset);


    //NB for a cube of N x*y planes, we regard z_index as belonging to the set
    // 1,2,....,N inclusive;
    //i.e. we conform to 1-origin indexing.
    void load_from_cube_FITS_stream(const FitsHeader& preloaded_fht,
                                    std::istream& is,
                                    const size_t hdr_begin_byte_offset,
                                    const size_t z_index
                                   );

    //-----------------------------------------------------------------------
    //co-ordinate grid access

    //Setter funcs:
    void initialize_CCD_grid_for_raw_data();
    void initialize_CCD_grid_to_specific_region(const CCD_BoxRegion& ccd_region);

    void initialize_CCD_grid_to_specific_offset_and_scale(
        const CCD_Position& ccd_low_corner,
        const double pixel_width_in_CCD_coords);

};

//===========================================================================
}//end namespace coela


#endif  /* CCD_IMAGE_H */

