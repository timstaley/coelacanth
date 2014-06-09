#ifndef COELA_MOSAIC_IMAGE_H
#define COELA_MOSAIC_IMAGE_H

#include "ccd_image.h"

namespace coela {

//MosaicImage:
/**
 A simple extension of CcdImage to bundle a 'mosaic' frame of reference.

Provides a couple of convenience routines to streamline handling of the
additional ImageGrid.

 */


template<typename data_type>
struct MosaicImage: public CCDImage<data_type> {
public:
    //------------------------------------------------
    //Data:
    //Inherits pix, CCD_grid from CCDImage
    ImageGrid<coordinate_types::mosaic> mosaic_grid;
    //------------------------------------------------

    MosaicImage();///< Null init.

    ///Cast a CcdImage to a MosaicImage (with blank mosaic_grid).
    ///(This is OK, since mosaic_grid can be checked with .is_initialized ).
    MosaicImage(const CCDImage<data_type>& rhs);

    void initialize_mosaic_grid_to_specific_region(
        const MosaicBoxRegion& mosaic_region,bool no_throw=false);

    void write_to_file(const std::string& filename,
                       FitsHeader additional_header_info = FitsHeader(),
                       const ArrayCompressionInfo=ArrayCompressionInfo::no_compression()
                      ) const;
};

} //end namespace



#endif  // MOSAIC_IMAGE_H


