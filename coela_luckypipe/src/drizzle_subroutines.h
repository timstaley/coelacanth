/*
 * File:   drizzle_subroutines.h
 * Author: ts337
 *
 * Created on 08 July 2011, 16:21
 */

#ifndef COELA_DRIZZLE_SUBROUTINES_H
#define COELA_DRIZZLE_SUBROUTINES_H
#include "coela_core/src/mosaic_image.h"
#include "gain_utils.h"
#include "ccd_dataset_info.h"
#include "drizzle_settings.h"
#include "camera_config_info.h"


namespace coela {
namespace pipeline {

void create_quick_drizzle_mosaic(const std::vector< MosaicImage<double> >& img_sums,
                                 const std::vector< MosaicImage<double> >& img_weights,
                                 MosaicImage<double>& full_output_sum,
                                 MosaicImage<double>& full_output_weights);

//=======================================================================================
MosaicBoxRegion get_mosaic_region(
    const std::vector<MosaicBoxRegion>& frame_output_regions);

MosaicImage<double> init_blank_mosaic_frame(const MosaicBoxRegion& mosaic_region,
        const double mosaic_frame_pixel_width);
//=======================================================================================
struct CcdDrizzleData {
    CcdDrizzleData(int CCD_id_):CCD_id(CCD_id_) {}
    int CCD_id;

    PixelRange raw_data_crop_box;

    CcdImage<float> combined_bias_frame;

    gain_utils::GainData gain_inf;

    //long_exposure_avg is normalised to photon levels
    MosaicImage<float> default_frame_weights;

    bool dark_current_map_available;
    CcdImage<float> dark_current_photon_excess,  dark_current_thresholded_excess;

    CcdBoxRegion input_CCD_region;
    CcdPosition output_CCD_low_corner;

    MosaicBoxRegion input_mosaic_region;
    MosaicBoxRegion output_mosaic_region;

    PixelRange drizzled_output_PixelRange;

    static CcdDrizzleData prep_drizzle_info(const CcdDatasetInfo&,
            const DrizzleSettings& ds,
            const CameraConfigInfo& cam_conf,
            const std::string& drizzle_dir);

    static const
    CcdDrizzleData& find_drizzle_info_for_ccd_id(int ccd_id,
            const vector<CcdDrizzleData>& drz_vec);

    static
    vector<CcdDrizzleData> prep_drizzle_data_vec(const DrizzleSettings& ds,
            const vector<CcdDatasetInfo> & datasets,
            const CameraConfigInfo& cam_conf,
            const MultiFrame& first_multi_frame,
            const std::string& drizzle_dir);

    static MosaicBoxRegion get_mosaic_region(
        const vector<CcdDrizzleData>& drz_inf_vec);
};
//=======================================================================================



}//namespace
}//namespace


#endif  /* DRIZZLE_SUBROUTINES_H */

