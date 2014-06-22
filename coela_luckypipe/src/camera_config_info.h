/*
 * File:   camera_config_info.h
 * Author: ts337
 *
 * Created on 15 March 2011, 13:01
 */

#ifndef COELA_CAMERA_CONFIG_INFO_H
#define COELA_CAMERA_CONFIG_INFO_H

#include "coela_core/src/cartesian_coords.h"
#include "ccd_dataset_info.h"

#include <string>
#include <vector>
using std::string;
using std::vector;

namespace coela {

//============================================================================================================
struct CcdCalibrationInfo {

    //Constructors
    CcdCalibrationInfo():
        precal_row_bias_frame_available(false), dark_current_frames_available(false),
        precal_row_bias_frame_path("none"), normalised_DC_frame_path("none"),
        thresholded_DC_frame_path("none"),
        milliseconds_per_frame_exposure(0) {}

    static vector<CcdCalibrationInfo> load_from_text_vec(const vector<string>&);
    friend std::ostream& operator<<(std::ostream& os,
                                    const vector<CcdCalibrationInfo>& ccd_prop);
    //-------------------------------------------------------
    //data members:
    int ccd_id;
    vector<size_t> bad_columns;
    PixelRange cropped_PixelRange;
    CcdBoxRegion crop_region;

    CcdBoxRegion default_temporal_debiasing_histogram_region;

    bool precal_row_bias_frame_available, dark_current_frames_available;
    string precal_row_bias_frame_path;
    string normalised_DC_frame_path, thresholded_DC_frame_path;

    double milliseconds_per_frame_exposure;
    //-------------------------------------------------------

    //Member functions
    //NB ignores anything outside the crop_region since that is all assumed to be rubbish anyway
    vector<CcdBoxRegion> get_bad_detector_regions() const;
};

//============================================================================================================

//TO DO - use this structure to include scale.
//struct CCD_mosaic_information{
//    int ccd_id;
//    double pixel_scale;
//    MosaicPixelShift mosaic_offset_for_CCD_origin;
//};

struct LensAndApertureInfo {
    void load_from_text_vec(const vector<string>& file_text);
    friend std::ostream& operator<<(std::ostream& os, const LensAndApertureInfo& lens_inf);

    //-------------------------------------------------------
    //data
    double telescope_outer_diameter, telescope_inner_diameter, nominal_pixel_scale_in_mas;
    ///pairs of (CCD id, offset)
    vector< std::pair<int, MosaicPixelShift> > ccd_mosaic_offsets;
    //-------------------------------------------------------

};


//============================================================================================================

///Represents an optical filter
///(As opposed to a data filter)
struct OpticalFilterInfo {
    static vector<OpticalFilterInfo> load_from_text_vec(const vector<string>& file_text);
    friend std::ostream& operator<<(std::ostream& os,
                                    const vector<OpticalFilterInfo>& filter_set);

    //-------------------------------------------------------
    int ccd_id;
    double central_wavelength_in_metres;
    string name;
    //-------------------------------------------------------
};

//============================================================================================================

class CameraConfigInfo {
public:
    CameraConfigInfo() {}
    CameraConfigInfo(const string& filename);

    void write_to_files(const std::string& camconf_filename,
                        std::string CCD_inf_rel_path=string(),
                        std::string lens_inf_rel_path=string(),
                        std::string filter_inf_rel_path=string());


    CcdCalibrationInfo& get_calibration_info_for_CCD_id(int ccd_id);
    const CcdCalibrationInfo& get_calibration_info_for_CCD_id(int ccd_id) const;

    OpticalFilterInfo& get_filter_info_for_ccd_id(int ccd_id);
    const OpticalFilterInfo& get_filter_info_for_ccd_id(int ccd_id) const;

    const MosaicPixelShift& get_mosaic_offset_for_ccd_id(int ccd_id) const;

    //-------------------------------------------------------
    //data
    bool simulated_data; ///< Informs that cleanup routines are not necessary
    std::string self_path; ///< Used for calculating paths of subfiles for ccd,lens,filter...

    vector<OpticalFilterInfo> filters;
    vector<CcdCalibrationInfo> CCD_vec;
    LensAndApertureInfo lens_inf;
    //-------------------------------------------------------

private:
    void load_from_text_vec(const vector<string>& file_text);
};

//============================================================================================================
}//end namespace coela


#endif  /* CAMERA_CONFIG_INFO_H */

