/*
 * File:   multi_ccd_filters.h
 * Author: ts337
 *
 * Created on 14 March 2011, 19:47
 */

#ifndef COELA_MULTI_CCD_FILTERS_H
#define COELA_MULTI_CCD_FILTERS_H

#include "coela_utility/src/file_buffer.h"
#include "camera_config_info.h"
#include "drizzle_settings.h"
#include "coela_core/src/ccd_image.h"
#include "gain_utils.h"
#include "frame_info.h"
#include "frame_registration.h"
#include "single_ccd_filters.h"
#include "drizzle_subroutines.h"

#include <tbb/pipeline.h>
#include <map>
#include <boost/filesystem.hpp>


namespace coela {
namespace pipeline {

//============================================================================================================================
namespace multi_CCD_filters {

//============================================================================================================================
struct DrizzleToken: public single_CCD_filters::FrameCleanupToken {
    DrizzleToken():drizzle_inf(INT_MAX) {clear();} //Set drizzle inf CCD id to something silly so it will get properly initialized on first run through the filters


    //These are placeholders to be reassigned.
    CcdImage<float> drizzle_input_weights;
    MosaicImage<float>
    drizzle_result_vals; //NB needs to be divided by weights image before it resembles the sky image
    CcdImage<float> drizzle_result_weights;

    CcdImage<float> thresholded_input;
    CcdImage<float> thresholded_drizzle_result_vals;

    //============
    //These get reset by clear():
    CcdPixelShift atmospheric_translation_shift;
    //=================

    //This does not get reset, since it's static per CCD (IIRC) and is reused if the token loads another frame with same CCD id.
    CcdDrizzleData drizzle_inf;

    void clear();

    static vector<DrizzleToken> initialize_token_vec(const size_t n_tokens);
};
//============================================================================================================================
class SelectiveLoadToBufferFilter: public tbb::filter {
public:
    SelectiveLoadToBufferFilter(
        const vector<FrameInfo>& GS_frame_list_,
        const vector<MultiFrame>& multi_frames_list_,
        const vector<CcdDrizzleData>& drz_inf_vec_,
        const CcdPosition& GS_target_Position
    );
    void* operator()(void*);
    void set_selection_limits(double low_limit, double high_limit);
    double upper_limit()const {return upper_selection_limit;}
    double lower_limit()const {return low_selection_limit;}
    size_t bytes_loaded() {return bytes_off_disk;}
private:
    void reset();//Called when setting selection limits

    const vector<FrameInfo>& GS_frame_vec;
    const vector<MultiFrame> multi_frames_vec;
    const vector<CcdDrizzleData> drz_inf_vec;
    const CcdPosition GS_tgt_CCD_posn;
//        const CcdBoxRegion GS_rgn;

    double low_selection_limit, upper_selection_limit;
    vector<FrameInfo>::const_iterator current_list_Position;
    size_t next_buffer, bytes_off_disk, CCDs_loaded_for_current_exposure;
    vector<DrizzleToken> token_buffers;
    MultiFrame current_exposure_frames;
};


//============================================================================================================================
class MultiframeCropDebiasFilter: public tbb::filter {
public:
    MultiframeCropDebiasFilter():filter(parallel) {}
    void* operator()(void*);
};
//============================================================================================================================
class NormalizationFilter: public tbb::filter {
public:
    NormalizationFilter():filter(parallel) {}
    void* operator()(void*);
};
//============================================================================================================================
class PhotonThresholdingFilter: public tbb::filter {
public:
    PhotonThresholdingFilter(const double threshold_level_in_photo_electrons):
        filter(parallel), threshold_in_photo_e_(threshold_level_in_photo_electrons) {}
    void* operator()(void*);
private:
    const double threshold_in_photo_e_;
};
//============================================================================================================================
class DarkCurrentSubtractionFilter: public tbb::filter {
public:
    DarkCurrentSubtractionFilter(bool thresholding_on):filter(parallel),
        thresholding(thresholding_on) {}
    void* operator()(void*);
private:
    const bool thresholding;
};
//============================================================================================================================
class CosmicRayDownweightingFilter: public tbb::filter {
public:
    CosmicRayDownweightingFilter(int pixel_padding_radius):filter(parallel),
        padding_radius(pixel_padding_radius) {}
    void* operator()(void*);
private:
    const int padding_radius;
};
//============================================================================================================================
class DrizzleFilter: public tbb::filter {
public:
    DrizzleFilter(const double drizzle_pixel_scale,
                   const double drizzle_pixel_fraction,
                   const bool thresholding_on_):
        filter(parallel),
        pixel_scale(drizzle_pixel_scale),pixel_frac(drizzle_pixel_fraction),
        thresholding_on(thresholding_on_) {}

    void* operator()(void*);
private:
    const double pixel_scale, pixel_frac;
    const bool thresholding_on;
};

//============================================================================================================================
//RawDataMosaicFilter:
/**
Only needed if 'raw data mosaics' are being produced, i.e. we are producing
a single mosaic image representing data from all 4 detectors in
one single exposure.

Could be parallelized further, but there is no need - it was a one-off special
to produce a big raw dataset for analysis in Matlab.
 */
class RawDataMosaicFilter: public tbb::filter {
public:
    RawDataMosaicFilter(const vector<CcdDrizzleData>& drz_inf_vec,
                           const string& output_folder,
                           const MosaicBoxRegion& full_output_region);

    void* operator()(void*);

private:
    MosaicImage<float> raw_mosaic_vals, raw_mosaic_weights;
    map<int,bool> CCDs_drizzled;
    const string output_folder;
    bool output_folder_created;
};

//============================================================================================================================
class OutputSummationFilter: public tbb::filter {
public:
    OutputSummationFilter(const vector<CcdDrizzleData>& drz_inf_vec,
                            const bool thresholding_on_);
    CcdImage<double> drizzled_image_for_CCD(const int ccd_id);
    CcdImage<double> drizzled_vals_for_CCD(const int ccd_id);
    CcdImage<double> drizzled_threshed_vals_for_CCD(const int ccd_id);
    CcdImage<double> drizzled_weights_for_CCD(const int ccd_id);
    void* operator()(void*);
    void reset();
private:
    bool ccd_id_valid(const int ccd_id) const;
    size_t get_vec_Position_for_ccd_id(const int ccd_id) const;
    const bool thresholding_on;
    vector< CcdImage<double> > summed_CCD_drizzle_val_imgs, summed_CCD_drizzle_weight_imgs;
    vector< CcdImage<double> > summed_threshed_vals;
//        vector<int> index_vec;
    std::map<int, size_t> index;

};

//============================================================================================================================
//============================================================================================================================
}//end namespace coela::multi_CCD_filters
}//end namespace
}//end namespace coela

#endif  /* MULTI_CCD_FILTERS_H */

