#include "../drizzle_subroutines.h"
#include "../image_cleanup.h"
#include <iostream>
#include <boost/filesystem.hpp>
#include <stdexcept>
namespace coela {
namespace pipeline {

using std::cout; using std::endl;
void create_quick_drizzle_mosaic(const std::vector< MosaicImage<double> >& img_sums,
                                 const std::vector< MosaicImage<double> >& img_weights,
                                 MosaicImage<double>& full_output_sum,
                                 MosaicImage<double>& full_output_weights)
{
    assert(img_sums.size() == img_weights.size());

    full_output_sum.pix.assign(0.0);
    full_output_weights.pix.assign(0.0);

    for (size_t i=0; i!=img_sums.size(); ++i) {
        const MosaicImage<double>& frm_sum = img_sums[i];
        const MosaicImage<double>& frm_weights = img_weights[i];

        double x_offset = frm_sum.mosaic_grid.image_outline_.low.x -
                          full_output_sum.mosaic_grid.image_outline_.low.x;
        double y_offset = frm_sum.mosaic_grid.image_outline_.low.y -
                          full_output_sum.mosaic_grid.image_outline_.low.y;

        assert(fmod(x_offset,1.0)==0)  ;
        assert(fmod(y_offset,1.0)==0)  ;
        assert(frm_sum.mosaic_grid.pixel_width_ ==
               full_output_sum.mosaic_grid.pixel_width_);
        assert(full_output_sum.mosaic_grid.image_outline_.contains_region(
                   frm_sum.mosaic_grid.image_outline_));

        PixelIndex output_pixel_offset(x_offset, y_offset);

        for (PixelIterator it(frm_sum.pix.range()); it!=it.end; ++it) {
            full_output_sum.pix(it+output_pixel_offset)+=frm_sum.pix(it);
            full_output_weights.pix(it+output_pixel_offset)+=frm_weights.pix(it);
        }
    }

    return;
}


//=======================================================================================
MosaicBoxRegion get_mosaic_region(
    const vector<MosaicBoxRegion>& frame_output_regions)
{
    MosaicBoxRegion full_mosaic_region(frame_output_regions.front());

    for (size_t i=0; i!=frame_output_regions.size(); ++i) {
        full_mosaic_region.enlarge_to_cover_region(frame_output_regions[i]);
    }

    full_mosaic_region.expand_to_pixel_boundaries();
    return full_mosaic_region;
}

MosaicImage<double> init_blank_mosaic_frame(const MosaicBoxRegion& mosaic_region,
        const double mosaic_frame_pixel_width)
{
    assert(fmod(mosaic_region.low.x, 1.0)==0.0);
    assert(fmod(mosaic_region.low.y, 1.0)==0.0);
    assert(fmod(mosaic_region.high.x, 1.0)==0.0);
    assert(fmod(mosaic_region.high.y, 1.0)==0.0);

    MosaicImage<double> mosaic;
    mosaic.pix = PixelArray2d<double>(mosaic_region.x_dim() / mosaic_frame_pixel_width,
                                      mosaic_region.y_dim()/ mosaic_frame_pixel_width,
                                      0.0);
    mosaic.initialize_mosaic_grid_to_specific_region(mosaic_region);
    return mosaic;

}

//=======================================================================================
CcdDrizzleData CcdDrizzleData::prep_drizzle_info(
    const CcdDatasetInfo& cdi,
    const DrizzleSettings& ds,
    const CameraConfigInfo& cam_conf,
    const string& drizzle_dir
)
{

    CcdDrizzleData drz_inf(cdi.ccd_id);
    CcdCalibrationInfo ccd_chars = cam_conf.get_calibration_info_for_CCD_id(cdi.ccd_id);

    ////========================   Load information abour mosaic layout, determine and initialise output frame  ===============================
    drz_inf.input_CCD_region = ccd_chars.crop_region;

    //-----------------------------------------------------------------------------------------------
    //Figure out the input CCD's mosaic placement:

    MosaicPixelShift CCD_mos_offset =cam_conf.get_mosaic_offset_for_ccd_id(cdi.ccd_id);
    drz_inf.input_mosaic_region = MosaicBoxRegion(CCD_mos_offset+MosaicPosition::origin,
                                  CCD_mos_offset+MosaicPosition::origin);
    //input_mosaic_region is now a null box at the mosaic Position corresponding to CcdPosition(0,0)

    drz_inf.input_mosaic_region.low +=
        MosaicPixelShift(ccd_chars.crop_region.low.x ,
                         ccd_chars.crop_region.low.y)/ ds.drizzle_scale_factor ;

    drz_inf.input_mosaic_region.high +=
        MosaicPixelShift(ccd_chars.crop_region.high.x ,
                         ccd_chars.crop_region.high.y) /ds.drizzle_scale_factor;

    cout<<"Prepping CCD "<<cdi.ccd_id<<", "
        <<"CCD offset: " << CCD_mos_offset<<", "
        <<"input mosaic region: " <<drz_inf.input_mosaic_region<<endl;

    drz_inf.output_mosaic_region = drz_inf.input_mosaic_region;
    drz_inf.output_mosaic_region.expand_to_pixel_boundaries();
    //Align pixel boundaries with integer detector values
    //This allows easy alignment of separately drizzled CCDs for quick concatenation
    cout<<"Prepping CCD "<<cdi.ccd_id<<", output mosaic region: "
        <<drz_inf.output_mosaic_region<<endl;


    drz_inf.drizzled_output_PixelRange = PixelRange(1,1,
                                         drz_inf.output_mosaic_region.x_dim(),
                                         drz_inf.output_mosaic_region.y_dim());

    MosaicPixelShift input_to_output_corner_mos_shift(
        drz_inf.output_mosaic_region.low - drz_inf.input_mosaic_region.low);

    CcdPixelShift input_to_output_corner_CCD_shift(
        input_to_output_corner_mos_shift.x * ds.drizzle_scale_factor,
        input_to_output_corner_mos_shift.y * ds.drizzle_scale_factor);

    drz_inf.output_CCD_low_corner =  drz_inf.input_CCD_region.low +
                                     input_to_output_corner_CCD_shift;

    //===================== load in the crop region, default frame weights: ===============================
    drz_inf.raw_data_crop_box = ccd_chars.cropped_PixelRange;

//    cropped_hdr.write_to_file("cropped_rgn.fits");
    drz_inf.default_frame_weights=image_cleanup::get_CCD_default_weight_map(ccd_chars);
    drz_inf.default_frame_weights.initialize_mosaic_grid_to_specific_region(
        drz_inf.input_mosaic_region);
    drz_inf.default_frame_weights.write_to_file(drizzle_dir+"default_weights_CCD"
            +string_utils::itoa(cdi.ccd_id)+".fits");

////========================   Load gain information if required for normalisation and frame to frame pedestal determination ===============================
    if (ds.normalisation_on) {
        map<int,long> data_histogram;

        string histogram_filename = cdi.most_recently_output_histogram;
        if (boost::filesystem::exists(histogram_filename)) {
            data_histogram= gain_utils::load_histogram_data_from_file(histogram_filename);
        }

        else throw std::runtime_error("No histogram found for CCD " +
                                          string_utils::itoa(cdi.ccd_id)+" at "+ histogram_filename);
        drz_inf.gain_inf= gain_utils::fit_CCD_histogram(data_histogram);
        cout<<"Gain inf for CCD " << cdi.ccd_id<<":"<<drz_inf.gain_inf<<endl;
    }

//
////========================  Load static debias frame ============================================
//    if (ls.get_info_for_CCD_id(ccd_id).column_debias_img_path!="none"){
//        drz_inf.bias_frame_available=true;
    if (!cam_conf.simulated_data) {
        cout<<"Loading column bias frame from:\n "
            <<cdi.most_recently_used_column_bias_frame<<endl;
        drz_inf.combined_bias_frame = CcdImage<float>(cdi.most_recently_used_column_bias_frame);

        if (ccd_chars.precal_row_bias_frame_available) {
            cout<<"Loading row bias frame from:\n "
                <<ccd_chars.precal_row_bias_frame_path<<endl;
            CcdImage<float> row_bias(ccd_chars.precal_row_bias_frame_path);
            drz_inf.combined_bias_frame.pix += row_bias.pix;
        }
    }
//    }
//    else drz_inf.bias_frame_available=false;
//
////========================  Prepare dark current correction             ============================================

    drz_inf.dark_current_map_available = ccd_chars.dark_current_frames_available;

    if (ccd_chars.dark_current_frames_available && ds.normalisation_on) {
        drz_inf.dark_current_photon_excess=         CcdImage<float>
                (ccd_chars.normalised_DC_frame_path);
        drz_inf.dark_current_thresholded_excess=    CcdImage<float>
                (ccd_chars.thresholded_DC_frame_path);
        cout<<"Loaded dark current calibration frames from:\n "
            <<ccd_chars.normalised_DC_frame_path<<" and\n"
            <<ccd_chars.thresholded_DC_frame_path<<endl;

    }

    return drz_inf;
}

//-------------------------------------------------------------------------------------------------------

const CcdDrizzleData&
CcdDrizzleData::find_drizzle_info_for_ccd_id(
    int ccd_id, const vector<CcdDrizzleData>& drz_vec)
{
    for (size_t i=0; i!=drz_vec.size(); ++i) {
        if (drz_vec[i].CCD_id==ccd_id) { return drz_vec[i]; }
    }
    throw std::runtime_error("CCD id does not match any in CCD_drizzle_info vector: "
                             + string_utils::itoa(ccd_id));
}

//-------------------------------------------------------------------------------------------------------

vector<CcdDrizzleData>
CcdDrizzleData::prep_drizzle_data_vec(const DrizzleSettings& ds,
        const vector<CcdDatasetInfo> & datasets,
        const CameraConfigInfo& cam_conf,
        const MultiFrame& first_multi_frame,
        const string& drizzle_dir
                                               )
{
    using namespace string_utils;
    vector<CcdDrizzleData> drz_inf_vec;

    for (size_t i=0; i!=first_multi_frame.synchronized_CCD_frames.size(); ++i) {
        const int ccd_id = first_multi_frame.synchronized_CCD_frames[i].ccd_id;
        cout<<"-------------------------------------------"<<endl;
        cout <<"Prepping drizzle info for ccd_id " <<ccd_id<<endl;

        //Now need to grab the relevant dataset info for this CCD id.
        //Currently implemented assuming we are only processing one run at a time.
        vector<CcdDatasetInfo> datasets_for_this_CCD_id;
        for (size_t i=0; i!=datasets.size(); ++i) {
            if (datasets[i].ccd_id == ccd_id) {
                datasets_for_this_CCD_id.push_back(datasets[i]);
            }
        }
        if (datasets_for_this_CCD_id.size()!=1) {
            throw std::runtime_error("CcdDrizzleData::prep_drizzle_data_vec --\n"
                                     "Error - none or more than one datasets with this CCD id");
        }

        CcdDrizzleData drz_inf=
            CcdDrizzleData::prep_drizzle_info(datasets_for_this_CCD_id.front(),
                    ds,
                    cam_conf,
                    drizzle_dir);

        drz_inf_vec.push_back(drz_inf);
    }
    return drz_inf_vec;
}

//-------------------------------------------------------------------------------------------------------
MosaicBoxRegion CcdDrizzleData::get_mosaic_region(
    const vector<CcdDrizzleData>& drz_inf_vec)
{
    vector<MosaicBoxRegion> frm_output_rgns;
    for (size_t i=0; i!=drz_inf_vec.size(); ++i) {
        frm_output_rgns.push_back(drz_inf_vec[i].output_mosaic_region);
    }
    return pipeline::get_mosaic_region(frm_output_rgns);
}


}//namespace
}//namespace
