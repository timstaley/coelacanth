#include "../multi_ccd_filters.h"

#include "../gain_utils.h"
#include "../image_cleanup.h"
#include "coela_core/src/drizzle.h"

using std::runtime_error;
using std::logic_error;
using std::cout;
using std::cerr;
using std::endl;

namespace coela {
namespace pipeline {
const double error_margin = 1e-6;



//=======================================================================================
//=======================================================================================

namespace multi_CCD_filters {
const int threaded_drizzle_debug=0;

//=======================================================================================
void drizzle_token::clear()
{
    single_CCD_filters::frame_cleanup_token::clear();
    atmospheric_translation_shift=CCD_PixelShift();
}

vector<drizzle_token> drizzle_token::initialize_token_vec(const size_t n_tokens)
{
    vector< drizzle_token> v(n_tokens);
    for (size_t i=0; i!=v.size(); ++i) {
        v[i].token_buffer_id=i;
    }
    return v;
}
//=======================================================================================
Quality_Selected_Load_to_Buffer_Filter::Quality_Selected_Load_to_Buffer_Filter(
    const vector<FrameInfo>& GS_frame_list_,
    const vector<MultiFrame>& multi_frames_list,
    const vector<CCD_specific_drizzle_inf>& drz_inf_vec_,
    const CCD_Position& GS_target_Position):
    tbb::filter(serial_in_order),
    GS_frame_vec(GS_frame_list_),
    multi_frames_vec(multi_frames_list.begin(), multi_frames_list.end()),
    drz_inf_vec(drz_inf_vec_),
    GS_tgt_CCD_posn(GS_target_Position),
//                GS_rgn(gs_rgn_),
    current_list_Position(GS_frame_vec.begin()),
    next_buffer(0), bytes_off_disk(0) ,
    CCDs_loaded_for_current_exposure(0),
    current_exposure_frames(-1)
{
    token_buffers=drizzle_token::initialize_token_vec(single_CCD_filters::lucky_n_tokens);
    assert(GS_frame_vec.front().quality_percentile_rank!=0.0);
    low_selection_limit=upper_selection_limit=0;
}

void Quality_Selected_Load_to_Buffer_Filter::reset()
{
    current_list_Position=GS_frame_vec.begin();
    next_buffer=0;
}


void Quality_Selected_Load_to_Buffer_Filter::set_selection_limits(double low_limit,
        double high_limit)
{
    low_selection_limit=low_limit;
    upper_selection_limit=high_limit;
    reset();
}

void* Quality_Selected_Load_to_Buffer_Filter::operator()(void*)
{
    if (threaded_drizzle_debug) { cerr<<"In load buffer "<<endl; }
//            cout <<"Selection limits: " <<low_selection_limit<<" , "<<upper_selection_limit<<endl;
    if (low_selection_limit==upper_selection_limit) { throw logic_error("limits not set"); }
    drizzle_token& d = token_buffers[next_buffer];
//            d.clear();//called by decommission filter
    next_buffer= (next_buffer+1) % single_CCD_filters::lucky_n_tokens;

    if (CCDs_loaded_for_current_exposure==0) {
        while (!(current_list_Position->quality_percentile_rank > low_selection_limit
                 && current_list_Position->quality_percentile_rank <= upper_selection_limit)
                && current_list_Position != GS_frame_vec.end()) {
            current_list_Position++;
        }
        if (current_list_Position==GS_frame_vec.end()) { return NULL; }

        current_exposure_frames=MultiFrame::pull_multi_frame_with_id(
                                    current_list_Position->corrected_frame_id, multi_frames_vec);
        assert(current_exposure_frames.multi_frame_has_been_drizzled==false);
        current_exposure_frames.multi_frame_has_been_drizzled=true;
    }

    d.frame_inf=
        current_exposure_frames.synchronized_CCD_frames[CCDs_loaded_for_current_exposure];

    d.atmospheric_translation_shift =
        CCD_PixelShift(GS_tgt_CCD_posn -
                       current_list_Position->guide_star_estimates.front().Position);

    string file_extension = string_utils::pull_extension(d.frame_inf.file_path);
//    d.fht = FitsHeader(d.frame_inf.file_path, d.frame_inf.header_byte_offset);
    if (d.frame_inf.byte_size==0 && (file_extension==".fits" || file_extension==".lcc")) {
        d.buffer.buffer_whole_file(d.frame_inf.file_path);
    } else {
        d.buffer.buffer_file_portion(d.frame_inf.file_path,
                                     d.frame_inf.header_byte_offset,
                                     d.frame_inf.byte_size);
    }

    bytes_off_disk+=d.buffer.data_vec().size();
    if (d.drizzle_inf.CCD_id!= d.frame_inf.ccd_id) {
        //FIXME: This has gotten to the stage of copying several bitmaps!
        //Would be more computationally efficient (although larger memory footprint) to give each token a full copy of the drz_inf_vec
        //Then just switch an index for each frame processed.
        d.drizzle_inf= CCD_specific_drizzle_inf::find_drizzle_info_for_ccd_id(d.frame_inf.ccd_id,
                       drz_inf_vec);
    }
    d.drizzle_input_weights = d.drizzle_inf.default_frame_weights;

    CCDs_loaded_for_current_exposure++;

    if (CCDs_loaded_for_current_exposure==
            current_exposure_frames.synchronized_CCD_frames.size()) {
        current_list_Position++;
        CCDs_loaded_for_current_exposure=0;
    }
    return &d;
}

//=======================================================================================


//=======================================================================================

void* Drizzle_Token_Crop_and_Debias_Filter::operator()(void* item)
{
    if (threaded_drizzle_debug) { cerr<<"In debias filter..."; }
    drizzle_token& d = *static_cast<drizzle_token*>(item);

    //FrameCropFilter
    d.img=CCDImage<float>::sub_image(d.img, d.drizzle_inf.raw_data_crop_box);

    assert(d.img.CCD_grid.image_outline_==d.drizzle_input_weights.CCD_grid.image_outline_);

    //Subtract uniform bias
    d.img.pix-=d.frame_inf.bias_pedestal;

    //Fix col-to-col variation
    d.img.pix-= d.drizzle_inf.combined_bias_frame.pix;
    if (threaded_drizzle_debug) { cerr<<"Done\n"; }
    return &d;
}

//=======================================================================================

void* normalization_filter::operator()(void* item)
{
    drizzle_token& d = *static_cast<drizzle_token*>(item);
    gain_utils::normalise_CCD_with_uniform_gain(d.img, d.drizzle_inf.gain_inf);
    return &d;
}

//=======================================================================================

void* Photon_Thresholding_Filter::operator()(void* item)
{
    if (threaded_drizzle_debug) { cerr<<"In thresholding filter"<<endl; }
    drizzle_token& d = *static_cast<drizzle_token*>(item);

    d.thresholded_input = gain_utils::threshold_bitmap(d.img,
                          threshold_in_photo_e_);

//    assert(d.thresholded_input.key_exists("THRESHED"));
    return &d;
}
//=======================================================================================
void* Dark_Current_Subtraction_Filter::operator()(void* item)
{
    if (threaded_drizzle_debug) { cerr<<"In DC subtraction filter"<<endl; }
    drizzle_token& d = *static_cast<drizzle_token*>(item);
    if (d.drizzle_inf.dark_current_map_available) {
        d.img.pix-=d.drizzle_inf.dark_current_photon_excess.pix;
        if (thresholding) { d.thresholded_input.pix -= d.drizzle_inf.dark_current_thresholded_excess.pix; }
    }
    return &d;
}

//=======================================================================================
void* Cosmic_Ray_Downweighting_Filter::operator()(void* item)
{
    drizzle_token& dt = *static_cast<drizzle_token*>(item);

    if (!dt.frame_inf.confirmed_cosmic_rays.empty()) {
        for (size_t i=0; i!=dt.frame_inf.confirmed_cosmic_rays.size(); ++i) {
            PixelIndex p = (
                               PixelPosition::pixel_centred_at(
                                   dt.img.CCD_grid.corresponding_pixel_Position(dt.frame_inf.confirmed_cosmic_rays[i])
                               )
                           );
            PixelRange blanking_box(p,p);
            blanking_box = PixelRange::pad(blanking_box,padding_radius);
            blanking_box = PixelRange::overlap(dt.drizzle_input_weights.pix.range(), blanking_box);
            for (PixelIterator it(blanking_box); it!=it.end; ++it) {
                dt.drizzle_input_weights.pix(it)=0;
            }
        }
    }

    return &dt;
}
//=======================================================================================
void* drizzle_filter::operator()(void* item)
{
    if (threaded_drizzle_debug) { cerr<<"In drizzle filter... "; }
    drizzle_token& d = *static_cast<drizzle_token*>(item);


    d.img.initialize_mosaic_grid_to_specific_region(
        d.drizzle_inf.input_mosaic_region, true);

    MosaicPixelShift fixed_Positioning_shift =
        d.drizzle_inf.input_mosaic_region.low
        - d.drizzle_inf.output_mosaic_region.low;

    //    assert( fabs(d.img.mosaic_grid.pixel_width_ - 1.0 ) <error_margin) ; // Not true, mosaic may have smaller pixel size than input
    //This should currently always be true
    assert(fabs(d.img.CCD_grid.pixel_width_ - 1.0) < error_margin);

    PixelShift drizzle_input_shift = d.img.mosaic_grid.corresponding_pixel_shift(
                                         fixed_Positioning_shift) +
                                     d.img.CCD_grid.corresponding_pixel_shift(d.atmospheric_translation_shift);


    if (threaded_drizzle_debug) { cerr<<"Ready to drizzle... "; }
    //reset the results files in case they have the wrong size, header vals etc. for the current CCD
    d.drizzle_result_vals.pix =
        PixelArray2d<float>(d.drizzle_inf.drizzled_output_PixelRange.x_dim(),
                            d.drizzle_inf.drizzled_output_PixelRange.y_dim(),0.0);
    d.drizzle_result_vals.initialize_mosaic_grid_to_specific_region(
        d.drizzle_inf.output_mosaic_region, true);
    d.drizzle_result_weights = d.drizzle_result_vals;

    if (thresholding_on) {
//        assert(d.thresholded_input.key_exists("THRESHED"));
        d.thresholded_drizzle_result_vals=CCDImage<float>(d.drizzle_result_vals);
//        d.thresholded_drizzle_results.set_keyword("THRESHED", d.thresholded_input.get_key_value("THRESHED"));

        drizzle::dual_translate_and_drizzle_frame(d.img.pix, d.thresholded_input.pix,
                d.drizzle_input_weights.pix,
                d.drizzle_result_vals.pix, d.thresholded_drizzle_result_vals.pix,
                d.drizzle_result_weights.pix,
                drizzle_input_shift,
                pixel_scale, pixel_frac
                                                 );
    } else drizzle::translate_and_drizzle_frame(d.img.pix, d.drizzle_input_weights.pix,
                d.drizzle_result_vals.pix, d.drizzle_result_weights.pix,
                drizzle_input_shift,
                pixel_scale, pixel_frac
                                                   );
    if (threaded_drizzle_debug) { cerr<<" done."<<endl; }
    return &d;
}
//=======================================================================================

raw_data_mosaic_filter::raw_data_mosaic_filter(
    const vector<CCD_specific_drizzle_inf>& drz_inf_vec,
    const string& output_folder_,
    const MosaicBoxRegion& full_output_region)
    :filter(serial_in_order),
     output_folder(output_folder_),
     output_folder_created(false)
{


    for (size_t i=0; i!=drz_inf_vec.size(); ++i) {
        CCDs_drizzled[drz_inf_vec[i].CCD_id]=false;
    }
    assert(CCDs_drizzled.size() == drz_inf_vec.size());

//    boost::filesystem::create_directories(output_folder);
    raw_mosaic_vals.pix = PixelArray2d<float>(full_output_region.x_dim(),
                          full_output_region.y_dim(), 0.0);
    raw_mosaic_vals.initialize_mosaic_grid_to_specific_region(full_output_region);
    raw_mosaic_weights = raw_mosaic_vals;
}

void* raw_data_mosaic_filter::operator()(void* item)
{
    if (threaded_drizzle_debug) { cerr<<"In raw mosaic filter"; }
    if (!output_folder_created) {
        boost::filesystem::create_directories(output_folder);
        output_folder_created=true;
    }
    drizzle_token& d = *static_cast<drizzle_token*>(item);

    //drizzle
    d.img.initialize_mosaic_grid_to_specific_region(d.drizzle_inf.input_mosaic_region);
//    d.debiased_input.write_to_file(output_folder+
//            "frame_id_"+string_utils::itoa(d.frame_inf.corrected_frame_id)+
//            "_CCD_id_"+string_utils::itoa(d.frame_inf.ccd_id)+".fits" );

//    cerr<<"CCD_id: "<< d.frame_inf.ccd_id<<" cam_conf.get_mosaic_offset_for_ccd_id(ccd_id) "<<cam_conf.get_mosaic_offset_for_ccd_id(ccd_id)<<endl;
//    cerr<<"d.drizzle_inf.input_mosaic_region.low " << d.drizzle_inf.input_mosaic_region.low<<endl;

    MosaicPixelShift fixed_Positioning_shift = d.drizzle_inf.input_mosaic_region.low -
            raw_mosaic_vals.mosaic_grid.image_outline_.low;

    assert(fabs(d.img.mosaic_grid.pixel_width_ - 1.0) <error_margin) ;
    //This should currently always be true
    assert(fabs(d.img.CCD_grid.pixel_width_ - 1.0) < error_margin);

    PixelShift drizzle_input_shift = d.img.mosaic_grid.corresponding_pixel_shift(
                                         fixed_Positioning_shift);
//    +d.debiased_input.CCD_grid.corresponding_pixel_shift(d.atmospheric_translation_shift); //no atmos shift applied
    drizzle_input_shift.x=ceil(drizzle_input_shift.x);
    drizzle_input_shift.y=ceil(drizzle_input_shift.y);


//    cerr<<"\nCCD_ID:"<<d.frame_inf.ccd_id<<  " Shift:"<<drizzle_input_shift<<endl;


    drizzle::translate_and_drizzle_frame(d.img.pix, d.drizzle_input_weights.pix,
                                         raw_mosaic_vals.pix, raw_mosaic_weights.pix,
                                         drizzle_input_shift,
                                         1.0, 0.01
                                        );

//    cerr<<"Drizzled"<<endl;

    CCDs_drizzled[d.frame_inf.ccd_id]=true;
    bool all_ccds_drizzled=true;

    for (map<int,bool>::const_iterator it = CCDs_drizzled.begin(); it!=CCDs_drizzled.end();
            ++it) {
        if (it->second==false) { all_ccds_drizzled = false; }
    }

    if (all_ccds_drizzled) {
        //output
        raw_mosaic_vals.write_to_file(output_folder+"mos_frame_id_"+string_utils::itoa(
                                          d.frame_inf.corrected_frame_id)+".fits");
        raw_mosaic_vals.pix.assign(0.0f);
        raw_mosaic_weights.pix.assign(0.0f);
        for (map<int,bool>::iterator it = CCDs_drizzled.begin(); it!=CCDs_drizzled.end(); ++it) {
            it->second=false;
        }
    }
    if (threaded_drizzle_debug) { cerr<<"Done\n"; }
    return &d;
}

//=======================================================================================
output_summation_filter::output_summation_filter(const vector<CCD_specific_drizzle_inf>&
        drz_inf_vec, const bool thresholding_on_):
    filter(serial_in_order), thresholding_on(thresholding_on_)
{

    for (size_t i=0; i!=drz_inf_vec.size(); ++i) {
        index[drz_inf_vec[i].CCD_id] = i ;
        CCD_specific_drizzle_inf drz_inf = drz_inf_vec[i];
        MosaicImage<double> output_img;
        output_img.pix = PixelArray2d<double>(drz_inf.drizzled_output_PixelRange.x_dim(),
                                              drz_inf.drizzled_output_PixelRange.y_dim(), 0.0);
        output_img.initialize_mosaic_grid_to_specific_region(drz_inf.output_mosaic_region);

        summed_CCD_drizzle_val_imgs.push_back(output_img);
        summed_CCD_drizzle_weight_imgs.push_back(output_img);

        if (thresholding_on) {
            summed_threshed_vals.push_back(output_img);
//            summed_threshed_vals.back().add_keyword("THRESHED",string_utils::ftoa(drz_inf_vec[i].precalculated_threshold_level),"Datacount Level at which data was thresholded");
        }
    }

}

void output_summation_filter::reset()
{
    for (size_t i=0; i!=summed_CCD_drizzle_val_imgs.size(); ++i) {
        summed_CCD_drizzle_val_imgs[i].pix.assign(0.0);
        summed_CCD_drizzle_weight_imgs[i].pix.assign(0.0);
        if (thresholding_on) { summed_threshed_vals[i].pix.assign(0.0); }
    }
}

bool output_summation_filter::ccd_id_valid(const int ccd_id) const
{
//            if (ccd_id>=(int)index_vec.size() || index_vec[ccd_id]==-1) return false; else return true;
    return index.find(ccd_id)!=index.end();
}

size_t output_summation_filter::get_vec_Position_for_ccd_id(const int ccd_id) const
{
    assert(ccd_id_valid(ccd_id));
    return (index.find(ccd_id))->second;
}

void* output_summation_filter::operator()(void* item)
{
    if (threaded_drizzle_debug) { cerr<<"In output sum filter"<<endl; }
    drizzle_token& d = *static_cast<drizzle_token*>(item);

    size_t vec_pos = get_vec_Position_for_ccd_id(d.frame_inf.ccd_id);

    summed_CCD_drizzle_val_imgs[vec_pos].pix+=d.drizzle_result_vals.pix;
    summed_CCD_drizzle_weight_imgs[vec_pos].pix+=d.drizzle_result_weights.pix; //could split this into 2 parallel filters...
    if (thresholding_on) { summed_threshed_vals[vec_pos].pix+=d.thresholded_drizzle_result_vals.pix; }
    return &d;
}

CCDImage<double> output_summation_filter::drizzled_vals_for_CCD(const int ccd_id)
{
    if (!ccd_id_valid(ccd_id)) {
        throw logic_error("This ccd not stored by summation filter: "+string_utils::itoa(
                              ccd_id));
    }
    return summed_CCD_drizzle_val_imgs[get_vec_Position_for_ccd_id(ccd_id)];
}

CCDImage<double> output_summation_filter::drizzled_threshed_vals_for_CCD(const int ccd_id)
{
    if (!ccd_id_valid(ccd_id)) {
        throw logic_error("This ccd not stored by summation filter: "+string_utils::itoa(
                              ccd_id));
    }
    return summed_threshed_vals[get_vec_Position_for_ccd_id(ccd_id)];
}

CCDImage<double> output_summation_filter::drizzled_weights_for_CCD(const int ccd_id)
{
    if (!ccd_id_valid(ccd_id)) {
        throw logic_error("This ccd not stored by summation filter: "+string_utils::itoa(
                              ccd_id));
    }
    return summed_CCD_drizzle_weight_imgs[get_vec_Position_for_ccd_id(ccd_id)];
}


}//end namespace coela::multi_CCD_filters
}//end namespace
}//end namespace coela
