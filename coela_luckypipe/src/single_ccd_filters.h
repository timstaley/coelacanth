/*
 * File:   single_ccd_filters.h
 * Author: ts337
 *
 * Created on 22 February 2011, 17:53
 */

#ifndef COELA_SINGLE_CCD_FILTERS_H
#define COELA_SINGLE_CCD_FILTERS_H

#include "coela_utility/src/file_buffer.h"
#include "coela_utility/src/histogram_container.h"
#include "coela_core/src/mosaic_image.h"
#include "frame_info.h"
#include "frame_registration.h"

#include <tbb/pipeline.h>
#include <boost/filesystem.hpp>


namespace coela {
namespace single_CCD_filters {
//=========================================================================================================
const size_t lucky_n_tokens = 8;

struct frame_cleanup_token {
    frame_cleanup_token() {  clear(); }

    //These get reset by "clear"
    FrameInfo frame_inf;
//    float original_bias_pedestal; //for outputting a log of variation.
    bool frame_is_good; //set to false if guide star no good due to cosmic ray proximity.
    bool cosmic_rays_confirmed_due_to_edge_case;
    bool cosmic_rays_rejected_temporally;
    bool cosmic_rays_rejected_spatially;
    vector<PixelIndex> cosmic_ray_candidates;


    //These are placeholders to be reassigned.
    HistogramContainer14bit raw_data_zeroing_histogram_space;
    FileBuffer buffer;
    FitsHeader fht;
    MosaicImage<float> img;

    //Gets set during initialisation
    size_t token_buffer_id; //Can be used to index into corresponding arrays of temporary data structs used only by a particular filter

    virtual void clear();

//  virtual ~frame_token(){} //Don't need this since we're still using static allocation in vectors via templates... but be aware if usage changed.

    template <class T>
    static vector<T> initialize_token_vec(const size_t n_tokens);
};

//----------------------------------------------------------------
//Template definition needs including with header....
template <class T>
vector<T> frame_cleanup_token::initialize_token_vec(const size_t n_tokens)
{
    vector< T > v(n_tokens);
    for (size_t i=0; i!=v.size(); ++i) {
        v[i].token_buffer_id=i;
    }
    return v;
}
//=========================================================================================================

class frame_predicate {
public:
    virtual bool operator()(const frame_cleanup_token& b) const=0;
    virtual ~frame_predicate() {}
};

class dummy_frame_predicate: public frame_predicate {
    bool operator()(const frame_cleanup_token&)const {return true;}
};
static dummy_frame_predicate dummy_predicate;

//=========================================================================================================
//filters...

//===============================================================================================
//Input and output filters.

//Initial buffering filter:
template < class T = frame_cleanup_token>
class Sequential_File_Buffer_Filter: public tbb::filter {
public:
    Sequential_File_Buffer_Filter(const vector<FrameInfo>& frame_info_vec,
                                  size_t n_frames_to_load=0);
    size_t bytes_loaded() {return bytes_off_disk;}
    size_t n_frames_set_to_process() const {return frm_inf_vec.size();}
    void* operator()(void*);
    void reset();
private:
    vector<FrameInfo> frm_inf_vec;
    vector<FrameInfo>::size_type current_posn_index;
    size_t next_buffer, bytes_off_disk, counter;
    vector<frame_cleanup_token*> token_buffer_ptrs;
    vector<T> token_buffers;

};
//---------------------------------------------------------------------------------------------------------------

//NB decompresses to an image<unsigned int>
///in case of uncompressed files, this simply loads the fits file from a buffer, which should happen very quickly
class Decompress_Filter: public tbb::filter {
public:
    Decompress_Filter():filter(parallel) {}
    void* operator()(void*);
};

//---------------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------------

//For outputting a small subset of files, based on a predicate:
class Frame_Output_Filter: public tbb::filter {
public:
    Frame_Output_Filter(const std::string& output_folder_,
                        const frame_predicate& predicate=dummy_predicate):
        tbb::filter(serial_in_order),
        output_folder(output_folder_),
        n_frames_output(0),
        frame_predicate_ref(predicate) {}
    void* operator()(void* item);
    size_t number_of_frames_output() const {return n_frames_output;}
private:
    const std::string output_folder;
    size_t n_frames_output;
    const frame_predicate & frame_predicate_ref;
};

//---------------------------------------------------------------------------------------------------------------

//For outputting every file:
class Frame_Sequential_Write_Filter: public tbb::filter {
public:
    Frame_Sequential_Write_Filter(const string& output_dir_):filter(serial_in_order),
        output_dir(output_dir_) {}
    void* operator()(void*);
private:
    const string output_dir;
};
//---------------------------------------------------------------------------------------------------------------

//To collect a new vector of frame_info with the guide star registration / cosmic ray info
class Frame_Info_Collection_Filter: public tbb::filter {
public:
    Frame_Info_Collection_Filter(size_t list_size): tbb::filter(
            serial_in_order) {frames.reserve(list_size);}
    void* operator()(void* item);
    vector<FrameInfo> collected_frames() {return frames;}
private:
    vector<FrameInfo> frames;
};

//---------------------------------------------------------------------------------------------------------------

//To ensure the circular buffer is used correctly:
class Serial_Decommission_Filter: public tbb::filter {
public:
    Serial_Decommission_Filter(): tbb::filter(serial_in_order) {}
    void* operator()(void* item);
};
//---------------------------------------------------------------------------------------------------------------

//=========================================================================================================

//General purpose filters that can be inserted pretty much anywhere along the pipeline.
//---------------------------------------------------------------------------------------------------------------
class File_Timestamp_Imprinting_Filter: public tbb::filter {
public:
    File_Timestamp_Imprinting_Filter(const FrameInfo& first_frame,
                                     const double milliseconds_per_frame_);
    void* operator()(void*);
private:
    time_t init_time;
    int init_frame_number;
    const double milliseconds_per_frame;

};

//---------------------------------------------------------------------------------------------------------------
class Timestamp_Load_Filter: public tbb::filter {
public:
    Timestamp_Load_Filter():filter(parallel) {}
    void* operator()(void*);
};

//---------------------------------------------------------------------------------------------------------------

class Frame_Count_Display_Filter: public tbb::filter {
public:
    Frame_Count_Display_Filter(bool display_frame_count_,
                               const frame_predicate& predicate=dummy_predicate,
                               const size_t flag_throw_limit=0):
        filter(serial_in_order), counter(0),
        display_frame_count(display_frame_count_),
        frame_predicate_ref(predicate),
        throw_limit(flag_throw_limit)
    {}
    void* operator()(void*);
    long n_counted() {return counter;}
    void reset() {counter=0;}
private:
    size_t counter;
    const bool display_frame_count;
    const frame_predicate& frame_predicate_ref;
    const size_t throw_limit;
};

//=========================================================================================================

//---------------------------------------------------------------------------------------------------------------

class FrameCropFilter: public tbb::filter {
public:
    FrameCropFilter(const PixelRange& crop_region):filter(parallel),crop_rgn(crop_region) {}
    void* operator()(void* item);
private:
    const PixelRange crop_rgn;
};

//---------------------------------------------------------------------------------------------------------------

class BiasDriftTracker: public tbb::filter {
public:
    BiasDriftTracker(const PixelRange& hist_region):
        filter(parallel), hist_region_box_(hist_region) {}
    void* operator()(void* item);   //NB ASSUMES PIXELS ARE STILL INTEGER VALUED
private:
    const PixelRange hist_region_box_;//, histogram_rgn;
};

//For simulated data which does not get debiased:
//class Float_To_Double_Filter: public tbb::filter {
//public:
//    Float_To_Double_Filter():filter(parallel){}
//    void* operator()( void* item ); //NB ASSUMES PIXELS ARE STILL INTEGER VALUED
//};

//---------------------------------------------------------------------------------------------------------------

//=========================================================================================================

//Filters that should be used on debiased data (i.e. they work on the image<float>)

//---------------------------------------------------------------------------------------------------------------
class Pixel_Time_Series_Record: public tbb::filter {
public:
    typedef std::pair< int, double> pixel_event; //frame id, pixel value
    Pixel_Time_Series_Record(const vector<PixelIndex>& pixels_to_watch,
                             const size_t n_frames_estimate);
    void* operator()(void* item);

    vector< std::pair< PixelIndex, vector<pixel_event> > > get_time_series() const;

private:
    vector< PixelIndex> pixel_indices;
    vector < vector< pixel_event > > pixval_series;
};
//---------------------------------------------------------------------------------------------------------------
class Col_Histogram_Gather: public tbb::filter {
public:
    Col_Histogram_Gather(const PixelRange& image_layout,
                         int y_range_low, int y_range_high);

    void* operator()(void*);
    void set_hist_min_value(const int);
    const vector<HistogramContainer14bit>& get_hists()const {return col_histograms;}
private:
    vector<HistogramContainer14bit> col_histograms;
    const int y_lo, y_hi;
};
//---------------------------------------------------------------------------------------------------------------
class RowHistogram_Gather: public tbb::filter {
public:
    RowHistogram_Gather(const PixelRange& image_layout,
                        int x_range_low);

    void* operator()(void*);
    void set_hist_min_value(const int);
    const vector<HistogramContainer14bit>& get_hists()const {return row_histograms;}
private:
    vector<HistogramContainer14bit> row_histograms;
    const int x_lo_;
};
//---------------------------------------------------------------------------------------------------------------
///Simply subtracts the value stored in frame_info::bias_pedestal_
class UniformDebias: public tbb::filter {
public:
    UniformDebias():filter(parallel) {}
    /*override*/void* operator()(void* item);
};
//---------------------------------------------------------------------------------------------------------------
class BiasFrameSubtractor: public tbb::filter {
public:
    BiasFrameSubtractor():filter(parallel) {}

//            BiasFrameSubtractor_Filter(const image<double>& image_to_subtract):filter(parallel), bias_frame(image_to_subtract){}
    /*override*/void* operator()(void* item);
    void set_bias_frame(const CCDImage<float>& bias_frame_input) {bias_frame=bias_frame_input;}
private:
    CCDImage<float> bias_frame;
};
//---------------------------------------------------------------------------------------------------------------
//to do: templatize (worth the effort?)
class Histogram14BitBuildFilter: public tbb::filter {
public:
    Histogram14BitBuildFilter(const PixelRange& hist_rgn, int minimum_hist_value):
        filter(serial_out_of_order), hist_box_(hist_rgn) {
        histogram_.set_min_value(minimum_hist_value);
    }

    void* operator()(void* item);
    HistogramContainer14bit histogram() const {return histogram_;}
private:
    HistogramContainer14bit histogram_;
    PixelRange hist_box_;
};
//---------------------------------------------------------------------------------------------------------------
//class Histogram10BitBuildFilter: public tbb::filter{
//public:
//    Histogram10BitBuildFilter(const PixelRange& hist_rgn, int minimum_hist_value):
//        filter(serial_out_of_order), hist_box_(hist_rgn){
//            histogram_.set_min_value(minimum_hist_value);
//        }
//
//    void* operator()( void* item );
//    HistogramContainer10bit histogram() const{return histogram_;}
//private:
//    HistogramContainer10bit histogram_;
//    PixelRange hist_box_;
//};
//---------------------------------------------------------------------------------------------------------------
struct registration_data_struct {
    vector<CCD_BoxRegion> gs_regions;
    psf_models::reference_psf ref_psf;
    double analysis_resample_factor, convolution_threshold;
};


class CrossCorrelator: public tbb::filter {
public:
    CrossCorrelator(
        const vector<CCD_BoxRegion>& gs_regions_init,
        const psf_models::reference_psf& kernel_init,
        const double input_resample_factor,
        const double convolution_threshold_factor,
        const double ray_proximity_limit_in_pixels,
        const bool use_parabola_fit_to_xcorr);

    void* operator()(void* item);
private:
    registration_data_struct reg_data_copies_[lucky_n_tokens];
    const bool use_parabola_fitting_;
    const double ray_proximity_limit_;
};
//---------------------------------------------------------------------------------------------------------------
class Dynamic_Row_Debias: public tbb::filter {
public:
    Dynamic_Row_Debias():filter(
            parallel) {} //NB CCD_id merely supplied for debugging purposes,
    //and to emphasise that this should not be used with a multi-CCD token stream!

    /*override*/void* operator()(void* item);
};
//---------------------------------------------------------------------------------------------------------------
class Frame_Summation: public tbb::filter {
public:
    Frame_Summation(const PixelRange& image_layout, const CCD_BoxRegion& image_region):
        tbb::filter(serial_out_of_order),
        n_frames_summed(0) {
        sum.pix = PixelArray2d<double>(image_layout.x_dim(), image_layout.y_dim(), 0.0);
        sum.initialize_CCD_grid_to_specific_region(image_region);
    }
    void* operator()(void* item);
    CCDImage<double> Sum() {return sum;}
    size_t N_Frames() {return n_frames_summed;}
    CCDImage<double> Avg();
private:
    CCDImage<double> sum;
    size_t n_frames_summed;
};
//---------------------------------------------------------------------------------------------------------------
class Region_Centroid: public tbb::filter {
public:
    Region_Centroid(const PixelRange& GS_box_, size_t list_size):
        tbb::filter(serial_in_order), GS_box(GS_box_) {frames.reserve(list_size);}
    void* operator()(void* item);
    vector<FrameInfo> collected_frames() {return frames;}
private:
    vector<FrameInfo> frames;
    const PixelRange GS_box;
};


//---------------------------------------------------------------------------------------------------------------

class Cosmic_Ray_Candidate_Detection: public tbb::filter {
public:
    Cosmic_Ray_Candidate_Detection(double raw_count_threshold):
        tbb::filter(parallel),
        threshold(raw_count_threshold),
        n_frames_with_candidates(0)
    {}
    void* operator()(void* item);
    size_t number_of_frames_with_pixel_above_threshold() const {return n_frames_with_candidates;}
private:
    const double threshold;
    size_t n_frames_with_candidates;

};

//---------------------------------------------------------------------------------------------------------------

class Cosmic_Ray_Confirmation: public tbb::filter {
public:
    Cosmic_Ray_Confirmation(const double frame_to_frame_reduction_factor,
                            const double pixel_spaced_reduction_factor):
        tbb::filter(serial_in_order), first_frame(true),
        temporal_reduction_factor(frame_to_frame_reduction_factor),
        spatial_reduction_factor(pixel_spaced_reduction_factor),
        n_frames_with_confirmed_rays(0)
    {}
    void* operator()(void* item);
    size_t number_of_frames_with_confirmed_rays() const {return n_frames_with_confirmed_rays;}
private:
    bool first_frame;
    FrameInfo prev_frame_inf;
    CCDImage<float> prev_bmp;
    const double temporal_reduction_factor, spatial_reduction_factor;
    size_t n_frames_with_confirmed_rays;
};

//---------------------------------------------------------------------------------------------------------------

class Cosmic_Ray_Diagnosis_Output: public tbb::filter {
public:
    Cosmic_Ray_Diagnosis_Output(const std::string& CCD_output_folder,
                                const std::string& cr_diag_subfolder):
        tbb::filter(serial_in_order), output_folder(CCD_output_folder+cr_diag_subfolder) {
        boost::filesystem::create_directories(output_folder+"/confirmed_clearly");
        boost::filesystem::create_directories(output_folder+"/confirmed_due_to_edge_case");
        boost::filesystem::create_directories(output_folder+"/rejected_spatially");
        boost::filesystem::create_directories(output_folder+"/rejected_temporally");
        boost::filesystem::create_directories(output_folder+"/rejected_both");
    }

    void* operator()(void* item);

    size_t number_of_frames_output() const {return n_frames_output;}
    vector<double> get_candidate_ray_pixel_peak_values() const { return confirmed_peak_vals; }
private:
    const std::string output_folder;
    size_t n_frames_output;
    vector<double> confirmed_peak_vals;
};
//---------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------
//=========================================================================================================







//=========================================================================================================
} //end namespace coela::single_CCD_filters
}//end namespace coela

#endif  /* SINGLE_CCD_FILTERS_H */

