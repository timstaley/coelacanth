/*
 * File:   frame_info.h
 * Author: ts337
 *
 * Created on 03 April 2009, 16:45
 */

#ifndef COELA_FRAME_INFO_H
#define COELA_FRAME_INFO_H

#include "file_info.h"
#include "frame_registration.h"
#include "coela_core/src/cartesian_coords.h"
#include <iostream>
#include<cassert>


namespace coela {
//=========================================================================================
//Forward declaration:
struct DrizzleSettings;


//=========================================================================================
/// struct frame_info:
/// Holds information about a single CCD exposure.
/// Inherits frame location information via "file_info" struct.

struct FrameInfo: public file_info {
    FrameInfo(): quality_percentile_rank(0.0), bias_pedestal(-42.f) {}
    FrameInfo(const file_info& basic): file_info(basic), quality_percentile_rank(0.0),
        bias_pedestal(-42.f) {};

    static std::list<FrameInfo> get_frames_list(const string& relative_path,
            const string& filestem,
            const string& extension,
            int CCD_id=-1, const string& index_output_dir="");

    static std::vector<FrameInfo> get_frames_vec(const string& relative_path,
            const string& filestem,
            const string& extension,
            int CCD_id=-1, const string& index_output_dir="");

    static std::list<FrameInfo> load_frames_list(const string filename);
    static std::vector<FrameInfo> load_frames_vec(const string filename);

    static FrameInfo load_from_formatted_text(const string& formatted_text);
    static std::list<FrameInfo> load_frames_list(const vector<string>& list_text);

    static void write_list_to_file(const std::list<FrameInfo>& frames, string filename);
    static void write_list_to_file(const std::vector<FrameInfo>& frames, string filename);


    //----------------------------------------------------------------------------------/
    //Data members:
    std::vector< frame_registration::gs_lock<CCD_Position> > guide_star_estimates;
    double quality_percentile_rank;
    std::vector<CCD_Position> confirmed_cosmic_rays;
    size_t cube_FITS_z_index;
    float bias_pedestal;
    //----------------------------------------------------------------------------------/

    static void convert_directory_paths_to(const std::string& new_dir_path,
                                           vector<FrameInfo>&);

    static std::vector<FrameInfo> merge_frame_vecs(
        const std::vector<FrameInfo>& input1,  const std::vector<FrameInfo>& input2);
    static std::list<FrameInfo> split_last_guide_star_into_separate_list(
        std::list<FrameInfo>& frames1);



    static bool first_star_signal_predicate(const FrameInfo& first, const FrameInfo& second);
    static bool first_star_location_vector_length_predicate(const FrameInfo& first,
            const FrameInfo& second);

    static void set_list_CCD_id_to(std::list<FrameInfo>& frames, int CCD_id);

    static void add_percentile_rankings_based_on_first_star(std::vector<FrameInfo>& frms);

    static void set_bias_pedestal_estimates(std::vector<FrameInfo>& frms,
                                            double bias_pedestal_estimate);
    static double get_median_bias_pedestal_estimate(const std::vector<FrameInfo>& frms);

    static FrameInfo get_frame_with_id(const int corrected_frame_id,
                                       const std::vector<FrameInfo>&);
    static std::vector<FrameInfo>& shift_frame_id_numbers(int shift, std::vector<FrameInfo>&);

    static frame_registration::gs_lock<CCD_Position> mean_guide_star_estimate(
        const std::vector<FrameInfo>&,
        int gs_index=0);

    static std::vector<double> pull_quality_estimates_for_gs(
        const std::vector<FrameInfo>&,
        int gs_index=0);
    static void normalise_guide_star_quality_estimates(std::vector<FrameInfo> & frame_vec);

    static std::vector<FrameInfo> combine_guide_star_estimates_with_weights(
        const std::vector<FrameInfo> & frame_vec,
        const std::vector<double> & gs_weights);

    static void set_blank_guide_star_positions(std::vector<FrameInfo>&);

private:
    static std::list<FrameInfo> convert_file_list_to_frame_list(const std::list<file_info>&
            files); ///<Initialise the camera_id depending on the folder Position
};
std::ostream& operator<<(std::ostream& os, const FrameInfo& frm);

struct nth_star_quality_predicate {
    nth_star_quality_predicate(size_t star_index_):star_index(star_index_) {}
    bool operator()(const FrameInfo& first, const FrameInfo& second) {
        assert(star_index<first.guide_star_estimates.size());
        assert(star_index<second.guide_star_estimates.size());
        return (first.guide_star_estimates[star_index].signal >
                second.guide_star_estimates[star_index].signal);
    }
    const size_t star_index;
};

std::vector<FrameInfo> frames_in_intersection_by_corrected_frame_id(
    const std::vector<FrameInfo>& copy_from_list,
    const std::vector<FrameInfo>& match_to_list);
std::vector<int> pull_corrected_frame_ids(const std::vector<FrameInfo>&);
std::vector<int> intersecting_integers(const std::vector<int>& , const std::vector<int>&);

struct MultiFrame {
    static std::vector<MultiFrame> collate_frames_present_for_all_ccds(
        const std::vector< std::vector<FrameInfo> >& lists);

    static MultiFrame pull_multi_frame_with_id(int frame_id, const std::vector<MultiFrame>&);

    static void correlate_timestamps_and_correct_frame_IDs(vector< vector<FrameInfo> >&
            lists);

    static std::vector<FrameInfo> retrieve_frames_for_CCD_id(int CCD_id,
            const std::vector<MultiFrame>&);

    static std::vector<MultiFrame> get_multi_frame_vec(const std::vector<FrameInfo>&
            gs_frm_vec,
            const std::vector< std::vector<FrameInfo> >& additional_frame_lists,
            bool perform_timestamp_correlation_and_correction);


    MultiFrame(int corrected_frame_id):frame_id(corrected_frame_id),
        multi_frame_has_been_drizzled(false) {}
    int frame_id;
    std::vector<FrameInfo> synchronized_CCD_frames;
    bool multi_frame_has_been_drizzled;

private:
    static const size_t correlation_size=50;
    static const size_t first_list_corr_offset=correlation_size;
};



struct cross_detector_GS_ests {
    vector<MosaicPosition> gs_estimates;
    int frame_id;
};

//======================================================================================================
}//end of namespace coela
#endif  /* _FRAME_INFO_H */

