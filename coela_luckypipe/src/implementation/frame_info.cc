#include "../frame_info.h"

#include <vector>
#include "coela_utility/src/string_utils.h"
#include "coela_utility/src/simple_serialization.h"
#include "coela_utility/src/microstats.h"


#include "coela_core/src/pixel_array_header.h"
//#include "../fits_header.h"
//#include "../list_analysis_code/list_analysis.h"
#include <boost/filesystem.hpp>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <coela_core/src/cartesian_coords.h>
using namespace std;
namespace coela {
using frame_registration::GuideStarLock;
using namespace string_utils;

void FrameInfo::write_list_to_file(const list<FrameInfo>& frames, string filename)
{
    std::ofstream outfile(filename.c_str(), ios::binary);
    try {
        if (outfile.is_open()) {
            outfile <<"#Lucky Frames Listing Version 1.0\n";
            outfile <<"# Number  ;   "
                    <<"Primary GS [; 2nd GS  ...] ; percentile_rank "
                    <<"bias_pedestal "
                    <<"file_path ccd_id "
                    <<"byte_offset byte_size "
                    <<"is_cube cube_index "
                    <<"header_frm_id corrected_frm_id timestamp; cosmic_rays \n" ;
            size_t counter=1;
            outfile.precision(7);
            for (list<FrameInfo>::const_iterator it(frames.begin()); it!=frames.end(); ++it) {
                outfile << counter++ <<";\t";
                outfile << *it <<"\n";
            }
            outfile.close();
        } else { throw runtime_error(filename); }
    } catch (runtime_error fname)
    {   cerr <<"Could not write to file: " << fname.what()<<", exiting"<<endl; exit(0); }
}

std::ostream& operator<<(std::ostream& os, const FrameInfo& frm)
{
    if (frm.guide_star_estimates.empty()) {
        os<<"0";
    } else {
        for (size_t i=0; i<frm.guide_star_estimates.size(); ++i) {
            os << frm.guide_star_estimates[i]<<" ";
        }
    }
    os<<";";
    os<< frm.quality_percentile_rank<<" "
      <<frm.bias_pedestal<<" "
      <<frm.file_path <<" "<<frm.ccd_id <<" "
      <<frm.header_byte_offset <<" "<<frm.byte_size<<" "
      <<frm.file_is_cube_FITS<<" "
      <<(frm.file_is_cube_FITS ? frm.cube_FITS_z_index : 0)<<" "
      <<frm.subfile_number<<" "
      <<frm.header_frame_id<<" "
      <<frm.corrected_frame_id<<" "
      <<frm.header_timestamp;
    os<<";";
    if (frm.confirmed_cosmic_rays.empty()) { os<<"0"; }
    else {
        for (size_t i=0; i!=frm.confirmed_cosmic_rays.size(); ++i) {
            os<<frm.confirmed_cosmic_rays[i]<<" ";
        }
    }

    return os;
}

FrameInfo FrameInfo::load_from_formatted_text(const string& formatted_text)
{
    FrameInfo frm;
    stringstream ss;
    vector<string> segments = tokenize(formatted_text,";");
    assert(segments.size()==3);

    if (segments.front()!=string("0")) {
        ss.str(segments.front());
        GuideStarLock<CcdPosition> gl;
        while (ss>>gl) { frm.guide_star_estimates.push_back(gl); }
    }
    ss.clear();

    ss.str(segments[1]);
    ss>>frm.quality_percentile_rank
      >>frm.bias_pedestal
      >>frm.file_path >>frm.ccd_id
      >>frm.header_byte_offset >> frm.byte_size
      >>frm.file_is_cube_FITS
      >>frm.cube_FITS_z_index
      >>frm.subfile_number
      >>frm.header_frame_id>>frm.corrected_frame_id>>frm.header_timestamp;
    if (frm.byte_size==0) { frm.derived_lcc_filename = string_utils::pull_filename(frm.file_path); }
    else { frm.derived_lcc_filename = string_utils::pull_filestem(frm.file_path)+"_"+itoa(frm.ccd_id)+"_" +itoa(frm.subfile_number)+".fits"; }

    ss.clear();
    if (segments.back() !=string("0")) {
        ss.str(segments.back());
        CcdPosition ray;
        while (ss>>ray) { frm.confirmed_cosmic_rays.push_back(ray); }
    }
    return frm;
}



std::list<FrameInfo> FrameInfo::load_frames_list(const vector<string>& list_text)
{
    list<FrameInfo> frames;
    if (list_text.front().find("#Lucky Frames Listing Version 1.0")!=0) { throw runtime_error("List text does not appear to be correct version frames listing (1.0)"); }
    size_t i=0;
    while (list_text[i][0]=='#') { ++i; }
    for (; i!=list_text.size(); ++i) {
        string::size_type comma_pos = list_text[i].find(';');
        if (comma_pos==string::npos || comma_pos==list_text[i].size()-1) { throw runtime_error("Line " +itoa(i) +"appears incorrectly formatted"); }
        string frm_text = list_text[i].substr(comma_pos
                                              +1); //snip off the counter number at the beginning
        frames.push_back(FrameInfo::load_from_formatted_text(frm_text));
    }
    return frames;
}





void FrameInfo::convert_directory_paths_to(const std::string& new_dir_path,
        vector<FrameInfo>& frames)
{
    string appended_path = new_dir_path;
    if (!new_dir_path.empty() && new_dir_path[new_dir_path.size()-1]!='/') { appended_path+="/"; }
    for (size_t i=0; i!=frames.size(); ++i) {
        frames[i].file_path = appended_path + string_utils::pull_filename(frames[i].file_path);
    }
}


std::vector<FrameInfo> FrameInfo::get_frames_vec(const string& relative_path,
        const string& filestem, const string& extension, int CCD_id,
        const string& index_output_dir)
{
    list<FrameInfo> loaded_list = get_frames_list(
                                      relative_path,
                                      filestem,
                                      extension,
                                      CCD_id,
                                      index_output_dir
                                  );
    return vector<FrameInfo> (loaded_list.begin(), loaded_list.end());
}

list<FrameInfo> FrameInfo::get_frames_list(const string& relative_path,
        const string& filestem, const string& extension, int CCD_id,
        const string& index_output_dir)
{
    if (!boost::filesystem::is_directory(relative_path)) {throw runtime_error("frame_info::get_frames_list: Dir: \""+relative_path+"\" does not exist"); }
    list<FileInfo> files;
    if (extension.find("lcz")!=string::npos) {
        files = lcz_utils::get_lcz_subfile_list(relative_path, index_output_dir);
        lcz_utils::sanitise_list(files);
    } else { files=FileInfo::get_image_file_list(relative_path, filestem, CCD_id, extension); }
    list<FrameInfo> frames(convert_file_list_to_frame_list(files));
    frames.sort(FrameInfo::frame_id_predicate);
    return frames;
}




void FrameInfo::write_list_to_file(const std::vector<FrameInfo>& frames, string filename)
{
    list<FrameInfo> temp_list(frames.begin(), frames.end());
    FrameInfo::write_list_to_file(temp_list, filename);
}


list<FrameInfo> FrameInfo::load_frames_list(const string fname)
{
    std::ifstream framesfile(fname.c_str(), ios::binary);
    try {
        vector<string> list_text = simple_serialization::convert_istream_to_string_vec(
                                       framesfile);
        framesfile.close();
        return load_frames_list(list_text);
    } catch (runtime_error e) {
        cerr <<"Could not open frames list at "<< fname<<": " << e.what()<<", exiting"<<endl;
        exit(0);
    }
}
std::vector<FrameInfo> FrameInfo::load_frames_vec(const string filename)
{
    list<FrameInfo> l1=load_frames_list(filename);
    return vector<FrameInfo>(l1.begin(),l1.end());
}




std::vector<FrameInfo> FrameInfo::merge_frame_vecs(
    const std::vector<FrameInfo>& input1,  const std::vector<FrameInfo>& input2)
{
//        if (frames1.size()!=frames2.size() ) {throw runtime_error("\nCannot merge lists of different sizes\n");}
    vector<FrameInfo> frames1(input1), frames2(input2);

    sort(frames1.begin(), frames1.end(), FileInfo::subfile_number_predicate);
    sort(frames2.begin(), frames2.end(), FileInfo::subfile_number_predicate);
//        frames1.sort(FileInfo::subfile_number_predicate);
//        frames2.sort(FileInfo::subfile_number_predicate);
//        size_t counter=0;
//        cout <<"\r"<<counter<<"...         ";  cout.flush();
    for (vector<FrameInfo>::iterator it1(frames1.begin()), it2(frames2.begin());
            it1!=frames1.end()&&it2!=frames2.end();
            ++it1, ++it2) {
        FrameInfo& frame1(*it1), frame2(*it2);
        if (frame1.file_path!=frame2.file_path) { throw runtime_error("List merge- filename mismatch"); }
        for (size_t nstar = 0; nstar!=frame2.guide_star_estimates.size(); ++nstar) {
            frame1.guide_star_estimates.push_back(frame2.guide_star_estimates[nstar]);
        }
    }
    return frames1;
}

std::list<FrameInfo> FrameInfo::split_last_guide_star_into_separate_list(
    std::list<FrameInfo>& frames1)
{
    list<FrameInfo> split_list;
    for (list<FrameInfo>::iterator it(frames1.begin()); it!=frames1.end(); ++it) {
        FrameInfo temp_info(*it);
        temp_info.guide_star_estimates.clear();
        temp_info.guide_star_estimates.push_back(it->guide_star_estimates.back());
        split_list.push_back(temp_info);
        it->guide_star_estimates.pop_back();
    }
    return split_list;
}



list<FrameInfo> FrameInfo::convert_file_list_to_frame_list(const list<FileInfo>& files)
{
    list<FrameInfo> frames;
    for (list<FileInfo>::const_iterator it(files.begin()); it!=files.end(); ++it) {
        FrameInfo current_frame(*it);
        current_frame.bias_pedestal=0.0;

        if (it->file_is_cube_FITS) {
            int frame_number, corrected_frame_id;
            if (frames.empty()) {
                frame_number=0;
                corrected_frame_id=0;
            } else {
                corrected_frame_id = frames.back().corrected_frame_id;
                frame_number = frames.back().subfile_number;
            }
            FitsHeader temp_fht(it->file_path, it->header_byte_offset);
            PixelCubeHeader cube_hdr(temp_fht);
            for (size_t slice_index=1; slice_index<=cube_hdr.z_dim(); ++slice_index) {
                current_frame.cube_FITS_z_index = slice_index;
                current_frame.subfile_number=++frame_number;
                current_frame.corrected_frame_id=++corrected_frame_id;
                frames.push_back(current_frame);
            }
        } else {
            frames.push_back(current_frame);
        }
    }
    return frames;
}


bool FrameInfo::first_star_signal_predicate(const FrameInfo& first,
        const FrameInfo& second)
{
    assert(first.guide_star_estimates.size()>=1);
    assert(second.guide_star_estimates.size()>=1);
    return (first.guide_star_estimates.front().signal >
            second.guide_star_estimates.front().signal);
}



bool FrameInfo::first_star_location_vector_length_predicate(const FrameInfo& first,
        const FrameInfo& second)
{
    CcdPosition origin(0.0, 0.0);
    return (coord_distance_squared(origin, (first.guide_star_estimates.front().Position))
            < coord_distance_squared(origin, second.guide_star_estimates.front().Position));
}

//        bool frame_info::frame_order_predicate(const frame_info& first, const frame_info& second){
//            return (first.frame_number < second.frame_number);
//        }

FrameInfo FrameInfo::get_frame_with_id(const int corrected_frame_id,
                                       const std::vector<FrameInfo>& frms)
{
    for (vector<FrameInfo>::const_iterator it =frms.begin(); it!=frms.end(); ++it) {
        if (it->corrected_frame_id == corrected_frame_id) { return *it; }
    }
    throw runtime_error("Corrected frame id not found in list: "+ string_utils::itoa(
                            corrected_frame_id));
}

void FrameInfo::set_list_CCD_id_to(std::list<FrameInfo>& frames, int CCD_id)
{
    for (list<FrameInfo>::iterator it=frames.begin(); it!=frames.end(); ++it) {
        it->ccd_id=CCD_id;
    }
}


void FrameInfo::add_percentile_rankings_based_on_first_star(std::vector<FrameInfo>& frms)
{
    if (frms.front().guide_star_estimates.size()<1) {
        throw runtime_error("Error: Trying to sort list by GS quality, without any guide stars registered");
    }
    sort(frms.begin(), frms.end(), FrameInfo::first_star_signal_predicate);
    for (size_t i = 0; i!=frms.size(); ++i) {
        frms[i].quality_percentile_rank  = (double)(frms.size() - i) /   frms.size();
    }
    //Rearrange list back into chronological order so we can perform correlation on multiple lists
    sort(frms.begin(), frms.end(), FileInfo::frame_id_predicate);
}

void FrameInfo::set_bias_pedestal_estimates(std::vector<FrameInfo>& frms,
        double bias_pedestal_estimate)
{
    for (size_t i = 0; i!=frms.size(); ++i) {
        frms[i].bias_pedestal = bias_pedestal_estimate;
    }
}
double FrameInfo::get_median_bias_pedestal_estimate(const std::vector<FrameInfo>& frms)
{
    vector<float> bias_pedestal_ests;
    for (size_t i = 0; i!=frms.size(); ++i) {
        bias_pedestal_ests.push_back(frms[i].bias_pedestal);
    }
    return vector_median(bias_pedestal_ests);
}

vector<FrameInfo>& FrameInfo::shift_frame_id_numbers(int shift,
        std::vector<FrameInfo>& frms)
{
    for (vector<FrameInfo>::iterator it=frms.begin(); it!=frms.end(); ++it) {
        it->corrected_frame_id+=shift;
    }
    return frms;
}

GuideStarLock<CcdPosition> FrameInfo::mean_guide_star_estimate(const std::vector<FrameInfo>&
        frm_vec,
        int gs_index)
{
    CcdPixelShift sum_shift(0.0,0.0);
    double sum_qual=0.0;
    for (size_t frm=0; frm!=frm_vec.size(); ++frm) {
        sum_shift += (frm_vec[frm].guide_star_estimates[gs_index].Position -
                      CcdPosition::origin);
        sum_qual += frm_vec[frm].guide_star_estimates[gs_index].signal;
    }
    sum_shift/=frm_vec.size();
    sum_qual/=frm_vec.size();
    return GuideStarLock<CcdPosition>(sum_shift+CcdPosition::origin, sum_qual);
}

std::vector<double> FrameInfo::pull_quality_estimates_for_gs(
    const std::vector<FrameInfo>& frms,
    int gs_index)
{
    vector<double> quals; quals.reserve(frms.size());
    for (size_t i=0; i!=frms.size(); ++i) {
        quals.push_back(frms[i].guide_star_estimates[gs_index].signal);
    }
    return quals;
}

void FrameInfo::normalise_guide_star_quality_estimates(std::vector<FrameInfo> & frame_vec)
{
    assert(!frame_vec.empty());

    vector< double > gs_medians;
    for (size_t gs=0; gs!=frame_vec.front().guide_star_estimates.size(); ++gs) {
        vector<double> gs_quals = pull_quality_estimates_for_gs(frame_vec, gs);
        gs_medians.push_back(vector_median(gs_quals));
    }
    assert(gs_medians.size() == frame_vec.front().guide_star_estimates.size());
    for (size_t i=0; i!=frame_vec.size(); ++i) {

        for (size_t gs=0; gs!=frame_vec[i].guide_star_estimates.size(); ++gs) {
            frame_vec[i].guide_star_estimates[gs].signal /= gs_medians[gs];
        }
    }

}

std::vector<FrameInfo> FrameInfo::combine_guide_star_estimates_with_weights(
    const std::vector<FrameInfo> & input_frms,
    const std::vector<double> & gs_weights)
{
    assert(!input_frms.empty());
    assert(input_frms.front().guide_star_estimates.size() == gs_weights.size());

    double weight_sum=0.0;
    for (size_t i=0; i!=gs_weights.size(); ++i) {
        weight_sum+=gs_weights[i];
    }

    vector<FrameInfo> output(input_frms);

    normalise_guide_star_quality_estimates(output);

    for (size_t i=0; i!=output.size(); ++i) {
        FrameInfo& frm = output[i];
        double sum_qual=0.0;
        CcdPixelShift sum_position(0,0);
        for (size_t star_num=0; star_num!=frm.guide_star_estimates.size();
                ++star_num) {
            sum_qual +=
                frm.guide_star_estimates[star_num].signal *
                gs_weights[star_num];
            sum_position +=
                (frm.guide_star_estimates[star_num].Position - CcdPosition::origin) *
                gs_weights[star_num];
        }
        sum_qual/=weight_sum;
        sum_position/=weight_sum;

        frm.guide_star_estimates.clear();
        frm.guide_star_estimates.push_back(
            GuideStarLock<CcdPosition>(sum_position +CcdPosition::origin, sum_qual));
    }
    return output;
}


void FrameInfo::set_blank_guide_star_positions(std::vector<FrameInfo>& frms)
{
    for (size_t i=0; i!=frms.size(); ++i) {
        frms[i].guide_star_estimates.clear();
        frms[i].guide_star_estimates.push_back(GuideStarLock<CcdPosition>(CcdPosition(0,0)));
    }

}

std::vector<int> pull_corrected_frame_ids(const std::vector<FrameInfo>& v)
{
    vector<int> ids; ids.reserve(v.size());
    for (size_t i=0; i!=v.size(); i++) {
        ids.push_back(v[i].corrected_frame_id);
    }
    return ids;
}

std::vector<FrameInfo> frames_in_intersection_by_corrected_frame_id(
    const std::vector<FrameInfo>& copy_from_list, const std::vector<FrameInfo>& match_to_list)
{
    vector<int> v1, v2;
    v1=pull_corrected_frame_ids(copy_from_list);
    v2=pull_corrected_frame_ids(match_to_list);

    vector<int> intersection_ids = intersecting_integers(v1,v2);

    vector<FrameInfo>
    intersection_frames; //intersection is now a subset of the ids from the "copy_from" list
    intersection_frames.reserve(intersection_ids.size());
    size_t intersection_index=0, copy_from_index=0;
    while (intersection_index!=intersection_ids.size()
            &&copy_from_index!=copy_from_list.size()) {
        if (intersection_ids[intersection_index]==
                copy_from_list[copy_from_index].corrected_frame_id) {
            intersection_frames.push_back(copy_from_list[copy_from_index]);
            intersection_index++;
        }
        copy_from_index++;
    }
    assert(intersection_frames.size()==intersection_ids.size());
    return intersection_frames;

}

std::vector<int> intersecting_integers(const std::vector<int>& v1,
                                       const std::vector<int>& v2)
{
    vector<int> intersection_ints; intersection_ints.reserve(min(v1.size(),v2.size()));
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
                     back_inserter(intersection_ints));
    return intersection_ints;
}

vector<clock_t> pull_frame_timestamps(const vector<FrameInfo>& frames)
{
    vector<clock_t> timestamps;
    timestamps.reserve(frames.size());
    for (size_t i=0; i!= frames.size(); ++i) {
        timestamps.push_back(frames[i].header_timestamp);
    }
    return timestamps;
}

void MultiFrame::correlate_timestamps_and_correct_frame_IDs(vector< vector<FrameInfo> >&
        lists)
{
    cout <<"Matching timestamps: "<<endl;
    if (lists.size()==1) { return; }

    assert(lists.front().size()>first_list_corr_offset+correlation_size);

    vector<FrameInfo> offset_frames_from_first_list(lists.front().begin()
            +first_list_corr_offset,
            lists.front().begin()+first_list_corr_offset +correlation_size) ;
    vector<clock_t> offset_timestamps_from_first_list = pull_frame_timestamps(
                offset_frames_from_first_list);

    for (size_t i=0; i!=lists.size(); ++i) {
        vector<FrameInfo> & current_list=lists[i];
        cout <<"First timestamp for this list: " << current_list.front().header_timestamp<<endl;
        double min_fit=1e10; //very large, should find a better min than this!
        int best_matching_shift=0;
        assert(current_list.size()>first_list_corr_offset+correlation_size*3);
        for (size_t curr_offset=0; curr_offset!=first_list_corr_offset+correlation_size*2;
                curr_offset++) {
            vector<FrameInfo> frames_from_current_list(current_list.begin()+curr_offset,
                    current_list.begin()+curr_offset+correlation_size);
            vector<clock_t> timestamps_from_current_list = pull_frame_timestamps(
                        frames_from_current_list);
            double fit_val =0;
            for (size_t t=0; t!=correlation_size;
                    ++t) fit_val+= (timestamps_from_current_list[t]-offset_timestamps_from_first_list[t])*
                                       (timestamps_from_current_list[t]-offset_timestamps_from_first_list[t]);
            if (fit_val<min_fit) {
                best_matching_shift=curr_offset;
                min_fit=fit_val;
            }
        }

        int shift  = (int)first_list_corr_offset - best_matching_shift ;

        cout<<"best offset " <<best_matching_shift<<endl;
        cout <<"List "<<i <<"(CCD "<<
             current_list.front().ccd_id<<") best matches the first frame (file 0) of list 0, at file number "
             << -shift <<endl;


        string timestr = file_utils::get_file_last_access_timestring(
                             current_list.front().file_path);
        cout <<"CCD: "<<current_list.front().ccd_id<<" zeroed at header id "
             <<(current_list.front().corrected_frame_id - shift)
             <<" ( header timestamp: "<<current_list.front().header_timestamp<<")"
             <<"( filesystem timestamp: "<<timestr<<")"
             <<endl;
        cout <<"Correcting frame ids with shift of " << -1*
             (current_list.front().corrected_frame_id - shift)<<endl;
        FrameInfo::shift_frame_id_numbers(-1*(current_list.front().corrected_frame_id - shift),
                                          current_list);
    }
}

vector<MultiFrame> MultiFrame::collate_frames_present_for_all_ccds(
    const vector< vector<FrameInfo> >& lists)
{
    if (lists.size()==1) {
        //create a vector of multi frames, each of size 1 frame.
        vector<MultiFrame> multi_frm_vec;
        multi_frm_vec.reserve(lists.at(0).size());
        for (size_t i=0; i!=lists.front().size(); ++i) {
            const vector<FrameInfo>& current_list = lists.front();
            MultiFrame tmp_mlt_frm(current_list[i].corrected_frame_id);
            tmp_mlt_frm.synchronized_CCD_frames.push_back(current_list[i]);
//                        tmp_mlt_frm.multi_frame_has_been_drizzled = false; //done by constructor
            multi_frm_vec.push_back(tmp_mlt_frm);
        }
        return multi_frm_vec;
    }
    //(else...)
    vector<int> intersecting_IDs = pull_corrected_frame_ids(lists.front());
    for (size_t list_num=1; list_num!=lists.size(); list_num++) {
        intersecting_IDs = intersecting_integers(intersecting_IDs,
                           pull_corrected_frame_ids(lists[list_num]));
    }
    //Ok, so intersection now holds the corrected frame IDs of the frames present for all lists.
    //Now to collate the multi-frame structures.

    assert(!intersecting_IDs.empty());

    //Create the MultiFrame holders with the correct frame IDs
    vector<MultiFrame> multi_frm_vec;
    multi_frm_vec.reserve(intersecting_IDs.size());
    for (size_t intersection_count=0; intersection_count!=intersecting_IDs.size();
            intersection_count++) {
        multi_frm_vec.push_back(MultiFrame(intersecting_IDs[intersection_count]));
    }

    //Cycle through the frame_lists appending the relevant frames to each MultiFrame holder
    //NB!!! Assumes frame_lists sorted in ascending corrected header id!
    for (size_t list_num=0; list_num!=lists.size(); list_num++) {
        const vector<FrameInfo> & current_list = lists[list_num];
        for (size_t i=0; i!=multi_frm_vec.size(); ++i) {
            size_t list_pos=0;
            while (current_list[list_pos].corrected_frame_id!=multi_frm_vec[i].frame_id
                    &&list_pos!=current_list.size()) { list_pos++; }
            if (list_pos==current_list.size()) { throw logic_error("Suspect lists not sorted correctly for \"Corrected_frame_id\" matching"); }
            multi_frm_vec[i].synchronized_CCD_frames.push_back(current_list[list_pos]);
        }
    }
    return multi_frm_vec;
}

MultiFrame MultiFrame::pull_multi_frame_with_id(int frame_id,
        const std::vector<MultiFrame>& frms_vec)
{
    for (size_t i=0; i!=frms_vec.size(); ++i) {
        if (frms_vec[i].frame_id==frame_id) { return frms_vec[i]; }
    }
    throw runtime_error("No matching id");
}

vector<FrameInfo> MultiFrame::retrieve_frames_for_CCD_id(int CCD_id,
        const std::vector<MultiFrame>& frms_vec)
{
    vector<FrameInfo> single_ccd_frms;
    assert(!frms_vec.empty());
    int relevant_frame_index=-1;
    for (size_t i=0; i!=frms_vec.front().synchronized_CCD_frames.size(); ++i) {
        if (frms_vec.front().synchronized_CCD_frames[i].ccd_id == CCD_id) { relevant_frame_index=i; }
    }
    if (relevant_frame_index==-1) { throw runtime_error("CCD id not found in multi-frame vec: "+ itoa(CCD_id)); }
    for (size_t i=0; i!=frms_vec.size(); ++i) {
        single_ccd_frms.push_back(frms_vec[i].synchronized_CCD_frames[relevant_frame_index]);
    }

    return vector<FrameInfo>(single_ccd_frms.begin(), single_ccd_frms.end());
}


vector<MultiFrame> MultiFrame::get_multi_frame_vec(const vector<FrameInfo>& gs_frm_vec,
        const vector< vector<FrameInfo> >& additional_frame_lists,
        bool perform_timestamp_correlation_and_correction
                                                  )
{
    vector< vector <FrameInfo> > all_frm_inf_vecs;
    all_frm_inf_vecs.push_back(gs_frm_vec);

    for (size_t i=0; i!=additional_frame_lists.size(); ++i) {
        all_frm_inf_vecs.push_back(additional_frame_lists[i]);
    }

    cout<<"Frame lists loaded;"<<endl;

    //Correlate timestamps in case the frame ids are out of step.
    if (all_frm_inf_vecs.front().size()>MultiFrame::correlation_size*4
            && perform_timestamp_correlation_and_correction==true) {
        // (only do this if using real data)
        cout<<"Correlating "<<all_frm_inf_vecs.size()<<" frame lists"<<endl;
        MultiFrame::correlate_timestamps_and_correct_frame_IDs(all_frm_inf_vecs);
    }

    //Collate the frame lists into  a vector of "multi-frames", structs for holding all the frame_infos relating to one exposure time frame.
    cout<<"Collating "<<all_frm_inf_vecs.size()<<" frame lists"<<endl;
    return MultiFrame::collate_frames_present_for_all_ccds(all_frm_inf_vecs);
}

}


