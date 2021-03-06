/*
 * File:   FileInfo.h
 * Author: ts337
 *
 * Created on 03 April 2009, 16:39
 */

#ifndef COELA_FileInfo_H
#define COELA_FileInfo_H

#include <string>
#include <list>
#include <vector>
#include <list>
#include <istream>
#include "coela_utility/src/file_utils.h"
#include "coela_core/src/fits_header.h"

//#include "../level1/lucky_file_buffer.h"
using std::string;
///This stores information about where a frame is stored;
///whether the frame is a single file or part of an LCM.
namespace coela {

//============================================================================
///struct FileInfo:
///Basic collection of information about an image file

///Can represent:
/// -1 image per file
/// -image files as subfiles of a tarball
/// -3d FITS cube.

/// NB this only employs the "FitsHeader" (very stable / unlikely to change
///interface), agnostic to details of header classes.

struct FileInfo { //I wonder if I should privatise any of this...
    FileInfo();

    ///constructor for a file NOT packed into an lcm/lcz custom format:
    static FileInfo basic_fits_file(const std::string& fname);

    //--------------------------------------------------------------------------------------------
    //data members:
    std::string file_path, derived_lcc_filename;
    size_t header_byte_offset, byte_size;
    bool file_is_cube_FITS;
    int ccd_id, header_frame_id, corrected_frame_id;
    int folder_number, subfile_number;
    unsigned short int total_files_in_lcm;
    clock_t header_timestamp;
    //--------------------------------------------------------------------------------------------

    void load_key_vals_from_header_table(const FitsHeader &);

    ///Appends all the files in the lcm_list to the subfiles list
    static void get_subfiles_info(std::list<string>& lcm_list,
            std::list<FileInfo>& subfiles, const std::string& extension);
    static void get_lcm_subfiles_info(const std::list<string>& lcm_list,
            std::list<FileInfo>& subfiles);

    static std::vector<size_t> parse_lcm_header(const std::string& filename);

    static bool filename_number_predicate(const string& first_filename,
            const string& second_filename);
    static bool subfile_number_predicate(const FileInfo& first,
            const FileInfo& second);

    /// Is first.frm_id < second.frm_id ?
    static bool frame_id_predicate(const FileInfo& first,
            const FileInfo& second);
    static bool subfile_numbers_equal(const FileInfo& first,
            const FileInfo& second);

    static std::list<FileInfo> get_image_file_list(const std::string& dir,
            const std::string& filestem, int CCD_id,
            const std::string& extension = "");
    static int read_num_files_in_lcm(const string& lcm_filename);

    static std::string generate_unpacked_lcc_filename(
            const std::string& filename, int file_num_in_lcm, int files_in_lcm);

    static FileInfo load_from_lcz_index_line(const std::string& index_line,
            const std::string& base_folder, const std::string& ccd_stem,
            const std::string& subfolder_stem, const std::string& filestem,
            const std::string& file_extension);

    static FileInfo load_from_lcz_index_line_segments(
            const vector<string>& line_segments, const std::string& base_folder,
            const std::string& ccd_stem, const std::string& subfolder_stem,
            const std::string& filestem, const std::string& file_extension);

    static void convert_paths_to_system_complete(std::list<FileInfo>&);

};
std::ostream& operator<<(std::ostream& os, const FileInfo&); ///<Output in index form.

//==================================================================================================================================

namespace lcz_utils {
string mangle_lcz_filename(const string& filename);

bool check_headers_match_index(const std::list<FileInfo>&);
int check_frames_stored_sequentially(const std::list<FileInfo>&);
int check_for_dropped_frames(std::list<FileInfo>&);
int check_for_duplicate_frame_ids(std::list<FileInfo>&);
void remove_duplicate_frame_ids(std::list<FileInfo>&,
        const bool rigorous_check = false);
void set_correct_frame_ids(std::list<FileInfo>&);
void sanitise_list(std::list<FileInfo>&);

std::list<FileInfo> get_lcz_subfile_list(const std::string& ccd_dir,
        const std::string& index_output_dir = "");

std::list<FileInfo> build_lcz_subfile_index(const std::list<string>& lcz_list);
bool empty_files_present(const std::list<string>& lcz_list);
std::list<FileInfo> parse_lcz_file(const std::string& lcz_filename);

std::list<FileInfo> parse_lcz_index(const string& index_filename,
        const string& data_root_dir = "");
bool check_lcz_index_is_standard_format(const string& index_filename);
std::list<int> get_lcz_nums(const std::string& index_filename);
std::list<int> get_lcz_nums(const std::vector<std::string>& index_file_text);
void write_lcz_index(const std::list<FileInfo>& subfiles,
        const string& index_filename);
void write_lcz_index(const std::list<FileInfo>& subfiles,
        const string& index_filename, const string& root_folder,
        const string& cam_folder_stem, const string& sub_folder_stem,
        const string& filestem);

//void get_lcz_subfiles_info(std::list<string>& lcm_list, std::list<FileInfo>& subfiles); //deprecated

}//end namespace lcz_utils

//==================================================================================================================================
}//end namespace coela

#endif  /* _FileInfo_H */

