/*
 * File:   CCD_dataset_info.h
 * Author: ts337
 *
 * Created on 15 March 2011, 12:01
 */

#ifndef COELA_CCD_DATASET_INFO_H
#define COELA_CCD_DATASET_INFO_H
#include <string>
#include <fstream>
#include <vector>

namespace coela {

struct CCD_DatasetInfo {
    CCD_DatasetInfo():ccd_id(-1), uniform_bias_pedestal_estimate(0.0) {}

    CCD_DatasetInfo(const std::string& dataset_info_filename);
    void write_to_file(const std::string& filename) const;


    int ccd_id;

    std::string CCD_inputdir, CCD_outputdir,
        in_filestem, out_filestem,
        extension,
        dataset_output_base_dir;

    std::string default_camera_config_file;
    std::string default_guide_star_region_file;
    std::string default_faint_histogram_region_file;
    double uniform_bias_pedestal_estimate;


    std::string most_recently_output_frame_list;
    std::string most_recently_output_histogram;
    std::string most_recently_used_column_bias_frame;
    std::string most_recently_used_row_bias_frame;


private:
    friend std::ostream& operator<<(std::ostream& os, const CCD_DatasetInfo& cdi);
    friend std::istream& operator>>(std::istream& is, CCD_DatasetInfo& cdi);

    void load_from_text_vec(const std::vector<std::string>&);

    //---------------------------------------------------------------
    //keys for serialization:
    static const std::string serialization_start_key, serialization_end_key;
    static const std::string ccd_id_key;
    static const std::string CCD_inputdir_key, CCD_outputdir_key,
           in_filestem_key, out_filestem_key,
           extension_key,
           dataset_output_base_dir_key;


    static const std::string default_camera_config_file_key;
    static const std::string default_guide_star_region_file_key;
    static const std::string default_faint_histogram_region_file_key;
    static const std::string uniform_bias_pedestal_estimate_key;


    static const std::string most_recently_output_frame_list_key;
    static const std::string most_recently_output_histogram_key;
    static const std::string most_recently_used_column_bias_frame_key;
    //---------------------------------------------------------------
};



}//end namespace coela
#endif  /* CCD_DATASET_INFO_H */

