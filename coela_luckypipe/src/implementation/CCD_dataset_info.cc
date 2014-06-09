#include "../CCD_dataset_info.h"
#include "coela_utility/src/string_utils.h"
#include "coela_utility/src/file_utils.h"
#include "coela_utility/src/simple_serialization.h"

#include <boost/filesystem.hpp>
#include <stdexcept>
#include <iostream>
using std::string;
using std::runtime_error;
using namespace std;

namespace coela {

CCD_DatasetInfo::CCD_DatasetInfo(const string& dataset_info_file)
{
    simple_serialization::generic_load_from_file_stream(dataset_info_file, *this);
}

void CCD_DatasetInfo::write_to_file(const string& filename) const
{
    simple_serialization::generic_output_to_file_stream(filename, *this);
}


const string CCD_DatasetInfo::serialization_start_key = "CCD_DatasetInfo Ver 1.0";
const string CCD_DatasetInfo::serialization_end_key = "End of CCD_DatasetInfo";
const string CCD_DatasetInfo::ccd_id_key = "CCD id number";
const string CCD_DatasetInfo::CCD_inputdir_key = "Raw data dir for this CCD";
const string CCD_DatasetInfo::CCD_outputdir_key = "Output dir for this CCD";
const string CCD_DatasetInfo::dataset_output_base_dir_key =
    "Output dir for this observation";
const string CCD_DatasetInfo::CCD_DatasetInfo::in_filestem_key = "Raw data filestem";
const string CCD_DatasetInfo::CCD_DatasetInfo::out_filestem_key =
    "Processed data filestem";
const string CCD_DatasetInfo::extension_key = "File extension";

const string CCD_DatasetInfo::default_camera_config_file_key =
    "Default camconf file location";
const string CCD_DatasetInfo::default_guide_star_region_file_key =
    "Default guide star region file location";
const string CCD_DatasetInfo::default_faint_histogram_region_file_key =
    "Default faint histogram region file location";
const string CCD_DatasetInfo::uniform_bias_pedestal_estimate_key =
    "Uniform bias pedestal estimate";


const string CCD_DatasetInfo::most_recently_output_frame_list_key =
    "Most recently output frame list";
const string CCD_DatasetInfo::most_recently_output_histogram_key =
    "Most recently output histogram";
const string CCD_DatasetInfo::most_recently_used_column_bias_frame_key =
    "Most recently used column bias frame";

std::ostream& operator<<(std::ostream& os, const CCD_DatasetInfo& cdi)
{
    using namespace simple_serialization;
    using namespace string_utils;
    os << CCD_DatasetInfo::serialization_start_key << simple_serialization::new_line;

    output_key_value_pair_to_stream(os, CCD_DatasetInfo::ccd_id_key, cdi.ccd_id);
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::CCD_inputdir_key, cdi.CCD_inputdir);
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::CCD_outputdir_key,
                                    cdi.CCD_outputdir);
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::dataset_output_base_dir_key,
                                    cdi.dataset_output_base_dir);
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::in_filestem_key, cdi.in_filestem);
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::out_filestem_key, cdi.out_filestem);
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::extension_key, cdi.extension);

    output_key_value_pair_to_stream(os, CCD_DatasetInfo::default_camera_config_file_key,
                                    cdi.default_camera_config_file);
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::uniform_bias_pedestal_estimate_key,
                                    cdi.uniform_bias_pedestal_estimate);

    output_key_value_pair_to_stream(os, CCD_DatasetInfo::default_guide_star_region_file_key,
                                    cdi.default_guide_star_region_file);
    output_key_value_pair_to_stream(os,
                                    CCD_DatasetInfo::default_faint_histogram_region_file_key,
                                    cdi.default_faint_histogram_region_file);
    os<<"\n";
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::most_recently_output_frame_list_key,
                                    cdi.most_recently_output_frame_list);
    output_key_value_pair_to_stream(os, CCD_DatasetInfo::most_recently_output_histogram_key,
                                    cdi.most_recently_output_histogram);
    output_key_value_pair_to_stream(os,
                                    CCD_DatasetInfo::most_recently_used_column_bias_frame_key,
                                    cdi.most_recently_used_column_bias_frame);

    os << CCD_DatasetInfo::serialization_end_key << simple_serialization::new_line;

    return os;
}

std::istream& operator>>(std::istream& is, CCD_DatasetInfo& cdi)
{
    using namespace simple_serialization;
    using namespace string_utils;

    vector<string> svec = convert_relevant_istream_section(is,
                          CCD_DatasetInfo::serialization_start_key,
                          CCD_DatasetInfo::serialization_end_key);

    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::ccd_id_key, cdi.ccd_id);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::CCD_inputdir_key,
                                    cdi.CCD_inputdir);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::CCD_outputdir_key,
                                    cdi.CCD_outputdir);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::dataset_output_base_dir_key,
                                    cdi.dataset_output_base_dir);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::in_filestem_key, cdi.in_filestem);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::out_filestem_key,
                                    cdi.out_filestem);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::extension_key, cdi.extension);

    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::default_camera_config_file_key,
                                    cdi.default_camera_config_file);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::uniform_bias_pedestal_estimate_key,
                                    cdi.uniform_bias_pedestal_estimate);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::default_guide_star_region_file_key,
                                    cdi.default_guide_star_region_file);
    get_keyed_value_from_string_vec(svec,
                                    CCD_DatasetInfo::default_faint_histogram_region_file_key,
                                    cdi.default_faint_histogram_region_file);

    get_keyed_value_from_string_vec(svec,
                                    CCD_DatasetInfo::most_recently_output_frame_list_key,
                                    cdi.most_recently_output_frame_list);
    get_keyed_value_from_string_vec(svec, CCD_DatasetInfo::most_recently_output_histogram_key,
                                    cdi.most_recently_output_histogram);
    get_keyed_value_from_string_vec(svec,
                                    CCD_DatasetInfo::most_recently_used_column_bias_frame_key,
                                    cdi.most_recently_used_column_bias_frame);
    return is;

}


}//end namespace coela
