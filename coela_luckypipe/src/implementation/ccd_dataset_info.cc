#include "../ccd_dataset_info.h"
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

CcdDatasetInfo::CcdDatasetInfo(const string& dataset_info_file)
{
    simple_serialization::generic_load_from_file_stream(dataset_info_file, *this);
}

void CcdDatasetInfo::write_to_file(const string& filename) const
{
    simple_serialization::generic_output_to_file_stream(filename, *this);
}


const string CcdDatasetInfo::serialization_start_key = "CcdDatasetInfo Ver 1.0";
const string CcdDatasetInfo::serialization_end_key = "End of CcdDatasetInfo";
const string CcdDatasetInfo::ccd_id_key = "CCD id number";
const string CcdDatasetInfo::CCD_inputdir_key = "Raw data dir for this CCD";
const string CcdDatasetInfo::CCD_outputdir_key = "Output dir for this CCD";
const string CcdDatasetInfo::dataset_output_base_dir_key =
    "Output dir for this observation";
const string CcdDatasetInfo::CcdDatasetInfo::in_filestem_key = "Raw data filestem";
const string CcdDatasetInfo::CcdDatasetInfo::out_filestem_key =
    "Processed data filestem";
const string CcdDatasetInfo::extension_key = "File extension";

const string CcdDatasetInfo::default_camera_config_file_key =
    "Default camconf file location";
const string CcdDatasetInfo::default_guide_star_region_file_key =
    "Default guide star region file location";
const string CcdDatasetInfo::default_faint_histogram_region_file_key =
    "Default faint histogram region file location";
const string CcdDatasetInfo::uniform_bias_pedestal_estimate_key =
    "Uniform bias pedestal estimate";


const string CcdDatasetInfo::most_recently_output_frame_list_key =
    "Most recently output frame list";
const string CcdDatasetInfo::most_recently_output_histogram_key =
    "Most recently output histogram";
const string CcdDatasetInfo::most_recently_used_column_bias_frame_key =
    "Most recently used column bias frame";

std::ostream& operator<<(std::ostream& os, const CcdDatasetInfo& cdi)
{
    using namespace simple_serialization;
    using namespace string_utils;
    os << CcdDatasetInfo::serialization_start_key << simple_serialization::new_line;

    output_key_value_pair_to_stream(os, CcdDatasetInfo::ccd_id_key, cdi.ccd_id);
    output_key_value_pair_to_stream(os, CcdDatasetInfo::CCD_inputdir_key, cdi.CCD_inputdir);
    output_key_value_pair_to_stream(os, CcdDatasetInfo::CCD_outputdir_key,
                                    cdi.CCD_outputdir);
    output_key_value_pair_to_stream(os, CcdDatasetInfo::dataset_output_base_dir_key,
                                    cdi.dataset_output_base_dir);
    output_key_value_pair_to_stream(os, CcdDatasetInfo::in_filestem_key, cdi.in_filestem);
    output_key_value_pair_to_stream(os, CcdDatasetInfo::out_filestem_key, cdi.out_filestem);
    output_key_value_pair_to_stream(os, CcdDatasetInfo::extension_key, cdi.extension);

    output_key_value_pair_to_stream(os, CcdDatasetInfo::default_camera_config_file_key,
                                    cdi.default_camera_config_file);
    output_key_value_pair_to_stream(os, CcdDatasetInfo::uniform_bias_pedestal_estimate_key,
                                    cdi.uniform_bias_pedestal_estimate);

    output_key_value_pair_to_stream(os, CcdDatasetInfo::default_guide_star_region_file_key,
                                    cdi.default_guide_star_region_file);
    output_key_value_pair_to_stream(os,
                                    CcdDatasetInfo::default_faint_histogram_region_file_key,
                                    cdi.default_faint_histogram_region_file);
    os<<"\n";
    output_key_value_pair_to_stream(os, CcdDatasetInfo::most_recently_output_frame_list_key,
                                    cdi.most_recently_output_frame_list);
    output_key_value_pair_to_stream(os, CcdDatasetInfo::most_recently_output_histogram_key,
                                    cdi.most_recently_output_histogram);
    output_key_value_pair_to_stream(os,
                                    CcdDatasetInfo::most_recently_used_column_bias_frame_key,
                                    cdi.most_recently_used_column_bias_frame);

    os << CcdDatasetInfo::serialization_end_key << simple_serialization::new_line;

    return os;
}

std::istream& operator>>(std::istream& is, CcdDatasetInfo& cdi)
{
    using namespace simple_serialization;
    using namespace string_utils;

    vector<string> svec = convert_relevant_istream_section(is,
                          CcdDatasetInfo::serialization_start_key,
                          CcdDatasetInfo::serialization_end_key);

    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::ccd_id_key, cdi.ccd_id);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::CCD_inputdir_key,
                                    cdi.CCD_inputdir);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::CCD_outputdir_key,
                                    cdi.CCD_outputdir);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::dataset_output_base_dir_key,
                                    cdi.dataset_output_base_dir);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::in_filestem_key, cdi.in_filestem);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::out_filestem_key,
                                    cdi.out_filestem);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::extension_key, cdi.extension);

    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::default_camera_config_file_key,
                                    cdi.default_camera_config_file);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::uniform_bias_pedestal_estimate_key,
                                    cdi.uniform_bias_pedestal_estimate);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::default_guide_star_region_file_key,
                                    cdi.default_guide_star_region_file);
    get_keyed_value_from_string_vec(svec,
                                    CcdDatasetInfo::default_faint_histogram_region_file_key,
                                    cdi.default_faint_histogram_region_file);

    get_keyed_value_from_string_vec(svec,
                                    CcdDatasetInfo::most_recently_output_frame_list_key,
                                    cdi.most_recently_output_frame_list);
    get_keyed_value_from_string_vec(svec, CcdDatasetInfo::most_recently_output_histogram_key,
                                    cdi.most_recently_output_histogram);
    get_keyed_value_from_string_vec(svec,
                                    CcdDatasetInfo::most_recently_used_column_bias_frame_key,
                                    cdi.most_recently_used_column_bias_frame);
    return is;

}


}//end namespace coela
