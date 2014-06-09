#include "../camera_config_info.h"

#include "coela_utility/src/string_utils.h"
#include "coela_utility/src/file_utils.h"
#include "coela_core/src/DS9Region.h"
#include "coela_core/src/cartesian_coords.h"
#include "coela_utility/src/simple_serialization.h"
//#include "../level2/file_info.h"
//#include "../level1/ds9_interface.h"
#include <boost/filesystem.hpp>
#include <ctime>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>



//using string_utils::atoi;
using namespace std;
namespace bfs=boost::filesystem;
//using namespace string_utils;
namespace coela {
using string_utils::atof;
using namespace string_utils;
using std::endl; using std::cout; using std::cerr;
using std::runtime_error; using std::logic_error;
using std::ifstream; using std::ostream;



//============================================================================================================
std::ostream& operator<<(std::ostream& os,
                         const vector<CCD_calibration_info>& ccd_chars_vec)
{
    os<<"#CCD characteristics information file:\n";
    os<<"#Version 1.1\n\n";
//        os<<"Crop raw CCD data?(true/false): " <<bool_to_string(ccd_inf.crop_CCD_data)<<"\n";
    os<<"Info present for CCD IDs : ";
    for (size_t i=0; i!=ccd_chars_vec.size(); ++i) {
        os<<ccd_chars_vec[i].ccd_id<<",";
    }
    os<<"\n\n";

    for (size_t i=0; i!=ccd_chars_vec.size(); ++i) {
        const CCD_calibration_info& ccd_chars=ccd_chars_vec[i];

        os << "Bad cols for CCD_ID " <<ccd_chars.ccd_id<<": \t";
        for (size_t i=0; i!=ccd_chars.bad_columns.size(); ++i) { os << ccd_chars.bad_columns[i]<<","; }
        os<<"\n";
        os<<"Crop region for CCD_ID "<<
          ccd_chars.ccd_id<<": \t"<<ccd_chars.cropped_PixelRange<<"\n";
        os<<"Default temporal debiasing histogram region for CCD_ID "
          << ccd_chars.ccd_id<<": \t"
          <<ccd_chars.default_temporal_debiasing_histogram_region<<"\n";

        os<<"Row Bias frame available for CCD_ID "<< ccd_chars.ccd_id<<": \t"
          << bool_to_string(ccd_chars.precal_row_bias_frame_available)<<"\n";
        os<<"Row Bias frame locale for CCD_ID "<< ccd_chars.ccd_id<<": \t"
          << ccd_chars.precal_row_bias_frame_path<<"\n";
        os<<"Dark current frame available for CCD_ID "<< ccd_chars.ccd_id<<": \t"
          << bool_to_string(ccd_chars.dark_current_frames_available)<<"\n";
        os<<"Normalised DC frame locale for CCD_ID "<< ccd_chars.ccd_id<<": \t"
          << ccd_chars.normalised_DC_frame_path<<"\n";
        os<<"Thresholded DC frame locale for CCD_ID "<< ccd_chars.ccd_id<<": \t"
          << ccd_chars.normalised_DC_frame_path<<"\n";
        os<<"Frame Exposure Time (milliseconds) for CCD_ID "<< ccd_chars.ccd_id<<": \t"
          << ccd_chars.milliseconds_per_frame_exposure<<"\n";
        os<<endl;

    }
    return os;
}

vector<CCD_calibration_info> CCD_calibration_info::load_from_text_vec(
    const vector<string>& file_text)
{
    using namespace string_utils;
    vector<CCD_calibration_info> ccd_vec;

    vector<string> ccd_id_strs=tokenize(parse_for_key_value(file_text,
                                        "Info present for CCD IDs"), ",");
    for (size_t i=0; i!=ccd_id_strs.size(); ++i) {

        CCD_calibration_info ccd_props;
        ccd_props.ccd_id=atoi(ccd_id_strs[i]);

        ccd_props.precal_row_bias_frame_available = string_to_bool(parse_for_key_value(file_text,
                "Row Bias frame available for CCD_ID "+ccd_id_strs[i]));
        ccd_props.precal_row_bias_frame_path=parse_for_key_value(file_text,
                                             "Row Bias frame locale for CCD_ID "+ccd_id_strs[i]);
        ccd_props.dark_current_frames_available = string_to_bool(parse_for_key_value(file_text,
                "Dark current frame available for CCD_ID "+ccd_id_strs[i]));
        ccd_props.normalised_DC_frame_path = parse_for_key_value(file_text,
                                             "Normalised DC frame locale for CCD_ID "+ccd_id_strs[i]);
        ccd_props.thresholded_DC_frame_path = parse_for_key_value(file_text,
                                              "Thresholded DC frame locale for CCD_ID "+ccd_id_strs[i]);
        ccd_props.milliseconds_per_frame_exposure = atof(parse_for_key_value(file_text,
                "Frame Exposure Time (milliseconds) for CCD_ID " +ccd_id_strs[i]));

        string bad_cols_str=parse_for_key_value(file_text, "Bad cols for CCD_ID "+ccd_id_strs[i]);
        vector<string> bad_col_strs = tokenize(bad_cols_str, ",");
        for (size_t k=0; k!=bad_col_strs.size(); ++k) { ccd_props.bad_columns.push_back(atoi(bad_col_strs[k])); }

        string crop_rgn_str=parse_for_key_value(file_text,
                                                "Crop region for CCD_ID "+ccd_id_strs[i]);
        stringstream ss(crop_rgn_str);
        ss>>ccd_props.cropped_PixelRange;
        ss.clear();



        string bg_rgn_str=parse_for_key_value(file_text,
                                              "Default temporal debiasing histogram region for CCD_ID "+ccd_id_strs[i]);
        ss.str(bg_rgn_str);
        ss >> ccd_props.default_temporal_debiasing_histogram_region;


        ccd_props.crop_region =
            CCD_BoxRegion(ccd_props.cropped_PixelRange.low.x-1, ccd_props.cropped_PixelRange.low.y-1,
                          ccd_props.cropped_PixelRange.high.x,ccd_props.cropped_PixelRange.high.y);

        ccd_vec.push_back(ccd_props);
    }
    return ccd_vec;
}

vector<CCD_BoxRegion> CCD_calibration_info::get_bad_detector_regions() const
{
    vector<CCD_BoxRegion> bad_regions;
    for (size_t i=0; i!=bad_columns.size(); ++i) {
        bad_regions.push_back(CCD_BoxRegion(bad_columns[i]-1.0,cropped_PixelRange.low.y-1.0,
                                            bad_columns[i], cropped_PixelRange.high.y));
    }
    return bad_regions;
}

//============================================================================================================


std::ostream& operator<<(std::ostream& os, const lens_and_aperture_info& lens_inf)
{
    os<<"#Lens and aperture information file \n";
    os<<"#Version 1.0\n\n";
    os<<"Telescope outer diameter: " <<ftoa(lens_inf.telescope_outer_diameter)<<"\n";
    os<<"Telescope inner diameter: " <<ftoa(lens_inf.telescope_inner_diameter)<<"\n";
    os<<"Pixel scale in milliarcseconds: " <<ftoa(lens_inf.nominal_pixel_scale_in_mas)<<"\n";
    os<<"\n#CCD sky configuration (in CCD pixel units)\n";
    for (size_t i=0; i!=lens_inf.ccd_mosaic_offsets.size(); ++i) {
        os<<"Offset for CCD_ID: " <<lens_inf.ccd_mosaic_offsets[i].first<<";\t"
          <<lens_inf.ccd_mosaic_offsets[i].second<<"\n";
    }
    return os;
}

void lens_and_aperture_info::load_from_text_vec(const vector<string>& file_text)
{
    telescope_outer_diameter = atof(parse_for_key_value(file_text,
                                    "Telescope outer diameter"));
    telescope_inner_diameter = atof(parse_for_key_value(file_text,
                                    "Telescope inner diameter"));
    nominal_pixel_scale_in_mas = atof(parse_for_key_value(file_text,
                                      "Pixel scale in milliarcseconds"));

    for (size_t i=0; i!=file_text.size(); ++i) {
        const string& line = file_text[i];
        if (line.find("Offset for CCD")!=string::npos) {
            vector<string> sections = tokenize(line,";");
            assert(sections.size()==2);
            int line_ccd_id = atoi(tokenize(sections.front(),":").back());
            MosaicPixelShift line_offset;
            stringstream ss;
            ss.str(sections.back());
            ss>>line_offset;
            ccd_mosaic_offsets.push_back(pair<int, MosaicPixelShift>(line_ccd_id,line_offset));
        }
    }
}

//============================================================================================================

vector<OpticalFilterInfo> OpticalFilterInfo::load_from_text_vec(
    const vector<string>& file_text)
{
    using namespace string_utils;
    vector<string> ccd_id_strs=tokenize(parse_for_key_value(file_text,
                                        "Info present for CCD IDs"), ",");

    vector<OpticalFilterInfo> filters;

    for (size_t i=0; i!=ccd_id_strs.size(); ++i) {
        OpticalFilterInfo fb;
        fb.ccd_id=atoi(ccd_id_strs[i]);
        fb.central_wavelength_in_metres = atof(parse_for_key_value(file_text,
                                               "Central filter wavelength for CCD ID "+ccd_id_strs[i]));
        fb.name = parse_for_key_value(file_text, "Filter name for CCD ID "+ccd_id_strs[i]);
        filters.push_back(fb);
    }
    return filters;
}

std::ostream& operator<<(std::ostream& os, const vector<OpticalFilterInfo>& filters)
{
    using namespace string_utils;

    os<<"#Filter set information file:\n";
    os<<"#Version 1.0\n\n";

    os<<"Info present for CCD IDs: ";
    for (size_t i=0; i!=filters.size(); ++i) {
        os<<filters[i].ccd_id<<",";
    }
    os<<"\n\n";
    for (size_t i=0; i!=filters.size(); ++i) {
        os<<"Central filter wavelength for CCD ID "<< filters[i].ccd_id << " : "
          << filters[i].central_wavelength_in_metres<<"\n";
        os<<"Filter name for CCD ID "<< filters[i].ccd_id << " : "
          << filters[i].name<<"\n";
    }

    return os;
}

//============================================================================================================
void CameraConfigInfo::write_to_files(const string& camconf_filename,
                                      string CCD_inf_rel_path,
                                      string lens_inf_rel_path,
                                      string filter_inf_rel_path
                                     )
{
    string camconf_namestem = bfs::path(
                                  string_utils::strip_file_extension(camconf_filename)).filename().string();
    bfs::path camconf_root = bfs::path(camconf_filename).parent_path();
    if (CCD_inf_rel_path.empty()) {
        CCD_inf_rel_path = camconf_namestem+ "_CCD_inf.txt";
    }
    if (lens_inf_rel_path.empty()) {
        lens_inf_rel_path =camconf_namestem+ "_lens_inf.txt";
    }
    if (filter_inf_rel_path.empty()) {
        filter_inf_rel_path =camconf_namestem+ "_filter_inf.txt";
    }

//        CCD_filename=boost::filesystem::system_complete(CCD_filename).string();
//        lens_filename=boost::filesystem::system_complete(lens_filename).string();
//        filter_filename=boost::filesystem::system_complete(filter_filename).string();

    ofstream camconf_file(camconf_filename.c_str(), ios::binary);
    if (!camconf_file.is_open()) {
        throw runtime_error(
            "camera_config_info - cannot write to file "+camconf_filename);
    }

    camconf_file<<"#Camera configuration redirection file \n";
    camconf_file<<"#Version 1.4\n\n";
    camconf_file<<"CCD info file : "<<CCD_inf_rel_path<<"\n";
    camconf_file<<"Lens info file : "<<lens_inf_rel_path<<"\n";
    camconf_file<<"Filter info file : "<<filter_inf_rel_path<<"\n";
    camconf_file<<"Simulated data : "<<
                string_utils::bool_to_string(simulated_data)<<"\n";

    camconf_file.close();

    bfs::path CCD_inf_abs_path = bfs::absolute(CCD_inf_rel_path, camconf_root);
    bfs::path lens_inf_abs_path = bfs::absolute(lens_inf_rel_path, camconf_root);
    bfs::path filter_inf_abs_path = bfs::absolute(filter_inf_rel_path, camconf_root);

    simple_serialization::generic_output_to_file_stream(CCD_inf_abs_path.string(), CCD_vec);
    simple_serialization::generic_output_to_file_stream(lens_inf_abs_path.string(), lens_inf);
    simple_serialization::generic_output_to_file_stream(filter_inf_abs_path.string(),
            filters);
}


void CameraConfigInfo::load_from_text_vec(const vector<string>& file_text)
{
//        camera_config_info cam_conf;
    assert(file_text.size()>1);
    if (file_text[1]!="#Version 1.4") {
        throw runtime_error(
            "Incorrect Cam conf file version: read\""+file_text[1]+"\"");
    }
    string CCD_inf_relpath = parse_for_key_value(file_text, "CCD info file");
    string lens_inf_relpath= parse_for_key_value(file_text, "Lens info file");
    string  filter_inf_relpath = parse_for_key_value(file_text, "Filter info file");
    simulated_data= string_utils::string_to_bool(parse_for_key_value(file_text,
                    "Simulated data"));

    bfs::path camconf_root = bfs::path(self_path).parent_path();
    bfs::path CCD_inf_abs_path = bfs::absolute(CCD_inf_relpath, camconf_root);
    bfs::path lens_inf_abs_path = bfs::absolute(lens_inf_relpath, camconf_root);
    bfs::path filter_inf_abs_path = bfs::absolute(filter_inf_relpath, camconf_root);

    CCD_vec = CCD_calibration_info::load_from_text_vec(
                  file_utils::convert_text_file_to_string_vec(CCD_inf_abs_path.string()));
    lens_inf.load_from_text_vec(
        file_utils::convert_text_file_to_string_vec(lens_inf_abs_path.string()));
    filters = OpticalFilterInfo::load_from_text_vec(
                  file_utils::convert_text_file_to_string_vec(filter_inf_abs_path.string()));

    return;
}

CameraConfigInfo::CameraConfigInfo(const string& filename)
{
    self_path = filename;
    load_from_text_vec(file_utils::convert_text_file_to_string_vec(filename));
}


CCD_calibration_info& CameraConfigInfo::get_calibration_info_for_CCD_id(int ccd_id)
{
    for (size_t i=0; i!=CCD_vec.size(); ++i) {
        if (CCD_vec[i].ccd_id==ccd_id) { return CCD_vec[i]; }
    }
    throw runtime_error("CameraConfigInfo::get_calibration_info_for_CCD_id ---\n"
                        "Cannot find CCD info for this id: "+itoa(ccd_id));
}

const CCD_calibration_info& CameraConfigInfo::get_calibration_info_for_CCD_id(
    int ccd_id) const
{
    for (size_t i=0; i!=CCD_vec.size(); ++i) {
        if (CCD_vec[i].ccd_id==ccd_id) { return CCD_vec[i]; }
    }
    throw runtime_error("CameraConfigInfo::get_calibration_info_for_CCD_id ---\n"
                        "Cannot find CCD info for this id: "+itoa(ccd_id));
}

OpticalFilterInfo& CameraConfigInfo::get_filter_info_for_ccd_id(int ccd_id)
{
    for (size_t i=0; i!=filters.size(); ++i) {
        if (filters[i].ccd_id==ccd_id) { return filters[i]; }
    }
    throw runtime_error("CameraConfigInfo::get_filter_info_for_ccd_id---\n "
                        "Cannot find filter info for this id: "+itoa(ccd_id));
}

const OpticalFilterInfo& CameraConfigInfo::get_filter_info_for_ccd_id(int ccd_id) const
{
    for (size_t i=0; i!=filters.size(); ++i) {
        if (filters[i].ccd_id==ccd_id) { return filters[i]; }
    }
    throw runtime_error("CameraConfigInfo::get_filter_info_for_ccd_id---\n "
                        "Cannot find filter info for this id: "+itoa(ccd_id));
}

const MosaicPixelShift& CameraConfigInfo::get_mosaic_offset_for_ccd_id(int ccd_id) const
{
    for (size_t i=0; i!=lens_inf.ccd_mosaic_offsets.size(); ++i) {
        if (ccd_id == lens_inf.ccd_mosaic_offsets[i].first) { return lens_inf.ccd_mosaic_offsets[i].second; }
    }
    throw runtime_error("camera_config_info::get_mosaic_offset_for_ccd_id - ccd_id not found");
}
//============================================================================================================
}//end namespace coela
