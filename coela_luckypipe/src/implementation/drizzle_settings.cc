
#include "../drizzle_settings.h"
#include "coela_utility/src/string_utils.h"
#include "coela_utility/src/file_utils.h"
#include "coela_utility/src/simple_serialization.h"
#include "coela_core/src/DS9Region.h"
#include "coela_core/src/cartesian_coords.h"

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
//using namespace string_utils;

using std::endl; using std::cout; using std::cerr;
using std::runtime_error; using std::logic_error;
using std::ifstream; using std::ostream;

namespace coela {
using namespace string_utils;


std::ostream& operator<<(std::ostream& os, const multi_ccd_region& c_rgn)
{
    os<<c_rgn.CCD_id <<"   " <<(CCD_BoxRegion)c_rgn;
    return os;
}

//============================================================================================================

DrizzleSettings::DrizzleSettings():
    normalisation_on(false),
    thresholding_on(true),
    post_drizzle_col_debias_on(true),
    create_drizzle_mosaic(false),
    drizzle_scale_factor(1.0), drizzle_pixel_fraction(0.45)
{}

DrizzleSettings::DrizzleSettings(const string& filename)
{
    parse_DrizzleSettings_text_vec(file_utils::convert_text_file_to_string_vec(filename));
}

vector<CCD_BoxRegion> multi_ccd_region::pull_regions_for_cam_id(int cam_id,
        vector<multi_ccd_region>& rgns)
{
    vector<CCD_BoxRegion> relevant_rgns;
    for (size_t i=0; i!=rgns.size(); ++i) {
        if (rgns[i].CCD_id==cam_id) {
            relevant_rgns.push_back(rgns[i]);
            rgns.erase(rgns.begin()+i);
            i--;
        }
    }
    return relevant_rgns;
}

void multi_ccd_region::push_regions_with_cam_id(vector<multi_ccd_region>&
        camera_region_list, const vector<CCD_BoxRegion>& ccd_rgns ,int cam_id)
{
    for (size_t i=0; i!=ccd_rgns.size(); ++i) {
        camera_region_list.push_back(multi_ccd_region(cam_id, ccd_rgns[i]));
    }
}


//    void lucky_settings::update_gs_regions(const string& long_exposure_filename, const int camera_id){
//           vector<CCD_BoxRegion> cam_regions  = multi_ccd_region::pull_regions_for_cam_id(camera_id, guide_regions);
//           cam_regions = ds9::get_CCD_regions(long_exposure_filename,"Update GS regions" ,cam_regions);
//            //FIXME Clear the guide region here to avoid duplicated - but must only clear those relevant to cam id
//           //i.e. guide_regions.clear_cameraXregions
//           for (size_t j=0; j< cam_regions.size(); j++)
//               guide_regions.push_back(multi_ccd_region(camera_id, cam_regions[j]));
//    }
void DrizzleSettings::write_to_file(const string& filename)
{
    simple_serialization::generic_output_to_file_stream<DrizzleSettings>(filename, *this);
}

string DrizzleSettings::create_output_dirs(
    const string& gs_frames_list_filename,
    const vector<FrameInfo>& gs_frame_list,
    const vector< MultiFrame >& multi_frame_list
) const
{
    assert(!gs_frame_list.empty());
    assert(multi_frame_list.size()==gs_frame_list.size());

    //------------------------------------------------------------------------------------
    //Copy in the name of the frames list, add some suffices....
    string drizzle_subdir=string_utils::pull_filename(
                              gs_frames_list_filename); //get rid of any directory bumpf
    drizzle_subdir=string_utils::strip_file_extension(drizzle_subdir);
    drizzle_subdir+="_drizzle";

    //Denote guiding CCD:
    drizzle_subdir+="_g";
    drizzle_subdir+=string_utils::itoa(gs_frame_list.front().ccd_id);

    //Denote all drizzled CCDs:
    const MultiFrame& mf = multi_frame_list.front();
    drizzle_subdir+="_d";
    {
        vector<int> drizzled_ccd_ids;
        for (size_t i=0; i!=mf.synchronized_CCD_frames.size(); ++i) {
            drizzled_ccd_ids.push_back(mf.synchronized_CCD_frames[i].ccd_id);
        }

        sort(drizzled_ccd_ids.begin(), drizzled_ccd_ids.end());
        for (size_t i=0; i!=drizzled_ccd_ids.size(); ++i) {
            drizzle_subdir+= string_utils::itoa(drizzled_ccd_ids[i]);
        }
    }
    //denote regular, or gain normalised
    if (!normalisation_on) { drizzle_subdir+= "_raw"; }

    drizzle_subdir+= "_" + string_utils::ftoa(1.0/drizzle_scale_factor) +"x/";

    //------------------------------------------------------------------------------------
    string full_path_drizzle_dir = output_base_folder+drizzle_subdir;
    cout <<"Creating dir "<<full_path_drizzle_dir<<endl;
    boost::filesystem::create_directory(full_path_drizzle_dir);

    for (size_t i=0; i!=mf.synchronized_CCD_frames.size(); ++i) {
        boost::filesystem::create_directory(full_path_drizzle_dir+"/CCD"+string_utils::itoa(
                                                mf.synchronized_CCD_frames[i].ccd_id));
        FrameInfo::write_list_to_file(
            MultiFrame::retrieve_frames_for_CCD_id(mf.synchronized_CCD_frames[i].ccd_id,
                    multi_frame_list),
            full_path_drizzle_dir+"/CCD"+string_utils::itoa(mf.synchronized_CCD_frames[i].ccd_id)
            +"/frames_list.txt");
    }

    return full_path_drizzle_dir;
}


string DrizzleSettings::get_CCD_output_folder_for_CCD_id(int ccd_id) const
{
    for (size_t i=0; i!=CCD_ids_present.size(); ++i) {
        if (CCD_ids_present[i]==ccd_id) { return output_base_folder+"/"+CCD_output_sub_folders[i]; }
    }
    throw logic_error("DrizzleSettings::get_CCD_output_folder_for_CCD_id - CCD id not found");
}

std::ostream& operator<<(std::ostream& os, const DrizzleSettings& ds)
{
    os <<"#Drizzle settings file\n"
       <<"#Version 1.0\n";

    os<< "Base Folder : "<<ds.output_base_folder <<"\n";

    os<< "CCDs present : ";
    for (size_t i=0; i!=ds.CCD_ids_present.size(); ++i) {
        os<<ds.CCD_ids_present[i]<<",";
    }
    os<<"\n";

    os<< "CCD output sub folders : ";
    for (size_t i=0; i!=ds.CCD_output_sub_folders.size(); ++i) {
        os<<ds.CCD_output_sub_folders[i]<<",";
    }
    os<<"\n";

    os <<"Normalisation on?  : " << string_utils::bool_to_string(ds.normalisation_on)<<"\n";
    os <<"Thresholding on? : " <<string_utils::bool_to_string(ds.thresholding_on)<<"\n";
    os <<"Perform post-drizzle col-debias?  : "<<string_utils::bool_to_string(
           ds.post_drizzle_col_debias_on)<<"\n";
    os <<"Create drizzle mosaics?  : "<<string_utils::bool_to_string(
           ds.create_drizzle_mosaic)<<"\n";

    os <<"Drizzle scale factor (typically 1.0, 0.5) : " <<ds.drizzle_scale_factor<<"\n";
    os <<"Drizzle pixfrac (typically 0.6) : " <<ds.drizzle_pixel_fraction<<"\n";

    os <<"Output image at cumulative percentiles : ";
    for (size_t i=0; i<ds.output_percentiles.size(); ++i) {
        os << ds.output_percentiles[i]<<",";
    }
    os <<"\n";
    return os;
}

void DrizzleSettings::parse_DrizzleSettings_text_vec(const vector<string>& file_text)
{
    using string_utils::parse_for_key_value;

    assert(!file_text.empty());
    output_base_folder = parse_for_key_value(file_text, "Base Folder");

    vector<string> ccd_id_strs = tokenize(parse_for_key_value(file_text,("CCDs present")),
                                          ",");
    for (size_t i=0; i!=ccd_id_strs.size(); ++i) {
        CCD_ids_present.push_back(atoi(ccd_id_strs[i]));
    }
    CCD_output_sub_folders = tokenize(parse_for_key_value(file_text,
                                      ("CCD output sub folders")),",");


    normalisation_on = string_utils::string_to_bool(parse_for_key_value(file_text,
                       "Normalisation on?"));
    thresholding_on =  string_utils::string_to_bool(parse_for_key_value(file_text,
                       "Thresholding on?"));
    post_drizzle_col_debias_on =  string_utils::string_to_bool(parse_for_key_value(file_text,
                                  "Perform post-drizzle col-debias?"));
    create_drizzle_mosaic = string_utils::string_to_bool(parse_for_key_value(file_text,
                            "Create drizzle mosaics?"));


    drizzle_scale_factor = atof(parse_for_key_value(file_text,("Drizzle scale factor")));
    drizzle_pixel_fraction = atof(parse_for_key_value(file_text,("Drizzle pixfrac")));



    vector<string> percentile_strings = tokenize(parse_for_key_value(file_text,
                                        ("Output image at cumulative percentiles")),",");
    for (size_t i=0; i!=percentile_strings.size(); ++i) {
        output_percentiles.push_back(atof(percentile_strings[i]));
    }

    std::sort(output_percentiles.begin(), output_percentiles.end());
}

}//end namespace coela
