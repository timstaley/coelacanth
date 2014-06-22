/*
 * File:   settings.h
 * Author: ts337
 *
 * Created on 01 September 2008, 15:07
 */

#ifndef COELA_SETTINGS_H
#define COELA_SETTINGS_H
#include "coela_core/src/cartesian_coords.h"
#include "ccd_dataset_info.h"
#include "frame_info.h"

#include <string>
#include <vector>
using std::string;
using std::vector;

namespace coela {

//============================================================================================================
///A CCD_region with a CCD ID tacked on:
struct MultiCcdRegion: public CcdBoxRegion {
    MultiCcdRegion() {}
    MultiCcdRegion(const int cam_id, const CcdBoxRegion& rgn): CcdBoxRegion(rgn),
        CCD_id(cam_id) {}
    int CCD_id;
    ///Erases the regions it pulls out
    static vector<CcdBoxRegion> pull_regions_for_cam_id(int cam_id,
            vector<MultiCcdRegion>& rgns);
    static void push_regions_with_cam_id(vector<MultiCcdRegion>& camera_region_list,
                                         const vector<CcdBoxRegion>& ccd_rgns ,int cam_id);
};
std::ostream& operator<<(std::ostream& os, const MultiCcdRegion& c_rgn);





struct DrizzleSettings { ///< Basically stores everything parsed from the lucky_settings.txt
    DrizzleSettings();
    DrizzleSettings(const string& lucky_settings_filename);
    void write_to_file(const string & filename);


    string output_base_folder;

    vector<int> CCD_ids_present;
    vector<string> CCD_output_sub_folders;

    bool normalisation_on;
    bool thresholding_on;
    bool post_drizzle_col_debias_on;
    bool create_drizzle_mosaic;

    double drizzle_scale_factor, drizzle_pixel_fraction;

    vector<double> output_percentiles;

    string create_output_dirs(
        const string& gs_frames_list_filename,
        const vector<FrameInfo>& gs_frame_list,
        const vector< MultiFrame >& multi_frame_list
    ) const;

    string get_CCD_output_folder_for_CCD_id(int ccd_id) const;

private:
    void parse_DrizzleSettings_text_vec(const vector<string>& file_text);
    friend std::ostream& operator<<(std::ostream& os, const DrizzleSettings& ds);

};

}//end namespace coela

#endif  /* _SETTINGS_H */

