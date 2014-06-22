#include <boost/python.hpp>
#include <boost/python/data_members.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../ccd_dataset_info.h"
#include "../drizzle_settings.h"


BOOST_PYTHON_MODULE(lucky_pipeline)
{
    using namespace boost::python;
    using namespace coela;

    class_<CcdDatasetInfo>("CcdDatasetInfo")
    .def(init<std::string>()) //Load from file(filename)
    .def("write", &CcdDatasetInfo::write_to_file)
    .def_readwrite("ccd_id",
                   &CcdDatasetInfo::ccd_id)
    .def_readwrite("CCD_inputdir",
                   &CcdDatasetInfo::CCD_inputdir)
    .def_readwrite("CCD_outputdir",
                   &CcdDatasetInfo::CCD_outputdir)
    .def_readwrite("in_filestem",
                   &CcdDatasetInfo::in_filestem)
    .def_readwrite("out_filestem",
                   &CcdDatasetInfo::out_filestem)
    .def_readwrite("extension",
                   &CcdDatasetInfo::extension)
    .def_readwrite("dataset_output_base_dir",
                   &CcdDatasetInfo::dataset_output_base_dir)

    .def_readwrite("default_camera_config_file",
                   &CcdDatasetInfo::default_camera_config_file)
    .def_readwrite("default_guide_star_region_file",
                   &CcdDatasetInfo::default_guide_star_region_file)
    .def_readwrite("default_histogram_region_file",
                   &CcdDatasetInfo::default_faint_histogram_region_file)
    .def_readwrite("most_recently_output_frame_list",
                   &CcdDatasetInfo::most_recently_output_frame_list)
    ;

    class_< std::vector<double> >("double_vec")
    .def(vector_indexing_suite< std::vector<double> >());

    class_< std::vector<int> >("int_vec")
    .def(vector_indexing_suite< std::vector<int> >());

    class_< std::vector<std::string> >("string_vec")
    .def(vector_indexing_suite< std::vector<std::string> >());


    class_<DrizzleSettings>("DrizzleSettings")
    .def(init<std::string>()) //Load from file(filename)
    .def("write", &DrizzleSettings::write_to_file)
    .def_readwrite("output_base_folder",
                   &DrizzleSettings::output_base_folder)
    .def_readwrite("CCD_ids_present",
                   &DrizzleSettings::CCD_ids_present)
    .def_readwrite("CCD_output_sub_folders",
                   &DrizzleSettings::CCD_output_sub_folders)
    .def_readwrite("normalisation_on",
                   &DrizzleSettings::normalisation_on)
    .def_readwrite("thresholding_on",
                   &DrizzleSettings::thresholding_on)
    .def_readwrite("post_drizzle_col_debias_on",
                   &DrizzleSettings::post_drizzle_col_debias_on)
    .def_readwrite("create_drizzle_mosaic",
                   &DrizzleSettings::create_drizzle_mosaic)
    .def_readwrite("drizzle_scale_factor",
                   &DrizzleSettings::drizzle_scale_factor)
    .def_readwrite("drizzle_pixel_fraction",
                   &DrizzleSettings::drizzle_pixel_fraction)
    .def_readwrite("output_percentiles",
                   &DrizzleSettings::output_percentiles)
    ;






}
