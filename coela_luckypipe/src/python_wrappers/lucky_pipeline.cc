#include <boost/python.hpp>
#include <boost/python/data_members.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../ccd_dataset_info.h"
#include "../drizzle_settings.h"


BOOST_PYTHON_MODULE(lucky_pipeline)
{
    using namespace boost::python;
    using namespace coela;

    class_<CCD_DatasetInfo>("CCD_DatasetInfo")
    .def(init<std::string>()) //Load from file(filename)
    .def("write", &CCD_DatasetInfo::write_to_file)
    .def_readwrite("ccd_id",
                   &CCD_DatasetInfo::ccd_id)
    .def_readwrite("CCD_inputdir",
                   &CCD_DatasetInfo::CCD_inputdir)
    .def_readwrite("CCD_outputdir",
                   &CCD_DatasetInfo::CCD_outputdir)
    .def_readwrite("in_filestem",
                   &CCD_DatasetInfo::in_filestem)
    .def_readwrite("out_filestem",
                   &CCD_DatasetInfo::out_filestem)
    .def_readwrite("extension",
                   &CCD_DatasetInfo::extension)
    .def_readwrite("dataset_output_base_dir",
                   &CCD_DatasetInfo::dataset_output_base_dir)

    .def_readwrite("default_camera_config_file",
                   &CCD_DatasetInfo::default_camera_config_file)
    .def_readwrite("default_guide_star_region_file",
                   &CCD_DatasetInfo::default_guide_star_region_file)
    .def_readwrite("default_histogram_region_file",
                   &CCD_DatasetInfo::default_faint_histogram_region_file)
    .def_readwrite("most_recently_output_frame_list",
                   &CCD_DatasetInfo::most_recently_output_frame_list)
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
