/*
 * File:   file_utils.cc
 * Author: ts337
 *
 * Created on 09 February 2011, 14:32
 */

#include "../file_utils.h"
#include "../string_utils.h"

#include "../simple_serialization.h"
#include <boost/filesystem.hpp>
#include <fstream>

namespace bfs=boost::filesystem;
using namespace std;
namespace coela {
namespace file_utils {

std::list<string> get_matching_files_recursively(const std::string& base_dir,
        const std::string& extension_input)
{
    //NB Boost extension seems to return with the dot included.
    //whereas I prefer to store without the dot
    string extension=extension_input;
    if (!extension.empty() && extension[0]!='.') { extension= "."+ extension; }

    list<string> matching_filenames;
    if (!bfs::exists(base_dir) || !bfs::is_directory(base_dir)) { throw runtime_error("Cannot find dir: \""+base_dir+"\""); }
    bfs::recursive_directory_iterator end_iter;
    for (bfs::recursive_directory_iterator dir_itr(base_dir); dir_itr != end_iter;
            ++dir_itr) {
        if (bfs::is_regular_file(*dir_itr) && bfs::extension(*dir_itr)==extension) {
            matching_filenames.push_back(dir_itr->path().string());
        }
    }
    matching_filenames.sort(compare_filenames);
    return matching_filenames;
}



bool compare_filenames(const string& first, const string& second)
{
    string stem1,stem2;
    stem1 = string_utils::strip_file_number_and_extension(first);
    stem2 =string_utils::strip_file_number_and_extension(second);
    if (stem1==stem2) {
        float num1, num2;
        num1 = string_utils::pull_file_number(first);
        num2 = string_utils::pull_file_number(second);
        return (num1 < num2);
    } else { return (stem1 <stem2); }
}

std::vector<string> convert_text_file_to_string_vec(const std::string& filename)
{
    if (!boost::filesystem::is_regular_file(filename)) {
        throw runtime_error("Error loading text file: does not appear to be a regular file (filename: \""
                            +filename+"\")");
    }
    ifstream text_file(filename.c_str());
    if (text_file.is_open()) {
        return simple_serialization::convert_istream_to_string_vec(text_file);
    } else { throw runtime_error("File not found (or error opening):" +filename); }
}

std::string get_file_last_access_timestring(const std::string& filename)
{
    std::time_t file_time = boost::filesystem::last_write_time(filename);
    string timestr((std::asctime(std::gmtime(&file_time))));
    timestr=timestr.substr(0, timestr.find("\n"));
    return timestr;
}



}//end namespace coela::file_utils
}//end namespace coela
