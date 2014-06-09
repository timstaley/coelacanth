/*
 * File:   file_utils.h
 * Author: ts337
 *
 * Created on 09 February 2011, 14:32
 */

#ifndef COELA_FILE_UTILS_H
#define COELA_FILE_UTILS_H
#include <list>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <typeinfo>

namespace coela {
namespace file_utils {

std::list<std::string> get_matching_files_recursively(const std::string& base_dir,
        const std::string& extension="");
bool compare_filenames(const std::string& first, const std::string& second);

std::string get_file_last_access_timestring(const std::string& filename);

//Deprecated (to do - remove!)
std::vector<std::string> convert_text_file_to_string_vec(const std::string& filename);


}//end namespace coela::file_utils
}//end namespace coela


#endif  /* FILE_UTILS_H */

