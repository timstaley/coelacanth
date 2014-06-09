/*
 * File:   simple_serialization_tools.h
 * Author: ts337
 *
 * Created on 20 June 2011, 12:30
 */

#ifndef COELA_SIMPLE_SERIALIZATION_TOOLS_H
#define COELA_SIMPLE_SERIALIZATION_TOOLS_H
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <typeinfo>
//#include <cstring>
#include <vector>
//#include <iostream>
#include "string_utils.h"

namespace coela {
namespace simple_serialization {

const std::string new_line ="\n";
const std::string key_value_separator = ":";
const std::string padded_key_value_separator = "  :  ";
const size_t max_line_size = 2048;

template<typename T>
std::ostream& output_key_value_pair_to_stream(std::ostream& os,
        const std::string& key, const T& value)
{
    os << key << padded_key_value_separator << value << new_line;
    return os;
}


/// Generic version, extracts just the first token after the separator
/// (see also string version defined in .cc file)
template< typename T >
void get_keyed_value_from_string_vec(const std::vector<std::string>& svec,
                                     const std::string& key, T& value)
{
    using namespace std;
    for (size_t i=0; i!=svec.size(); ++i) {
        if (svec[i].find(key) != string::npos) {
            std::stringstream ss(svec[i].substr(svec[i].find(key_value_separator) + 1));
            ss >> value;
            return;
        }
    }
    throw std::runtime_error("get_keyed_value_from_stream error: Key not found:" + key);
}


/// This is partially an error handling function,
/// e.g. in case we accidentally feed a large binary file to a function expecting a short text file
/// We want fast error reporting (and no buffer overflows)
/// so we effectively limit the read size to  max_lines * max_line_size.

/// It is also used to aid flexible serialization -
/// by grabbing a vector of string, we can loop over the stream contents as many times as we like
/// Without worrying about the nature of the original stream.

std::vector<std::string> convert_istream_to_string_vec(std::istream& is,
        const size_t max_lines_to_read=1000);

std::vector<std::string> convert_relevant_istream_section(std::istream& is,
        const std::string& start_key, const std::string& end_key,
        const size_t max_lines_to_read=1000);


template<class T>
void generic_output_to_file_stream(const std::string& filename, const T& object)
{
    std::ofstream outfile(filename.c_str(), std::ios::binary);
    //Requires typeinfo:
    if (!outfile.is_open())
        throw std::runtime_error(std::string("Error writing class of type ")+
                                 typeid(object).name() + " to file" + filename+" , could not open file.");
    //Simpler version
//    if (!outfile.is_open()) throw std::runtime_error("Error writing class to file" + filename);
    try {
        outfile<<object;
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("Error writing class of type ")+ typeid(
                                     object).name()
                                 + " to file \"" + filename+"\", caught exception which reads:\n"
                                 + e.what());
    }
    outfile.close();

}

template<class T>
void generic_load_from_file_stream(const std::string& filename, T& object)
{
//    T& object = *object_ptr;
    std::ifstream infile(filename.c_str(), std::ios::binary);
    //Requires typeinfo:
    if (!infile.is_open())
        throw std::runtime_error(std::string("Error loading class of type ")+
                                 typeid(object).name() + " from file" + filename+" , could not open file.");
    //Simpler version
//    if (!outfile.is_open()) throw std::runtime_error("Error writing class to file" + filename);
    try {
        infile>> (object);
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("Error loading class of type ")+ typeid(
                                     object).name()
                                 + " from file \"" + filename+"\", caught exception which reads:\n"
                                 + e.what());
    }
    infile.close();

}



//=================================================================================
} //End of namespace coela::simple_serialization
}


#endif  /* SIMPLE_SERIALIZATION_TOOLS_H */

