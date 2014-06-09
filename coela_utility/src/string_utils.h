/*
 * File:   string_utils.h
 * Author: ts337
 *
 * Created on 24 November 2008, 17:00
 */

#ifndef COELA_STRING_UTILS_H
#define COELA_STRING_UTILS_H
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <iomanip>

namespace coela {
namespace string_utils {

inline int atoi(const std::string& a) {return std::atoi(a.c_str());} //returns 0 if a is empty.
inline long atol(const std::string& a) {return std::atol(a.c_str());} //returns 0 if a is empty.
inline double atof(const std::string& a) {return std::atof(a.c_str());}

template<class T>
std::string num_to_string(const T&);
std::string itoa(long int x);
inline std::string itoa(int x) {return itoa((long int)x);}
inline std::string itoa(size_t x) {return itoa(static_cast<long int>(x));}
std::string itoa(float x)
; //undefined, declared in an attempt to prevent cock-ups at compile time (better than finding them at run-time!).
std::string ftoa(double x, size_t precision=12)
; //no harm in promoting an int to a float here, so ftoa(int) is left implicitly defined.
bool string_to_bool(const std::string&);
std::string bool_to_string(const bool);


std::string pull_filename(const std::string&
                          path); //returns only the part of the file path below any containing directories
std::string pull_extension(const std::string& path);
std::string strip_file_extension(const std::string&
                                 filename); // removes everything after and including the final '.'
std::string strip_file_number_and_extension(const std::string& filename);
std::string pull_filestem(const std::string& filename);
std::string pull_file_number_string(const std::string& filename);
std::string strip_lead_tail_chars(const std::string& input,
                                  const std::string& strip_chars);

std::string  pull_camera_id_string(const std::string& filename);

///tokenize: strips any chars from the delimiters input from the beginning, middle and end of a string, returning a vector of the bits inbetween
std::vector<std::string> tokenize(const std::string& str,
                                  const std::string& delimiters=" ");
///Same as above, but then returns the substrings with any spaces stripped.
std::vector<std::string> tokenize_and_strip_spaces(const std::string& str,
        const std::string& delimiters=" ");

std::string find_and_replace(const std::string& input, const std::string& search_text,
                             const std::string& subst_text);
std::string find_and_erase(const std::string& input, const std::string& erase_text);

//FIXME - this duplicates functionality in settings.cc
std::string parse_for_key_value(const std::vector<std::string>& text_lines,
                                const std::string& key);

int pull_arroyo_emitter_number(const std::string& filename);
float pull_file_number(const std::string& filename);

//while atoi() is defined in stdlib, we have to do this one ourselves!
inline std::string ftoa(double x, size_t precision)
{
    std::stringstream ss;
    ss<<std::setprecision(precision);
    ss<<x;
    return ss.str();
}

inline std::string itoa(long int
                        x)   //while atoi() is defined in stdlib, we have to do this one ourselves!
{
    std::stringstream ss;
    ss<<x;
    return ss.str();
}

template<class T>
std::string num_to_string(const T& num)
{
    std::stringstream ss;
    ss<<num;
    return ss.str();
}

} //end namespace coela::string_utils
}//end namespace coela

#endif  /* _STRING_UTILS_H */

