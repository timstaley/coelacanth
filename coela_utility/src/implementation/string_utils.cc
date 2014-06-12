
#include <vector>
#include <stdexcept>
//#include <iostream>
#include "../string_utils.h"
#include<iostream>

namespace coela {
using namespace std;
namespace string_utils {

bool string_to_bool(const std::string& s)
{
    if (s=="true") { return true; }
    return false;
}
std::string bool_to_string(const bool b)
{
    if (b) { return "true"; }
    return "false";
}

std::string pull_filename(const std::string& path)
{
    return tokenize(path, "/").back();
}

std::string pull_extension(const std::string& path)
{
    return path.substr(path.find_last_of('.'));
}
std::string strip_file_extension(const std::string& fname)
{
    return fname.substr(0, fname.find_last_of('.'));
}

string  strip_file_number_and_extension(const string& filename)
{
    return filename.substr(0, filename.find_last_of('_'));
}

std::string strip_lead_tail_chars(const std::string& input,
                                  const std::string& strip_chars)
{
    string copy(input);
    if (copy.find_first_not_of(strip_chars)==string::npos) { return ""; }
    copy = copy.substr(copy.find_first_not_of(strip_chars));
    copy = copy.substr(0, copy.find_last_not_of(strip_chars)+1);
    return copy;
}

string  pull_file_number_string(const string& filename)
{
    string filestem(strip_file_extension(filename));
    string number(filestem.substr(filestem.find_last_of('_') + 1));
    return number;
}

string  pull_camera_id_string(const string& filename)
{
    string filestem(strip_file_extension(filename));  //blah/[cam]/[subfolder]/run8009_0_1.lcz
    string filenum_stripped(filestem.substr(0,
                                            filestem.find_last_of('_')));     //contains /.../run8009_0
    string cam_id(filenum_stripped.substr(filenum_stripped.find_last_of('_') + 1));
    return cam_id;
}

float pull_file_number(const std::string& filename)
{
    string number_str = pull_file_number_string(filename);
    return atof(number_str);
}

std::string pull_filestem(const std::string& filename,
                          const bool greedy)
{
    string filestem = filename.substr(filename.find_last_of('/')+1);
    if (greedy == true){
      filestem = filestem.substr(0, filestem.find_last_of('_'));
    }else{
      filestem = filestem.substr(0, filestem.find('_'));
    }
    return filestem;
}

std::string find_and_replace(const std::string& input, const std::string& search_text,
                             const std::string& subst_text)
{
    string copy(input);
    string::size_type Position_marker;
    while ((Position_marker=copy.find(search_text))!=string::npos) {
        copy.erase(Position_marker, search_text.size());
        copy.insert(Position_marker, subst_text);
    }
    return copy;
}

std::string find_and_erase(const std::string& input, const std::string& erase_text)
{
    string copy(input);
    string::size_type Position_marker;
    while ((Position_marker=copy.find(erase_text))!=string::npos) {
        copy.erase(Position_marker, erase_text.size());
    }
    return copy;
}

int pull_arroyo_emitter_number(const std::string& filename)
{
    size_t linepos = filename.find("e-");
    linepos = filename.find('_', linepos); //find the first underscore after wavelength.
    string tail=filename.substr(linepos+1);
    tail = tail.substr(0, tail.find('_')); //snip off everything that isn't the emitter number
    return atoi(tail);
}

std::string parse_for_key_value(const std::vector<std::string>& text_lines,
                                const std::string& key)
{
//            string value;
    for (size_t i=0; i!=text_lines.size(); ++i) {
        const string& this_line = text_lines[i];
        if (this_line.find(key)!=string::npos) {
            return strip_lead_tail_chars(this_line.substr(this_line.find(':') + 1),
                                         " \t");    //strip any surrounding spaces or tabs
        }
    }
//            return value;
    throw runtime_error("Text Key not found in settings file: \""+key+"\" not found");
}

std::vector<std::string> tokenize(const std::string& str, const std::string& delimiters)
{
    // Discard delimiters at beginning of string.
    string::size_type start_token = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type end_token     = str.find_first_of(delimiters, start_token);
    vector<string> tokens;
    while (string::npos != end_token || string::npos != start_token) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(start_token, end_token - start_token));
        // Skip delimiters.  Note the "not_of"
        start_token = str.find_first_not_of(delimiters, end_token);
        // Find next "non-delimiter"
        end_token = str.find_first_of(delimiters, start_token);
    }
    return tokens;
}

std::vector<std::string> tokenize_and_strip_spaces(const std::string& str,
        const std::string& delimiters)
{
    vector<string> tokens = tokenize(str, delimiters);
    for (size_t i=0; i!=tokens.size(); ++i) {
        strip_lead_tail_chars(tokens[i]," ");
    }
    return tokens;
}




}//end namespace string_utils
}
