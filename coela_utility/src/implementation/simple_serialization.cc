#include <limits.h>

#include "../simple_serialization.h"
#include "../string_utils.h"

#include <cstring>
#include <stdexcept>

using std::string;
using std::vector;
namespace coela {
namespace simple_serialization {

std::vector<std::string> convert_istream_to_string_vec(std::istream& is,
        const size_t max_lines_to_read)
{
    if (is) {
        char line_buffer[max_line_size];
        vector<string> text;
        size_t lines_read=0;
        while (is.getline(line_buffer, max_line_size) && lines_read!=max_lines_to_read) {
            text.push_back(line_buffer);
            if (!text.back().empty() && text.back()[text.back().size()-1]=='\r') {
                text.back().resize(text.back().size()-1);
            }
        }
        return text;
    } else { throw std::runtime_error("Could not convert istream"); }
}

std::vector<std::string> convert_relevant_istream_section(std::istream& is,
        const std::string& start_key, const std::string& end_key,
        const size_t max_lines_to_read)
{

    std::vector<string> relevant_lines;
    char line_buffer[max_line_size];

    size_t lines_read=0;
    bool lines_relevant=false;
    while (is.getline(line_buffer, max_line_size) && lines_read!=max_lines_to_read) {
        if (!lines_relevant) {
            if (strncmp(line_buffer, start_key.c_str() , start_key.size()) == 0) {
                lines_relevant=true;
            }
        } else if (strncmp(line_buffer, end_key.c_str() , end_key.size()) == 0) {
            return relevant_lines;
        } else {
            relevant_lines.push_back(line_buffer);
            string& last_line = relevant_lines.back();
            if (!last_line.empty() && last_line[last_line.size()-1]=='\r') {
                last_line.resize(last_line.size()-1);
            }
        }
    }
    if (!lines_relevant)
        throw std::runtime_error(
            "get_relevant_istream_section: start key not found:" + start_key);
    throw std::runtime_error("get_relevant_istream_section: end key not found:" + end_key);
}

///String version, returns all the line after the separator rather than just the first token.
template<>
void get_keyed_value_from_string_vec<std::string>(const std::vector<std::string>& svec,
        const std::string& key, std::string& value)
{
    using namespace std;
    for (size_t i=0; i!=svec.size(); ++i) {
        if (svec[i].find(key) != string::npos) {
            value =
                string_utils::strip_lead_tail_chars(
                    svec[i].substr(svec[i].find(key_value_separator) + 1) , " ");
            return;
        }
    }
    throw std::runtime_error("get_keyed_value_from_stream error: Key not found:" + key);
}


//=================================================================================
} //End of namespace coela::simple_serialization
}

