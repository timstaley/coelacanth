/*
 * File:   lucky_file_buffer.h
 * Author: ts337
 *
 * Created on 07 August 2009, 13:51
 */

#ifndef COELA_FILE_BUFFER_H
#define COELA_FILE_BUFFER_H

#include <string>
#include <vector>

using std::string;
///This stores information about where a frame is stored; whether the frame is a single file or part of an LCM.
namespace coela {
//=========================================================================================================

class FileBuffer {
public:
    FileBuffer(size_t initial_size=4*1024*1024); //4MB file size.

    void buffer_whole_file(const string&
                           filename); //finds length of file and calls "buffer_file_portion" with appropriate params
    void buffer_file_portion(const string& filename, size_t start_offset,
                             size_t bytes_to_buffer); //Does all the real work and record keeping

    const std::vector<char>& data_vec() const {return data;} //Currently only need read only access, although this might change
    string associated_filename() const {return recorded_filename;}


private:
    string recorded_filename;
    std::vector<char> data;
};
//=========================================================================================================
}//end namespace coela


#endif  /* _FILE_BUFFER_H */

