
#include "../file_buffer.h"
#include "../string_utils.h"

#include <fstream>
#include <stdexcept>
#include <iostream>
#include <cassert>
using namespace std;
namespace coela {



FileBuffer::FileBuffer(size_t initial_size)
{
    data.reserve(initial_size);
}


void FileBuffer::buffer_whole_file(const string& filename)
{
    std::ifstream infile(filename.c_str(), ios::binary | ios::ate);
    if (infile.is_open() && infile.tellg()!=-1) {
        size_t bytes_to_buffer=infile.tellg();
        infile.close();
        buffer_file_portion(filename, 0, bytes_to_buffer);
    } else { throw runtime_error(filename); }
}

void FileBuffer::buffer_file_portion(const string& filename, size_t start_offset,
                                     size_t bytes_to_buffer)
{
    recorded_filename="";
    try {
        std::ifstream infile(filename.c_str(),
                             ios::binary |
                             ios::ate);   //NB open file with stream Position set to end - this allows us to check that the file is long enough that it probably contains the data the header promises

        if (!infile.is_open() ||
                (infile.tellg() < (ptrdiff_t)(start_offset+bytes_to_buffer))) {
            throw runtime_error(filename);
        }

        data.resize(bytes_to_buffer);

        infile.seekg(start_offset);
        infile.read((&data[0]), bytes_to_buffer);
        infile.close();
        recorded_filename=filename;
    } catch (std::bad_alloc)
    {   throw runtime_error(" Could not allocate buffer of size " +  string_utils::itoa(bytes_to_buffer) + " elements");    }
    catch (std::runtime_error fname)
    {   throw  std::runtime_error(string("\nFatal Error buffering file: ") + fname.what()) ;   }
}


}//end namespace coela
