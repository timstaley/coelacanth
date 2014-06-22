/*
 * File:   fits_header.h
 * Author: ts337
 *
 * Created on 07 February 2011, 14:18
 */

#ifndef COELA_FITS_HEADER_TABLE_H
#define COELA_FITS_HEADER_TABLE_H

#include "coela_utility/src/ordered_map.h"
#include "coela_utility/src/file_buffer.h"
#include "fits_header_conventions.h"

#include <istream>
#include <ostream>

namespace coela {

//fwd declaration
struct file_info;

//=======================================================================================
///Fits header table:
//Stores an ordered  map of "key/value/comments" (see includes)
//implements conventions of FITS header table formatting, and provides serialization in this form
//NB  will throw for overlong entries or disallowed characters in any field

class FitsHeader {
public:

    FitsHeader() {}

    ///load from stream
    FitsHeader(std::istream& is, const size_t hdr_begin_byte_offset) {Position_stream_and_load(is, hdr_begin_byte_offset);}
    ///load from file
    FitsHeader(const std::string& filename, const size_t hdr_begin_byte_offset=0);
    ///load from buffered file
    FitsHeader(const FileBuffer& , const size_t hdr_begin_byte_offset);

    void write_to_file(const std::string& filename);

    void merge_with_prefix_table(const FitsHeader& prefix);

    void clear() {entry_table.clear();}

    void add_comment(string comment);
    ///Throws if key is pre-existing
    void add_keyword(const string& keyword, const string& value="",
                     const string& key_comment="");
    void set_keyword(const string& keyword, const string& value);
    string get_key_value(const string& keyword) const;
    bool key_exists(const string& keyword) const {   return entry_table.key_exists(keyword);}
    void remove_key(const string& keyword) {entry_table.remove_key(keyword);}
    size_t size() const {return entry_table.size();}
    size_t header_file_length_in_bytes() const;

    fits_header_conventions::fits_imagetype FITS_imagetype() const;

    //Serialization:
    friend std::ostream& operator<<(std::ostream& , const FitsHeader&);
    friend std::istream& operator>>(std::istream& , FitsHeader &);

    bool operator==(const FitsHeader& rhs) const { return (entry_table == rhs.entry_table); }

private:
    void Position_stream_and_load(std::istream&, const size_t hdr_begin_byte_offset);
    OrderedMap entry_table;
    string make_keyword_line(const key_value_comment&) const;
    void load_from_string_vec(const vector<string>& header_text);
};


//=======================================================================================
}//end namespace coela
#endif  /* FITS_HEADER_TABLE_H */

