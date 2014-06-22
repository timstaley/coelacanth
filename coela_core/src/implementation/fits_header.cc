/*
 * File:   FitsHeader.cc
 * Author: ts337
 *
 * Created on 07 February 2011, 14:18
 */

//#include "coela_utility/src/my_debug_config.h"
#include "coela_core/src/fits_header.h"

#include "coela_utility/src/string_utils.h"
#include "coela_core/src/fits_header_conventions.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <iostream> //debug only

using std::runtime_error;

namespace coela {

//=================================================================================================
///Fits header table:

FitsHeader::FitsHeader(const std::string& filename, const size_t byte_offset)
{
    std::ifstream file_stream(filename.c_str() , std::ios::binary);
    if (!file_stream.is_open()) { throw runtime_error("FitsHeader - error loading from file: "+filename); }
    Position_stream_and_load(file_stream, byte_offset);
    file_stream.close();
}

void FitsHeader::Position_stream_and_load(std::istream& is, const size_t byte_offset)
{
    is.seekg(0, std::ios::end);
    if (byte_offset> (size_t)is.tellg() || is.tellg()==-1) {throw runtime_error("header stream is empty or truncated");}
    is.seekg(byte_offset, std::ios_base::beg);
    is >> *this;
}


FitsHeader::FitsHeader(const FileBuffer& buf, const size_t hdr_begin_byte_offset)
{
    ///This is a pain - we have to convert an array of raw data back to a stream somehow, to use the standard operator...
    //Go about it by grabbing just the header part, then using a stringstream.
    assert(!buf.data_vec().empty());
    assert(hdr_begin_byte_offset<buf.data_vec().size());

    vector<char>::const_iterator current_Position(buf.data_vec().begin() +
            hdr_begin_byte_offset);
    string line;
    while (line.find("END ")  != 0 &&   current_Position < buf.data_vec().end()) {
        line= string(current_Position, current_Position+ fits_header_conventions::line_max);
        current_Position+=fits_header_conventions::line_max;
    }
    if (current_Position>=buf.data_vec().end()) throw
        runtime_error("Problem loading header from buffer; apparently undersize; "+
                      string_utils::itoa(buf.data_vec().size()) +" < " + string_utils::itoa(
                          current_Position-buf.data_vec().begin()) +"\n"
                      +"Byte offset: "+string_utils::itoa(hdr_begin_byte_offset)
                      +"\nFilename: " + buf.associated_filename()
                      +"\n First line:\n"
                      +string(buf.data_vec().begin()+hdr_begin_byte_offset,
                              buf.data_vec().begin()+hdr_begin_byte_offset+ fits_header_conventions::line_max)
                     );
    std::stringstream buffer_stream(string(buf.data_vec().begin()+hdr_begin_byte_offset,
                                           current_Position));

    buffer_stream >> *this;
}

//FitsHeader::FitsHeader(const FileInfo& fi){
//    *this=FitsHeader(fi.file_path, fi.byte_offset);
//}

void FitsHeader::write_to_file(const std::string& filename)
{
    std::ofstream outfile(filename.c_str(), std::ios::binary);
    if (!outfile.is_open()) {
      throw runtime_error(
          "FitsHeader::save_to_file - error writing to file: "+ filename);
    }

    outfile << *this;

    outfile.close();
}

void FitsHeader::merge_with_prefix_table(const FitsHeader& prefix)
{
    entry_table.merge_with_prefix_map(prefix.entry_table);
}

void FitsHeader::add_comment(string comment_input)
{
    //Need to split input into lines of 80 - sizeof("COMMENT ")length.
    const size_t comment_max= fits_header_conventions::line_max -
                              fits_header_conventions::comment_prefix.size();
    const string& prefix = fits_header_conventions::comment_prefix;
    while (comment_input.size() >comment_max) {
        entry_table.add_key(KeyValueComment(prefix+comment_input.substr(0, comment_max)));
        comment_input.erase(0,comment_max);
    }
    entry_table.add_key(KeyValueComment(prefix+comment_input));
}

void FitsHeader::add_keyword(const string& keyword, const string& value,
                             const string& key_comment)
{
    if (keyword.size() > fits_header_conventions::key_max
            || value.size() > fits_header_conventions::value_max
            || key_comment.size() > fits_header_conventions::kv_comment_max) {
        throw std::runtime_error("FitsHeader::add_keyword - One of the table entries is too long for the FITS format");
    }
    entry_table.add_key(KeyValueComment(keyword,value, key_comment));
}

void FitsHeader::set_keyword(const string& keyword, const string& value)
{
    if (keyword.size() > fits_header_conventions::key_max
            || value.size() > fits_header_conventions::value_max) {
        throw std::runtime_error("FitsHeader::set_keyword - One of the table entries is too long for the FITS format");
    }
    entry_table[keyword]=value;
}

string FitsHeader::get_key_value(const string& keyword) const
{
    if (entry_table.key_exists(keyword)) {
        return entry_table[keyword];
    } else { throw std::runtime_error("FitsHeader::get_key_value -  Key does not exist: "+ keyword); }
}




size_t FitsHeader::header_file_length_in_bytes() const
{
    //+1 for END line not stored. then /36 lines per card, +1 since remainder goes to a whole card. 2880 bytes per card
    int rpc = fits_header_conventions::rows_per_card;
    if ((entry_table.size() +1) % rpc ==0) { return ((entry_table.size() +1) /rpc) *fits_header_conventions::bytes_per_card; }
    else {
        return ((entry_table.size() +1) /rpc  + 1) * fits_header_conventions::bytes_per_card;
    }
}

fits_header_conventions::fits_imagetype FitsHeader::FITS_imagetype() const
{
    return (fits_header_conventions::fits_imagetype)string_utils::atoi(
               get_key_value("BITPIX"));
}

string FitsHeader::make_keyword_line(const KeyValueComment& entry) const
{
    string line;
    if (!entry.key.empty()) {
        line = entry.key;
        //pad to col 9
        while (line.size()<fits_header_conventions::key_max) { line +=' '; }
        line += "= ";  //cols 9 & 10
        //right justify value to col 30
        while (line.size() + entry.value.size()
                < fits_header_conventions::value_col_max) {
            line += ' ';
        }
        line += entry.value;
        if (!entry.comment.empty()) {
            line +=" / ";
            line +=entry.comment;
        }
    } else {
        line = entry.comment;
    }

    line=line.substr(0,
                     fits_header_conventions::line_max); //chop off any excess chars //NB shouldn't need this anymore but may as well leave it in.
    //pad to width 80 if underlength;
    while (line.size() != fits_header_conventions::line_max) { line +=' ' ; }
    return line;
}


std::ostream& operator<<(std::ostream& os, const FitsHeader& hdr)
{
    for (size_t i=0; i!=hdr.entry_table.size(); ++i) {
        os << hdr.make_keyword_line(hdr.entry_table[i]);
    }
    size_t lines_in_header=hdr.entry_table.size();
    string line = "END";
    //pad to width 80;
    while (line.size() != fits_header_conventions::line_max) { line +=' ' ; }
    os << line;
    lines_in_header++;

    //Length of a FITS header must be 36 lines of 80 cols, need to pad
    string spaces(fits_header_conventions::line_max, ' ');  //a string of 80 spaces
    //Each header chunk must be 36 lines long (80 chars per line):
    while (lines_in_header % fits_header_conventions::rows_per_card  != 0) {
        os << spaces; lines_in_header++;
    }
    return os;
}

std::istream& operator>>(std::istream& is, FitsHeader & hdr)
{
    hdr.clear();
    vector<string> header_text;
    char line_chars[85];
//    is.seekg(byte_offset, ios_base::beg); //NB must implement this under calling code !!!

    //First, split the stream into lines of length 80 chars, these will be easier to process
    string line;
    unsigned int line_max = fits_header_conventions::line_max;

    //Run a check on the first line to make sure we have the right kind of file:
    is.get(line_chars, line_max
           +1);  //gets 'up to but not including' the Nth char, apparently.
    line=line_chars;
    header_text.push_back(line);
    if (line.find(fits_header_conventions::standard_first_ten_chars) != 0) {
        throw runtime_error("File does not conform to implemented FITS header standard; got:\n"
                            +line+"\n");
    }

    do {
        is.get(line_chars, line_max
               +1);  //gets 'up to but not including' the Nth char, apparently.
        line=line_chars;
//        assert(line.size()==line_max);
        header_text.push_back(line);
        //NB the END keyword comes on a line of it's own:
    } while ((line.find("END ")!=0) && !is.bad());
    if (is.bad()) { throw std::runtime_error("std::istream& operator>>(std::istream& is, FitsHeader &) - Problem loading header: istream not returning nicely"); }
    hdr.load_from_string_vec(header_text);
    return is;
}



void FitsHeader::load_from_string_vec(const vector<string>& header_text)
{
    string::size_type line_Position;

    for (vector<string>::size_type l=0; l<header_text.size() ; ++l) {
        string line=header_text[l];  //simply for readability.
        KeyValueComment key_to_be_added("");
        string keyword, value;


        //Parse for COMMENTs first since it's possible to have an equals operator in a comment:
        if (line.find("COMMENT ") == 0) {
            key_to_be_added=KeyValueComment(line) ;
//          if (fits_hdr_debug) cout <<"comment loaded" <<endl;
            entry_table.add_key(key_to_be_added);

        } else if (line.find("HIERARCH ") == 0) {
            //Special case; VLT/NACO images have multiple COMMENT like entries prefaced "HIERARCH":
            key_to_be_added=KeyValueComment(line) ;
//          if (fits_hdr_debug) cout <<"comment loaded" <<endl;
            entry_table.add_key(key_to_be_added);

        } else if ((line_Position = line.find('=')) != string::npos)  {
            //(if it's a line with a keyword / value pair)
//            assert(line.find("=")==key_max);
            line_Position++; //move past '=' to string of spaces...
            line_Position = line.find_first_not_of(' ', line_Position); //move to start of value

            if (line_Position != string::npos) {
                //(if it's not blank!)
                value = line.substr(line_Position);// the rest of the line after and including the value
            } else { value =" "; }

            keyword = line.substr(0, line.find('='));
            keyword = keyword.substr(0, keyword.find(' '));

            string key_comment;
            //comment comes after the forward slash and a space:
            if (value.find('/')!=string::npos) { key_comment  = value.substr(value.find('/') +2); }
            value=value.substr(0, value.find('/'));
//              if (value[value.size()-1]==" ") value=value.substr(0,value.size()-1);
            if (value.find_last_not_of(' ')+1 < value.size()) { value=value.substr(0, value.find_last_not_of(' ')+1); }
            key_to_be_added=KeyValueComment(keyword, value, key_comment) ;
            //NB END keyword not stored - this is appended in write to file.
            entry_table.add_key(key_to_be_added);
        }



    }
//    if (fits_hdr_debug) cout << "header parsed" <<endl;
}


///end Fits header table.
//=======================================================================================



}//end namespace coela

