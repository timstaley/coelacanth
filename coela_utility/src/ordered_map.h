/*
 * File:   fits_header.h
 * Author: ts337
 *
 * Created on 06 November 2008, 13:24
 */

#ifndef COELA_ORDERED_MAP_H
#define COELA_ORDERED_MAP_H
#include <map>
#include <vector>
#include <string>

//I consider these to be so ubiquitous that a header "using" decl. is harmless
using std::string;
using std::vector;


namespace coela {

//=====================================================================================
//Simple struct to hang together a line containing key, value, and comment:
struct key_value_comment {

    //Constructor for a key / value pair
    key_value_comment(const string& input_key,
                      const string& input_value,
                      const string& input_comment=""):
        key(input_key), value(input_value), comment(input_comment) {}

    //Constructor for a comment line
    key_value_comment(const string& input_comment):comment(input_comment) {}


    //Data:
    string key, value, comment;

    bool operator==(const key_value_comment& rhs) const {
        return (key==rhs.key && value==rhs.value && comment == rhs.comment);
    }
};

//=====================================================================================


//=====================================================================================
///Ordered map:  (actually it's OrderedMap<string,string> - i.e. we could easily templatize it, but have no need currently.)
///the best way to get at key/values is through a map -
//but in this case (fits_header table) we want to keep the keys in a particular (Custom, non-alphabetical) order,
//so we actually store them in a vector and provide a key-linked index number for each.
//NB therefore "remove" operations are very inefficient, fortunately uncommon for FITS header purposes.

///This is actually a fairly generalised "key/value" map, which also happens to record an ordering

///It has one quirk- in general there can only be one copy of each key,
///However the user is allowed any number of empty key entries, to record comments.

///NB: Calling / derived class must enforce any desired characteristics
// e.g.
// * character limits on the length of text,
// * lack of special characters e.g. newlines.


class OrderedMap {
public:

    OrderedMap() {table.reserve(100);}
    void clear();
    size_t size() const {return table.size();}

    void add_key(const key_value_comment&); ///<NB throws if the key already exists
    void remove_key(const string& key); ///< If key does not exist, no side effects.
    bool key_exists(const string& key) const;

    void merge_with_prefix_map(const OrderedMap& prefix_map);

    ///Accessors via index key string:
    string & operator[](const string&
                        key); //If the key does not exist, it gets added and a reference to the empty value string is returned
    const string& operator[](const string& key)
    const; //Const version - if the key does not exist, throw an exception

    ///Accessor via ordering number:
    const key_value_comment operator[](const size_t line_number) const {return table[line_number];}

    bool operator==(const OrderedMap& rhs) const {
        return (table == rhs.table);
    }

private:
//----------------------------------------------------------------
    //Data:
    vector<key_value_comment> table;
    std::map <string, size_t> index_map;
//----------------------------------------------------------------

    void update_index_map();
};



//=====================================================================================
}//end namespace coela
#endif  /* _ORDERED_MAP_H */

