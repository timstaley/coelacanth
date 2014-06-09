/*
 * File:   histogram_container.h
 * Author: ts337
 *
 * Created on 30 September 2009, 16:40
 */


#ifndef COELA_HISTOGRAM_CONTAINER_H
#define COELA_HISTOGRAM_CONTAINER_H
#include <vector>
#include <map>
using std::vector;

namespace coela {

//TO DO - this class needs re-writing.

///Wraps a vector so that essentially we use it as an array of pairs - the indices corresponds to a histogram value, the indexed memory corresponds to the counts
struct HistogramContainer14bit: public vector<long> {
    HistogramContainer14bit(): vector<long>(fourteen_bit_range+1, 0), pixels_counted(0),
        zero_index_value(0), max_index_value(fourteen_bit_range), min_val_set(false) {}

    void write_to_file(const std::string& filename) const;

    void count_int_value(int ivalue);
    void count_rounded_value(double fvalue);
    void pedantic_count_value(int
                              ivalue); ///<Will throw an exception if the given value cannot be stored in the histogram
    static const int fourteen_bit_range = 16384;
    long pixel_count() const {return pixels_counted;}
    int mode() const;
    void clear();
    void set_min_value(int);
    int storage_range_min_value()const {return zero_index_value;}
    int storage_range_max_value()const {return max_index_value;}
    int minimum_has_been_set() const {return min_val_set;}
    long& get_count(int value);
    const long& get_count(int value)const;
    long get_value(int index) const {return index+storage_range_min_value();}

    std::map<int,long> convert_to_value_count_map() const;


private:
//        int mode_ignoring_zero() const;
    long pixels_counted;
    int zero_index_value;
    int max_index_value;
    bool min_val_set;
};

//================================================================================
struct HistogramContainer10bit: public vector<long> {
    HistogramContainer10bit(): vector<long>(ten_bit_range+1, 0), pixels_counted(0),
        zero_index_value(0), max_index_value(ten_bit_range), min_val_set(false) {}

    void write_to_file(const std::string& filename) const;

    void count_int_value(int ivalue);
    void count_rounded_value(double fvalue);
    void pedantic_count_value(int
                              ivalue); ///<Will throw an exception if the given value cannot be stored in the histogram
    static const int ten_bit_range = 1024;
    long pixel_count() const {return pixels_counted;}
    int mode() const;
    void clear();
    void set_min_value(int);
    int storage_range_min_value()const {return zero_index_value;}
    int storage_range_max_value()const {return max_index_value;}
    int minimum_has_been_set() const {return min_val_set;}
    long& get_count(int value);
    const long& get_count(int value)const;
    long get_value(int index) const {return index+storage_range_min_value();}

    std::map<int,long> convert_to_value_count_map() const;


private:
//        int mode_ignoring_zero() const;
    long pixels_counted;
    int zero_index_value;
    int max_index_value;
    bool min_val_set;
};

}//end namespace coela


#endif  /* _HISTOGRAM_CONTAINER_H */

