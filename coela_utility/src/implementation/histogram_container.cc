#include <stdexcept>
//#include "../level1/my_debug_config.h"
#include "../histogram_container.h"
#include <cassert>
#include <fstream>
#include <cmath>
using namespace std;
namespace coela {

int fast_rint_alternative(const double f)
{
    int i = (int) f;
    if (f >= 0.0) {
        return ((f-i) >= 0.5) ? (i + 1) : (i);
    } else {
        return (-f+i >= 0.5) ? (i - 1) : (i);
    }
}

void HistogramContainer14bit::set_min_value(int min_value_)
{
    clear();
    zero_index_value=min_value_;
    max_index_value = zero_index_value + size()-1; //The value stored at [size()-1]
    min_val_set=true;
}

void HistogramContainer14bit::count_int_value(int value)
{
    //simply drop the value if out of range:
    if (value<storage_range_min_value() || value > storage_range_max_value()) { return; }
    ++pixels_counted;
    ((*this)[value - storage_range_min_value()])++;
}

void HistogramContainer14bit::count_rounded_value(double fvalue)
{
//        count_value((int)rint(fvalue));
    count_int_value(fast_rint_alternative((fvalue)));
}

void HistogramContainer14bit::pedantic_count_value(int value)
{
//        if (value<min_hist_val() || value > max_hist_val() ) return;
    if (value<storage_range_min_value() || value > storage_range_max_value()) { throw runtime_error("Out of range value for this histogram container"); }
    count_int_value(value);
}

long& HistogramContainer14bit::get_count(int value)
{
    assert(value>=storage_range_min_value() && value<=storage_range_max_value());
    return (*this)[value - storage_range_min_value()];
}
const long& HistogramContainer14bit::get_count(int value) const
{
    assert(value>=storage_range_min_value() && value<=storage_range_max_value());
    return (*this)[value - storage_range_min_value()];
}

int HistogramContainer14bit::mode() const
{
    int mode(0), max_count(0);
    for (size_t i=0; i<size(); ++i) {
        if ((*this)[i]>max_count) {
            mode= get_value(i);
            max_count=(*this)[i];
        }
    }
    return mode;
}

void HistogramContainer14bit::clear()
{
    pixels_counted=0;
    for (size_t i=0; i<size(); ++i) {
        (*this)[i]=0;
    }
}

map<int,long> HistogramContainer14bit::convert_to_value_count_map() const
{
    map<int,long> VC_map;
    for (size_t i=0; i!=this->size(); ++i) {
        if ((*this)[i]) { VC_map[get_value(i)]=(*this)[i]; }
    }
    return VC_map;
}


void HistogramContainer14bit::write_to_file(const string& filename) const
{
    std::ofstream datafile(filename.c_str());
    for (size_t i=0; i!=size(); ++i) {
        if ((*this)[i]) { datafile << (*this).get_value(i) <<"  " <<(*this)[i] <<"\n"; }
    }
    datafile.close();
}
//================================================================================
void HistogramContainer10bit::set_min_value(int min_value_)
{
    clear();
    zero_index_value=min_value_;
    max_index_value = zero_index_value + size()-1; //The value stored at [size()-1]
    min_val_set=true;
}

void HistogramContainer10bit::count_int_value(int value)
{
    //simply drop the value if out of range:
    if (value<storage_range_min_value() || value > storage_range_max_value()) { return; }
    ++pixels_counted;
    ((*this)[value - storage_range_min_value()])++;
}

void HistogramContainer10bit::count_rounded_value(double fvalue)
{
//        count_value((int)rint(fvalue));
    count_int_value(fast_rint_alternative((fvalue)));
}

void HistogramContainer10bit::pedantic_count_value(int value)
{
//        if (value<min_hist_val() || value > max_hist_val() ) return;
    if (value<storage_range_min_value() || value > storage_range_max_value()) { throw runtime_error("Out of range value for this histogram container"); }
    count_int_value(value);
}

long& HistogramContainer10bit::get_count(int value)
{
    assert(value>=storage_range_min_value() && value<=storage_range_max_value());
    return (*this)[value - storage_range_min_value()];
}
const long& HistogramContainer10bit::get_count(int value) const
{
    assert(value>=storage_range_min_value() && value<=storage_range_max_value());
    return (*this)[value - storage_range_min_value()];
}

int HistogramContainer10bit::mode() const
{
    int mode(0), max_count(0);
    for (size_t i=0; i<size(); ++i) {
        if ((*this)[i]>max_count) {
            mode= get_value(i);
            max_count=(*this)[i];
        }
    }
    return mode;
}

void HistogramContainer10bit::clear()
{
    pixels_counted=0;
    for (size_t i=0; i<size(); ++i) {
        (*this)[i]=0;
    }
}

map<int,long> HistogramContainer10bit::convert_to_value_count_map() const
{
    map<int,long> VC_map;
    for (size_t i=0; i!=this->size(); ++i) {
        if ((*this)[i]) { VC_map[get_value(i)]=(*this)[i]; }
    }
    return VC_map;
}


void HistogramContainer10bit::write_to_file(const string& filename) const
{
    std::ofstream datafile(filename.c_str());
    for (size_t i=0; i!=size(); ++i) {
        if ((*this)[i]) { datafile << (*this).get_value(i) <<"  " <<(*this)[i] <<"\n"; }
    }
    datafile.close();
}

}//end NS lucky.
