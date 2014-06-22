#include "../ds9_region.h"
#include "coela_utility/src/string_utils.h"

#include <fstream>
#include <stdexcept>

#include <algorithm>
#include <iostream>


using namespace std;

namespace coela {
namespace ds9 {

//========================================================================================

DS9Region::DS9Region(const std::string& type, const std::string& ref_frame,
                     const std::string& sub_type):
    region_type(type), coord_ref_frame(ref_frame), sub_type_comment(sub_type) {}

std::ostream& operator<<(std::ostream& os, const DS9Region& r)
{
    os<<r.region_type<<"(";

    for (size_t j=0; j!=r.parameters.size(); ++j) {
        os << r.parameters[j];
        if (j!=r.parameters.size()-1) { os<<","; }
    }

    os<<")" << (r.sub_type_comment.empty() ? "" : " ")
      <<r.sub_type_comment;
    return os;
}



std::istream& operator>>(std::istream& is, DS9Region& r)
{
    r.parameters.clear();


    std::stringbuf buf;
    is.get(buf,'(');
    r.region_type=buf.str();
    char c=is.get();
    assert(c=='(');
    do {
        double p;
        is>>p;
        r.parameters.push_back(p);
        is>>c;
    } while (c==',');
    assert(c==')');

    if (is.peek()==' ') {
        is.get();
        buf.str("");
        is.get(buf);
        r.sub_type_comment=buf.str();
    }


    return is;
}
//-----------------------------------------------------------------------------------------------
template<class T>
DS9Region DS9Region::point(const Position<T>& posn, const std::string marker_type_string)
{
    //to do - assert marker type string is a valid type
    DS9Region r("point", posn.ds9_reference_frame_string(), "# point="+marker_type_string);
    r.parameters.reserve(2);
    r.parameters.push_back(posn.x);
    r.parameters.push_back(posn.y);
    return r;
}

template DS9Region DS9Region::point(const Position<coordinate_types::pixels>& ,
                                    const std::string);
template DS9Region DS9Region::point(const Position<coordinate_types::CCD>& ,
                                    const std::string);
template DS9Region DS9Region::point(const Position<coordinate_types::mosaic>&,
                                    const std::string);

template<class T>
DS9Region DS9Region::box(const RectangularRegion<T>& rgn)
{
    DS9Region r("box", rgn.ds9_reference_frame_string(), "");
    r.parameters.reserve(4);
    r.parameters.push_back(rgn.centre().x);
    r.parameters.push_back(rgn.centre().y);
    r.parameters.push_back(rgn.x_dim());
    r.parameters.push_back(rgn.y_dim());
    return r;
}

template DS9Region DS9Region::box(const RectangularRegion<coordinate_types::pixels>&);
template DS9Region DS9Region::box(const RectangularRegion<coordinate_types::CCD>&);
template DS9Region DS9Region::box(const RectangularRegion<coordinate_types::mosaic>&);

vector<DS9Region> DS9Region::load_regions_from_file(const std::string& filename)
{

    ifstream rgnfile(filename.c_str(), std::ios::binary);
    if (!rgnfile.is_open())
        throw std::runtime_error("DS9Region::load_regions_from_file - "
                                 "Unable to open DS9Region file \""+ filename +"\" for reading");

    //skip comments
//    char c;

    while (rgnfile.peek()=='#') {
        rgnfile.ignore(std::numeric_limits<streamsize>::max(), '\n');
    }

    string s;
    rgnfile>>s;
    if (s!="global") {
        throw runtime_error("DS9Region::load_regions_from_file-"
                            "Expected \"global\" keyword, not found in file "
                            + filename + " - got "+s);
    }
    //no parsing of global settings currently implemented, YAGNI
    rgnfile.ignore(std::numeric_limits<streamsize>::max(), '\n');

    string coord_ref_frame_for_vector;
    rgnfile>>coord_ref_frame_for_vector;

    if (!coord_ref_frame_for_vector.empty()) {
        vector<string> co_types = coordinate_types::all_ds9_strings();
        if (find(co_types.begin(), co_types.end(), coord_ref_frame_for_vector)==co_types.end()) {
            throw runtime_error(
                "DS9Region::load_regions_from_file-"
                "Unrecognised co-ord type keyword, found in file " + filename);
        }
    }


    rgnfile.ignore(std::numeric_limits<streamsize>::max(), '\n');

    vector<DS9Region> rgn_vec;
    if (coord_ref_frame_for_vector.empty()) { return rgn_vec; }

    while (rgnfile.peek()!=EOF) {
        if (rgnfile.peek()!='#') {
            DS9Region r;
            r.coord_ref_frame=coord_ref_frame_for_vector;
            rgnfile>>r;
            rgn_vec.push_back(r);
        } else if (rgnfile.peek()=='#') {
            throw runtime_error("DS9Region::load_regions_from_file-"
                                "Unrecognized DS9Region type"); //YAGNI (see wikipedia )
        }

        while (rgnfile.peek()=='\n') {
            rgnfile.get();
        }
    }

    return rgn_vec;
}

void DS9Region::save_regions_to_file(const std::string& filename,
                                     const std::vector<DS9Region>& rgn_vec)
{

    string global_settings=
        "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0";
    string assoc_filename="un_implemented";

    string vec_coord_type;

    if (!rgn_vec.empty()) {
        vec_coord_type = rgn_vec.front().coord_ref_frame;
        for (size_t i=0; i!=rgn_vec.size(); ++i) {
            if (rgn_vec[i].coord_ref_frame != vec_coord_type)
                throw runtime_error("DS9Region::save_regions_to_file - "
                                    "cannot save vector of regions with mixed co-ord types");

        }
    }

    ofstream rgnfile(filename.c_str(), std::ios::binary);
    if (!rgnfile.is_open())
        throw std::runtime_error("DS9Region::save_regions_to_file - "
                                 "Unable to open DS9Region file \""+ filename +"\" for writing");

    rgnfile<< "# Region file format: DS9 version 4.1\n";
    rgnfile <<"# Filename: " +assoc_filename +"\n";
    rgnfile <<global_settings +"\n";
    if (!vec_coord_type.empty()) { rgnfile <<vec_coord_type +"\n"; }

    rgnfile.precision(12);
    for (size_t i=0; i!=rgn_vec.size(); ++i) {
        rgnfile<<rgn_vec[i]<<endl;
    }
    rgnfile<<endl;
    rgnfile.close();
}

template<class ref_frame>
DS9Region::operator Position<ref_frame>() const
{
    if (region_type!="point") {
        throw logic_error("DS9Region::operator Position<T> "
                          "- Cannot convert this DS9Region to a point");
    }
    if (coord_ref_frame!=ref_frame::ds9_string()) {
        throw logic_error("DS9Region::operator Position<T> - "
                          "Cannot convert this ds9 point to Position:"
                          +ref_frame::ds9_string());
    }
    assert(parameters.size()==2);
    return Position<ref_frame>(parameters[0] , parameters[1]);
}

template DS9Region::operator Position<coordinate_types::pixels>() const;
template DS9Region::operator Position<coordinate_types::CCD>() const;
template DS9Region::operator Position<coordinate_types::mosaic>() const;

template<class ref_frame>
DS9Region::operator RectangularRegion<ref_frame>() const
{
    if (region_type!="box") {
        throw logic_error("DS9Region::operator RectangularRegion<T> "
                          "- Cannot convert this DS9Region to a rectangular region");
    }
    if (coord_ref_frame!=ref_frame::ds9_string()) {
        throw logic_error("DS9Region::operator RectangularRegion<T> "
                          "- Cannot convert this ds9 box to coord type:"
                          +ref_frame::ds9_string());
    }
    assert(parameters.size()>=4);
    Position<ref_frame> centre(parameters[0] , parameters[1]);
    TwoVector<ref_frame> diagonal(parameters[2] , parameters[3]);
    diagonal/=2.0;

    return RectangularRegion<ref_frame>(centre-diagonal, centre+diagonal);
}
template DS9Region::operator RectangularRegion<coordinate_types::pixels>() const;
template DS9Region::operator RectangularRegion<coordinate_types::CCD>() const;




//=============================================================================

template<class std_rgn_type>
std::vector< std_rgn_type > load_vector_from_file(const std::string& filename)
{
    vector<DS9Region> ds9_rgns = DS9Region::load_regions_from_file(filename);
    vector< std_rgn_type > posns; posns.reserve(ds9_rgns.size());
    for (size_t i=0; i!=ds9_rgns.size(); ++i) {
        posns.push_back(std_rgn_type(ds9_rgns[i]));
    }
    return posns;
}
template std::vector< PixelPosition > load_vector_from_file(const std::string& filename);
template std::vector< CCD_Position > load_vector_from_file(const std::string& filename);
template std::vector< MosaicPosition > load_vector_from_file(const std::string& filename);

template std::vector< PixelBoxRegion > load_vector_from_file(const std::string& filename);
template std::vector< CCD_BoxRegion > load_vector_from_file(const std::string& filename);

//=============================================================================
} //end of namespace coela::ds9
} //end of namespace coela
