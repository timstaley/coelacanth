/*
 * File:   ds9_region.h
 * Author: ts337
 *
 * Created on 24 November 2008, 15:23
 */

#ifndef COELA_DS9_REGIONS_H
#define COELA_DS9_REGIONS_H
#include "cartesian_coords.h"
#include <vector>
#include <string>

namespace coela {
namespace ds9 {

//========================================================================================
/**

A ds9 representation of any co-ordinate point or region;
generally stored as a line of text in a .reg file.
All constructors are via static functions ("Named constructor idiom")
since this helps to distinguish the different cases.
*/
class DS9Region {
public:
    DS9Region() {}

    ///Construct 'point' marker from specified Position.
    template<class T>
    static DS9Region point(const Position<T>&,
                           const std::string marker_type_string="x");

    ///Construct box from RectangularRegion
    template<class T>
    static DS9Region box(const RectangularRegion<T>&);

    static std::vector<DS9Region> load_regions_from_file(
        const std::string& filename);

    static void save_regions_to_file(const std::string& filename,
                                     const std::vector<DS9Region>&);

    ///Cast to Position
    template<class ref_frame>
    operator Position<ref_frame>() const;

    ///Cast to RectangularRegion
    template<class ref_frame>
    operator RectangularRegion<ref_frame>() const;


private:
    friend std::ostream& operator<<(std::ostream& , const DS9Region&);
    friend std::istream& operator>>(std::istream& , DS9Region&);
    //------------------------------------------------------------------------------------------
    //Data members:
    std::vector<double> parameters;
    std::string region_type, coord_ref_frame, sub_type_comment;
    //------------------------------------------------------------------------------------------

    //handy string constructor
    DS9Region(const std::string& type, const std::string& ref_frame,
              const std::string& sub_type);
};


/**
 Convenience function: Loads DS9 regions and converts to internal std_rgn_type.
 e.g. load a DS9 region file listing points, convert them all to
 'Position' type.
 */

template<class std_rgn_type>
std::vector< std_rgn_type > load_vector_from_file(const std::string& filename);

} //end of namespace coela::ds9
} //end of namespace coela

#endif  /* _DS9_REGIONS_H */

