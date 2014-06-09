/*
 * File:   serialization_body_includes.h
 * Author: ts337
 *
 * Created on 20 April 2011, 15:28
 */

#ifndef COELA_MY_BOOST_SERIALIZATION_BODY_INCLUDES_H
#define COELA_MY_BOOST_SERIALIZATION_BODY_INCLUDES_H

#include "boost_serialization_header_includes.h"

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include<fstream>
#include<sstream>


namespace coela {
namespace boost_serialization {

template <class T>
void save_serializable_object_to_ostream(std::ostream& ss, const T& object)
{
    boost::archive::xml_oarchive oa(ss);
    oa << BOOST_SERIALIZATION_NVP(object);
}

template <class T>
T load_serializable_object_from_istream(std::istream& ss)
{
    T object;
    boost::archive::xml_iarchive ia(ss);
    ia >> BOOST_SERIALIZATION_NVP(object);
    return object;
}


template <class T>
void save_serializable_object_to_file(const std::string& filename, const T& object)
{
    std::ofstream ofs(filename.c_str());
    assert(ofs.good());
    boost::archive::xml_oarchive oa(ofs);
    oa << BOOST_SERIALIZATION_NVP(object);
}

template <class T>
T load_serializable_object_from_file(const std::string& filename)
{
    T object;
    std::ifstream ifs(filename.c_str());
    assert(ifs.good());
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(object);
    return object;
}


template <class T>
void save_serializable_vector_to_file(const std::string& filename,
                                      const std::vector<T>& v)
{
    std::ofstream ofs(filename.c_str());
    assert(ofs.good());
    boost::archive::xml_oarchive oa(ofs);
    oa << BOOST_SERIALIZATION_NVP(v);
}

template <class T>
std::vector<T> load_serializable_vector_from_file(const std::string& filename)
{
    std::vector<T> obs_vec;
    std::ifstream ifs(filename.c_str());
    assert(ifs.good());
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(obs_vec);
    return obs_vec;
}

}//end namespace coela::boost_serialization
}//end namespace coela


#endif  /* MY_BOOST_SERIALIZATION_BODY_INCLUDES_H */

