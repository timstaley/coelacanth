/*
 * File:   histogram_Container_tests.cc
 * Author: ts337
 *
 * Created on 07 February 2011, 12:41
 */
#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include "../boost_serialization_body_includes.h"
#include <boost/filesystem.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <iostream>
#include <sstream>
using namespace coela;
using namespace std;

SUITE(Serialization)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(basic_variable_serialization) {
        double x=12.45;
        stringstream ss;
        boost_serialization::save_serializable_object_to_ostream(ss,x);
        double y =   boost_serialization::load_serializable_object_from_istream<double>(ss);
        CHECK_EQUAL(x,y);
    }


    class DummyClass {
    public:
        double x,y;
    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int /* file_version */) {
            ar  & BOOST_SERIALIZATION_NVP(x)
            & BOOST_SERIALIZATION_NVP(y)
            ;
        }

    };
//
    TEST(basic_class_serialization) {
        DummyClass d;
        d.x=1.2;
        d.y=1.3;

        DummyClass e;

        stringstream ss;

//        boost::archive::text_oarchive oa(ss);
//        oa <<  BOOST_SERIALIZATION_NVP(d);
//        boost::archive::text_iarchive ia(ss);
//        ia >> BOOST_SERIALIZATION_NVP(e);

        boost_serialization::save_serializable_object_to_ostream(ss,d);
        e = boost_serialization::load_serializable_object_from_istream<DummyClass>(ss);

        CHECK_EQUAL(d.x, e.x);
        CHECK_EQUAL(d.y, e.y);



    }
}
