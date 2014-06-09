
#include <UnitTest++/UnitTest++.h>
#include "../simple_serialization.h"
#include "../string_utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <exception>
#include <sstream>
#include <boost/filesystem.hpp>

using namespace coela;
using namespace std;


extern string lucky_lib_test_resources_dir;

SUITE(simple_serialization)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(convert_relevant_istream_section_and_whole_istream_conversion) {
        using namespace simple_serialization;

        vector<string> content;
        content.push_back("load of stuff");
        content.push_back("morestuff");
        content.push_back("cheese");
        const string start_mark= "START", end_mark="END";

        stringstream s1,s2;
        s1<<" rubbish! rubbish!"<<endl;
        s1 << start_mark << endl;
        for (size_t i=0; i!=content.size(); ++i) {
            s1<<content[i]<<endl;
            s2<<content[i]<<endl;
        }
        s1 << end_mark << endl;
        s1<<" rubbish! rubbish!"<<endl;

        vector<string> svec1 = convert_relevant_istream_section(s1, start_mark, end_mark);
        vector<string> svec2 = convert_istream_to_string_vec(s2);

        CHECK_EQUAL(content.size(), svec1.size());
        CHECK_EQUAL(content.size(), svec2.size());

        for (size_t i=0; i!=content.size(); ++i) {
            CHECK_EQUAL(content[i], svec1[i]);
            CHECK_EQUAL(content[i], svec2[i]);
        }

    }

    TEST(Key_value_pair_line_serialization) {
        using namespace simple_serialization;
        string key1="cheese types? ", value1= "cheddar, edam-stilton-hybrid, brie";
        string key2="some number...";
        double value2 = 42.42;
        stringstream s;


        stringstream ss;
        ss <<"blah blah blah stuff and things"<<endl;
        output_key_value_pair_to_stream(ss, key1, value1);
        output_key_value_pair_to_stream(ss, key2, value2);
        ss<<"more stuff\n";

//        cout<<"SS:\n"<<ss.str()<<endl;

        vector<string> svec = convert_istream_to_string_vec(ss);
//
//        for (size_t i=0; i!=svec.size(); ++i){
//            cout<<i<<":" <<svec[i]<<endl;
//        }

        string value1_output;
        double value2_output;

        get_keyed_value_from_string_vec(svec, key1, value1_output);
        CHECK_EQUAL(value1.size(),value1_output.size());
        CHECK_EQUAL(value1, value1_output);

//        cout<<"Value1 returned:"<<value1_output<<endl;
        get_keyed_value_from_string_vec(svec, key2, value2_output);
        CHECK_EQUAL(value2, value2_output);
    }





}
