#include <UnitTest++/UnitTest++.h>
#include <sstream>
#include "../FitsHeader.h"
#include "coela_utility/src/string_utils.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <boost/filesystem.hpp>

using namespace coela;
using namespace std;

SUITE(Fits_Header_Table)
{

    string test_suite_output_dir;
    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"FitsHeader\" unit tests running ***" <<endl;
        test_suite_output_dir = string(UnitTest::CurrentTest::Details()->suiteName) + "_tests/";
        boost::filesystem::create_directories(test_suite_output_dir);

    }

    struct basic_hdr_table {
        basic_hdr_table() {
            tbl.add_keyword("SIMPLE","T", "mandatory standard prefix");
            tbl.add_keyword("BITPIX", string_utils::itoa(fits_header_conventions::FLOATIMG));
            tbl.add_keyword("A_KEY","SOMEVAL","Example key value pair");
            tbl.add_keyword("B_KEY","SOMEVAL","Another key value pair");
            tbl.add_comment("Short comment string");
            tbl.add_keyword("c_key","Some_Other_Val","Another key value pair");
            tbl.add_comment("Very long comment string, we expect this comment to be "
                            "split onto several lines ................ that ought to do it.");
        }

        FitsHeader tbl;
    };

    TEST_FIXTURE(basic_hdr_table, Basic_Functionality) {
        FitsHeader& tbl_ref=tbl;
        CHECK_EQUAL("Some_Other_Val", tbl.get_key_value("c_key"));
        CHECK_EQUAL(tbl_ref.FITS_imagetype(), fits_header_conventions::FLOATIMG);
    }

    TEST_FIXTURE(basic_hdr_table, Serialization) {
        std::stringstream ss1, ss2;
        ss1<<tbl;

        CHECK_EQUAL(tbl.header_file_length_in_bytes(), ss1.str().size());

        FitsHeader tbl_copy;
        ss1>>tbl_copy;

//      CHECK( tbl == tbl_copy); //precise Operator comparison fails  due to whitespace mismatch in comments

        CHECK_EQUAL(tbl.size(), tbl_copy.size());
        CHECK(tbl_copy.key_exists("A_KEY"));
        CHECK_EQUAL(tbl_copy.get_key_value("A_KEY"), "SOMEVAL");
        CHECK(tbl_copy.key_exists("B_KEY"));
        CHECK(tbl_copy.key_exists("c_key"));
        CHECK_EQUAL(tbl_copy.get_key_value("c_key"), "Some_Other_Val");

        ss1.str("");

        ss1<<tbl;
        ss2<<tbl_copy;


//      cout<<"\n\n\nSS1\n";
//      cout<<ss1.str()<<endl;
//      cout<<"\n\n\nSS2\n";
//      cout<<ss2.str()<<endl;
//      cout<<"\n\n\n\n";

        CHECK_EQUAL(ss1.str(), ss2.str());
    }

    TEST_FIXTURE(basic_hdr_table, File_Read_and_Write) {
        tbl.write_to_file("temp_file_hdr1.fits");

        FitsHeader tbl_copy1("temp_file_hdr1.fits");
        std::stringstream ss1, ss2;
        ss1<<tbl;
        ss2<<tbl_copy1;
        CHECK_EQUAL(ss1.str(), ss2.str());

        tbl_copy1.write_to_file("temp_file_hdr2.fits");
        FitsHeader tbl_copy2("temp_file_hdr2.fits");
        //Equality op should work now that it's been output once already so whitespace is standardized
        CHECK(tbl_copy1==tbl_copy2);
    }

    TEST_FIXTURE(basic_hdr_table, Load_from_buffer) {
        tbl.write_to_file("temp_file_hdr1.fits");
        FitsHeader tbl_copy1("temp_file_hdr1.fits");
        tbl_copy1.write_to_file("temp_file_hdr2.fits");
        FileBuffer b;
        b.buffer_whole_file("temp_file_hdr2.fits");
        FitsHeader tbl_loaded_from_buffer(b, 0);

        CHECK(tbl_copy1 == tbl_loaded_from_buffer);

    }

    TEST(no_file_exception) {
        string no_file_name =test_suite_output_dir+"non_existent.lcz";
//       FitsHeader fht(no_file_name);
//       cerr<<"Eh?\n"<<fht<<endl;

        CHECK_THROW(FitsHeader fht(no_file_name), std::runtime_error);
    }

    TEST(bad_file_exception) {
        string bad_file_name =test_suite_output_dir+"bad_file.lcz";
        ofstream bad_file(bad_file_name.c_str());
        bad_file<<"fnUJIDHSDVHDIHVIDFHVCOISH IFDSIFHS\\ FDISHFIASDNFAIHDFASIOHIOSAHIDHASIHFIHGUIDHIUGHSDUGFUISD";
        bad_file.close();

        CHECK_THROW(FitsHeader fht(bad_file_name), std::runtime_error);
    }

}




