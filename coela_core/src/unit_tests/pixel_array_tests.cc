
#include <UnitTest++/UnitTest++.h>
#include "../pixel_array2d.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <exception>
#include <boost/filesystem.hpp>
#include <typeinfo>

using namespace coela;
using namespace std;


extern string lucky_lib_test_resources_dir;

SUITE(pixel_array_pixel_access)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    TEST(null_constructor) {
        PixelArray2d<float> null_float_array;
        PixelArray2d<double> null_double_array;
    }



    const double test_array_init_value = 5.0;
    struct basic_float_pixel_array {
        basic_float_pixel_array() {
            test_array= PixelArray2d<float>(10,5, 5.0);
        }
        PixelArray2d<float> test_array;

    };

    TEST_FIXTURE(basic_float_pixel_array, SimpleConstructFromBits_float) {
        //check dims
        CHECK_EQUAL(10, test_array.range().x_dim());
        CHECK_EQUAL(5, test_array.range().y_dim());

        //check pixel values
        CHECK_EQUAL(test_array_init_value, test_array(5,1));


    }

    TEST_FIXTURE(basic_float_pixel_array, sum) {
        CHECK_CLOSE(test_array.range().n_pix()*test_array_init_value, test_array.sum() , 1e-6);
    }

    TEST(double_pixel_array_fractional_init) {
        PixelArray2d<double> pixarray(10,10, 1.0/9.0);
        double val = 1.0/9.0;
        CHECK_EQUAL(val, pixarray(1,1));
        CHECK_CLOSE(val*pixarray.range().n_pix(), pixarray.sum(), 1e-6);
    }

    TEST_FIXTURE(basic_float_pixel_array, pixel_iteration_for_float_array) {
        size_t n_its=0;
        for (PixelIterator pi(test_array.range()); pi!=pi.end; ++pi) {
            CHECK_EQUAL(5.0, test_array(pi));
            test_array(pi)-=3.0;
            n_its++;
        }

        CHECK_EQUAL(test_array.range().n_pix(), n_its);

        for (PixelIterator pi(test_array.range()); pi!=pi.end; ++pi) {
            CHECK_EQUAL(2.0,test_array(pi));
        }

    }

    TEST_FIXTURE(basic_float_pixel_array, copy_constructor_float_array) {
        PixelArray2d<float> array2(test_array);
        CHECK(test_array== array2);
    }

    TEST_FIXTURE(basic_float_pixel_array, float_array_assignment) {
        PixelArray2d<float> array2=test_array;
        CHECK(test_array==array2);
    }


    TEST_FIXTURE(basic_float_pixel_array, subarray_constructor_float_array) {
        float x=0.0;
        for (PixelIterator pi(test_array.range()); pi!=pi.end; ++pi) {
            test_array(pi)=x;
            x+=0.1;
        }
        PixelRange subarray_box(1,1, 5,5);
        PixelArray2d<float> subarray = PixelArray2d<float>::sub_array(test_array, subarray_box);

        CHECK(subarray_box.n_pix()==subarray.range().n_pix());

        for (PixelIterator pi_orig(subarray_box), pi_sub(subarray.range());
                pi_orig!=pi_orig.end;
                ++pi_orig, ++pi_sub) {
            CHECK(test_array(pi_orig)== subarray(pi_sub));
        }

    }


    TEST_FIXTURE(basic_float_pixel_array, factor_multiply) {
        double factor=2.0;
        test_array*=factor;
        PixelArray2d<float> twice_val(test_array.range().x_dim(),
                                      test_array.range().y_dim(),
                                      test_array_init_value*factor);

        CHECK(twice_val == test_array);
    }

    TEST_FIXTURE(basic_float_pixel_array, promotion_to_double_type) {
        PixelArray2d<double> promoted(test_array) ;
        for (PixelIterator it(promoted.range()); it!=it.end; ++it) {
            CHECK_CLOSE((double)test_array(it), promoted(it), 1e-6);
        }
    }

    TEST_FIXTURE(basic_float_pixel_array, conversion_between_types) {
        PixelArray2d<double> promoted(test_array);

        promoted.check_header_type_matches_array_type();
        PixelArray2d<int> i(5,5, 3);
//        PixelArray2d<int> f_to_i(test_array);
        PixelArray2d<float> i_to_f(i);

//        PixelArray2d<double> i_to_d(i);
//        PixelArray2d<short int> i_to_si(i);

        PixelArray2d<float> demoted(promoted);
//

//        for (PixelIterator it(promoted.range()); it!=it.end; ++it){
//            CHECK_CLOSE((double)test_array(it), demoted(it), 1e-6);
//        }
    }

    TEST(subtraction_of_double_From_array_of_int) {

        PixelArray2d<int> init(10,10, 5);

        PixelArray2d<double> debiased = init - 2.5;
        for (PixelIterator i(debiased.range()); i!=i.end; ++i) {
            CHECK_EQUAL(2.5, debiased(i));
        }

    }

    TEST_FIXTURE(basic_float_pixel_array, min_max_value_grabs) {
        PixelArray2d<float>& pix=test_array;

        Pixel<float> max(PixelIndex(pix.range().x_dim()/2, pix.range().y_dim()/2),
                         1.563465*1e5);
        Pixel<float> min(PixelIndex(pix.range().x_dim()/2 - 1, pix.range().y_dim()/2),
                         -1.343234*1e5);
        pix(max)=max.value;
        pix(min)=min.value;

        CHECK_EQUAL(PixelIndex(max), pix.max_PixelIndex());
        CHECK_EQUAL(PixelIndex(min), pix.min_PixelIndex());
        CHECK_EQUAL(min.value, pix.min_val());
        CHECK_EQUAL(max.value, pix.max_val());


    }

    //TO DO - test image arithmetic
}


SUITE(restricted_functions_eg_stream_load_and_write)
{

    string test_suite_output_dir= string(UnitTestSuite::GetSuiteName()) + "_tests/";

    TEST(Run_Standard_Suite_Setup) {
        cout << "*** \""<< UnitTestSuite::GetSuiteName() <<"\" unit tests running ***" <<endl;
        boost::filesystem::create_directories(test_suite_output_dir);
    }

    //Create a fixture to be tested:
    const double test_array_init_value = 5.0;
    struct basic_float_pixel_array {
        PixelArray2d<float> test_array;
        basic_float_pixel_array() {
            test_array= PixelArray2d<float>(10,5, 5.0);
        }
    };


    TEST_FIXTURE(basic_float_pixel_array, data_serialization_for_float_no_compression) {

        //NB fht serialization, PixelArrayHeader serialization, both unit tested elsewhere,

        ArrayCompressionInfo compression_inf = ArrayCompressionInfo::no_compression();
        stringstream ss;
        test_array.write_data_to_stream(ss, compression_inf);

        PixelArray2d<float> a2;

        a2.hdr = PixelArray2dHeader(test_array.bitpix(),
                                    test_array.range().x_dim(), test_array.range().y_dim());

        size_t byte_offset=0;
        a2.seek_and_load_data_from_stream(ss, byte_offset, compression_inf);
        CHECK(test_array==a2);
    }

    TEST_FIXTURE(basic_float_pixel_array, file_access_for_float) {
        string filename(test_suite_output_dir +"test_output_file.fits");

        test_array.write_to_file(filename);
        PixelArray2d<float> a2 = PixelArray2d<float>::load_from_file(filename);

        CHECK(test_array.range() == a2.range());
        CHECK(test_array== a2);

        FitsHeader fht(filename);
        CHECK(fht.get_key_value("BZERO")=="0.0");

    }

    TEST_FIXTURE(basic_float_pixel_array, buffered_file_access_for_float) {
        string filename(test_suite_output_dir +"test_output_file.fits");
        test_array.write_to_file(filename);

        PixelArray2d<float> reg_load=
            PixelArray2d<float>::load_from_file(filename);

        FileBuffer buf;
        buf.buffer_whole_file(filename);
        PixelArray2d<float> buf_load=
            PixelArray2d<float>::load_from_buffer(buf);

        CHECK(reg_load==buf_load);
    }

    TEST(file_access_for_short) {
        PixelArray2d<uint16_t> test_array(10,5,25.0);
        string filename(test_suite_output_dir +"test_output_short_file.fits");

        test_array.write_to_file(filename);

        PixelArray2d<uint16_t> a2= PixelArray2d<uint16_t>::load_from_file(filename);

        CHECK(test_array.range() == a2.range());

        CHECK(test_array==a2);

        CHECK_THROW(PixelArray2d<uint32_t>::load_from_file(filename), std::runtime_error);
        CHECK_THROW(PixelArray2d<double>::load_from_file(filename), std::runtime_error);
    }



    TEST(load_from_bitstream) {
        string bitstream_compressed_img_filename(lucky_lib_test_resources_dir+"bit_file.lcz");
        string float_img_filename(lucky_lib_test_resources_dir+"old_code_float_copy.fits");

        PixelArray2d<uint16_t> lcz_short_img =
            PixelArray2d<uint16_t>::load_from_file(bitstream_compressed_img_filename);
        PixelArray2d<uint32_t> lcz_int_img =
            PixelArray2d<uint32_t>::load_from_file(bitstream_compressed_img_filename);
        PixelArray2d<float> lcz_float_img =
            PixelArray2d<float>::load_from_file(bitstream_compressed_img_filename);
        PixelArray2d<double> lcz_double_img =
            PixelArray2d<double>::load_from_file(bitstream_compressed_img_filename);

        lcz_short_img.write_to_file(test_suite_output_dir +"decompressed_short_img.fits");
        lcz_int_img.write_to_file(test_suite_output_dir +"decompressed_int_img.fits");
        lcz_float_img.write_to_file(test_suite_output_dir +"decompressed_float_img.fits");
        lcz_double_img.write_to_file(test_suite_output_dir +"decompressed_double_img.fits");

        PixelArray2d<float> float_img =
            PixelArray2d<float>::load_from_file(float_img_filename);


        for (PixelIterator i(float_img.range()); i!=i.end; ++i) {
            CHECK_EQUAL(lcz_float_img(i), float_img(i));
        }
    }

    TEST(load_from_bitstream_buffer) {

        string bitstream_compressed_img_filename(lucky_lib_test_resources_dir+"bit_file.lcz");
        FileBuffer buf; buf.buffer_whole_file(bitstream_compressed_img_filename);

        PixelArray2d<uint16_t> from_file=
            PixelArray2d<uint16_t>::load_from_file(bitstream_compressed_img_filename);
        PixelArray2d<uint16_t> from_buf=
            PixelArray2d<uint16_t>::load_from_buffer(buf, 0);

        CHECK(from_file == from_buf);

    }

    TEST(load_from_lcc_bytestream) {
        string lcc_byte_stream_compressed_img_filename(
            lucky_lib_test_resources_dir+"lcc_file.lcc");
        string lcc_byte_stream_decompressed_img_filename(
            lucky_lib_test_resources_dir+"decompressed_lcc_short_file.fits");

        PixelArray2d<int16_t> lcc_short_img=
            PixelArray2d<int16_t>::load_from_file(lcc_byte_stream_compressed_img_filename);
        PixelArray2d<int32_t> lcc_uint_img=
            PixelArray2d<int32_t>::load_from_file(lcc_byte_stream_compressed_img_filename);
        for (PixelIterator i(lcc_short_img.range()); i!=i.end; ++i) {
            CHECK_EQUAL(lcc_short_img(i), lcc_uint_img(i));
        }
//        cout<<"Writing to file:"<<test_suite_output_dir +"decompress_lcc_to_short_img.fits"<<endl;
        lcc_short_img.write_to_file(test_suite_output_dir +"decompress_lcc_to_short_img.fits");
        lcc_uint_img.write_to_file(test_suite_output_dir +"decompress_lcc_to_uint_img.fits");

    }

    TEST(load_from_lcc_buffer) {
        string lcc_byte_stream_compressed_img_filename(
            lucky_lib_test_resources_dir+"lcc_file.lcc");

        FileBuffer buf; buf.buffer_whole_file(lcc_byte_stream_compressed_img_filename);

        PixelArray2d<uint16_t> from_file =
            PixelArray2d<uint16_t>::load_from_file(lcc_byte_stream_compressed_img_filename);
        PixelArray2d<uint16_t> from_buf=
            PixelArray2d<uint16_t>::load_from_buffer(buf, 0);

        CHECK(from_file == from_buf);

    }


}
