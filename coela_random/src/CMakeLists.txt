project(coela_random)

find_package(Boost REQUIRED)
find_package(RNGStreams REQUIRED)
find_package(UNURAN REQUIRED)
find_package(UnitTestpp REQUIRED)

#--------------------------------------------------------------------
#library
file(GLOB sources implementation/*.cc)
file(GLOB headers *.h)
add_library(coela_random ${sources})
    target_link_libraries(coela_random
	     ${RNGStreams_LIBRARIES} ${UNURAN_LIBRARIES} 
	     boost_system boost_filesystem boost_serialization
	    )

    install(TARGETS coela_random
	    ARCHIVE DESTINATION ${LIB_FOLDER}
	    LIBRARY DESTINATION ${LIB_FOLDER}
    )

    install(FILES ${headers}
        DESTINATION ${INCLUDE_FOLDER}
    )
#--------------------------------------------------------------------
#Unit tests
file(GLOB test_sources unit_tests/*.cc)
add_executable(test_coela_random ${test_sources} run_coela_random_unit_tests.cc)
    target_link_libraries(test_coela_random 
                            coela_random ${UnitTestpp_LIBRARIES} )
    install(TARGETS test_coela_random DESTINATION ${TEST_FOLDER})

add_dependencies(check test_coela_random)

message(STATUS "TEST_DIR=${TEST_FOLDER}")
file(MAKE_DIRECTORY ${TEST_FOLDER})

add_test(NAME runtest_coela_random
        WORKING_DIRECTORY ${TEST_FOLDER}
        COMMAND  test_coela_random)


