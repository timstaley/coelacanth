project(analysis) 

find_package(UnitTestpp REQUIRED)
find_package(Minuit2 REQUIRED)
find_package(Boost REQUIRED)
find_package(RNGStreams REQUIRED)
find_package(UNURAN REQUIRED)

#--------------------------------------------------------------------
#library
file(GLOB sources implementation/*.cc)
file(GLOB headers *.h)

add_library(coela_analysis ${sources})
    target_link_libraries(coela_analysis	 
                        coela_core
                        coela_random
                        ${Minuit2_LIBRARIES}
                        )
    install(TARGETS coela_analysis
	    ARCHIVE DESTINATION ${LIB_FOLDER}
	    LIBRARY DESTINATION ${LIB_FOLDER}
    )
    install(FILES ${headers}
        DESTINATION ${INCLUDE_FOLDER}
    )

#--------------------------------------------------------------------
#Unit tests

file(GLOB test_sources unit_tests/*.cc)
add_executable(test_coela_analysis ${test_sources} run_coela_analysis_unit_tests.cc)
    target_link_libraries(test_coela_analysis 
                        coela_analysis 
                        ${UnitTestpp_LIBRARIES})
    install(TARGETS test_coela_analysis DESTINATION ${TEST_FOLDER})

add_dependencies(check test_coela_analysis)
file(MAKE_DIRECTORY ${TEST_FOLDER})
add_test(NAME runtest_coela_analysis
        WORKING_DIRECTORY ${TEST_FOLDER}
        COMMAND  test_coela_analysis)
