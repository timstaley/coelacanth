project(coela_utility) #set project name to enable sub-folders for unit test outputs.

find_package(Boost REQUIRED)
find_package(UnitTestpp REQUIRED)


#--------------------------------------------------------------------
#library
file(GLOB sources implementation/*.cc)
file(GLOB headers *.h)
add_library(coela_utility ${sources})
    install(TARGETS coela_utility
	    ARCHIVE DESTINATION ${LIB_FOLDER}
	    LIBRARY DESTINATION ${LIB_FOLDER}
    )

    install(FILES ${headers}
        DESTINATION ${INCLUDE_FOLDER}
    )

#--------------------------------------------------------------------
#Unit tests
file(GLOB test_sources unit_tests/*.cc)
add_executable(test_coela_utility ${test_sources} run_coela_utility_unit_tests.cc)
    target_link_libraries(test_coela_utility 
                        coela_utility 
                        ${UnitTestpp_LIBRARIES}
                        boost_system boost_filesystem)
    install(TARGETS test_coela_utility DESTINATION ${TEST_FOLDER})

add_dependencies(check test_coela_utility)
file(MAKE_DIRECTORY ${TEST_FOLDER})

add_test(NAME runtest_coela_utility
        WORKING_DIRECTORY ${TEST_FOLDER}
        COMMAND  test_coela_utility)
