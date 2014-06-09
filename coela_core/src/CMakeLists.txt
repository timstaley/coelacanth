project(coela_core) 

find_package(Boost REQUIRED)
find_package(UnitTestpp REQUIRED)

#--------------------------------------------------------------------
#library
file(GLOB sources implementation/*.cc)
file(GLOB headers *.h)

add_library(coela_core ${sources})
    target_link_libraries(coela_core	 
                        coela_utility
                        boost_system boost_filesystem
                        )
    install(TARGETS coela_core
	    ARCHIVE DESTINATION ${LIB_FOLDER}
	    LIBRARY DESTINATION ${LIB_FOLDER}
    )
    install(FILES ${headers}
        DESTINATION ${INCLUDE_FOLDER}
    )

#--------------------------------------------------------------------
#Unit tests
file(GLOB test_sources unit_tests/*.cc)
add_executable(test_coela_core ${test_sources} run_coela_core_unit_tests.cc)
    target_link_libraries(test_coela_core 
                        coela_core 
                        ${UnitTestpp_LIBRARIES})
    install(TARGETS test_coela_core DESTINATION ${TEST_FOLDER})

add_dependencies(check test_coela_core)
file(MAKE_DIRECTORY ${TEST_FOLDER})
add_test(NAME runtest_coela_core
        WORKING_DIRECTORY ${TEST_FOLDER}
        COMMAND  test_coela_core ${PROJECT_SOURCE_DIR}/../test_resources)

