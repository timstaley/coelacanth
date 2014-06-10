project(coela_core) 

find_package(Boost REQUIRED)
find_package(UnitTestpp REQUIRED)

#--------------------------------------------------------------------
#library
file(GLOB sources implementation/*.cc)
file(GLOB headers *.h)

coela_add_library(coela_core ${sources})

target_link_libraries(coela_core	 
                    coela_utility
                    boost_system boost_filesystem
                    )

install(FILES ${headers}
  DESTINATION ${INCLUDE_FOLDER}
)

#--------------------------------------------------------------------
#Unit tests
file(GLOB test_sources unit_tests/*.cc)

set(TEST_RESOURCE_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/../test_resources)

coela_add_test_exec(test_coela_core 
    ${test_sources} run_coela_core_unit_tests.cc
)
target_link_libraries(test_coela_core 
                    coela_core 
                    ${UnitTestpp_LIBRARIES}
)

