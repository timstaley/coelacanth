project(coela_random)

find_package(Boost REQUIRED)
find_package(RNGStreams REQUIRED)
find_package(UNURAN REQUIRED)
find_package(UnitTestpp REQUIRED)

#--------------------------------------------------------------------
#library
file(GLOB sources implementation/*.cc)
file(GLOB headers *.h)
coela_add_library(coela_random ${sources})
target_link_libraries(coela_random
     ${RNGStreams_LIBRARIES} ${UNURAN_LIBRARIES} 
     boost_system boost_filesystem boost_serialization
    )
install(FILES ${headers}
    DESTINATION ${INCLUDE_FOLDER}
)
#--------------------------------------------------------------------
#Unit tests
file(GLOB test_sources unit_tests/*.cc)
coela_add_test_exec(test_coela_random 
    ${test_sources} run_coela_random_unit_tests.cc
    )
target_link_libraries(test_coela_random 
                    coela_random ${UnitTestpp_LIBRARIES} 
                    )
