project(coela_utility) #set project name to enable sub-folders for unit test outputs.

find_package(Boost REQUIRED)
find_package(UnitTestpp REQUIRED)


#--------------------------------------------------------------------
#library
file(GLOB sources implementation/*.cc)
file(GLOB headers *.h)
coela_add_library(coela_utility ${sources})

install(FILES ${headers}
    DESTINATION ${INCLUDE_FOLDER}
)

#--------------------------------------------------------------------
#Unit tests
file(GLOB test_sources unit_tests/*.cc)

coela_add_test_exec(test_coela_utility 
        ${test_sources}
        run_coela_utility_unit_tests.cc
        )
target_link_libraries(test_coela_utility 
                    coela_utility 
                    ${UnitTestpp_LIBRARIES}
                    boost_system boost_filesystem
                    )

