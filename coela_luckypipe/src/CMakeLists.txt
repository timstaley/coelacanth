project(pipeline) 

find_package(UnitTestpp REQUIRED)
find_package(Minuit2 REQUIRED)
find_package(Boost REQUIRED)
find_package(RNGStreams REQUIRED)
find_package(UNURAN REQUIRED)
find_package(TBB REQUIRED)

#target_link_libraries(pipeline	 
#                    general analysis boost_date_time
#                    optimized ${TBB_LIBRARIES}
#                    debug ${TBB_DEBUG_LIBRARIES}
#                    )

#--------------------------------------------------------------------
#library
file(GLOB sources implementation/*.cc)
file(GLOB headers *.h)

coela_add_library(coela_luckypipe ${sources})
target_link_libraries(coela_luckypipe 
                    coela_analysis boost_date_time
                    optimized ${TBB_LIBRARIES}
                    debug ${TBB_DEBUG_LIBRARIES}
                    )

install(FILES ${headers}
    DESTINATION ${INCLUDE_FOLDER}
)

#--------------------------------------------------------------------
#add_python_wrapper(python-pipeline lucky_pipeline 
#                    python_wrappers/lucky_pipeline.cc)
#target_link_libraries(python-pipeline 
#                      boost_python pipeline)

#--------------------------------------------------------------------
#Unit tests

file(GLOB test_sources unit_tests/*.cc)
add_executable(test_coela_luckypipe ${test_sources} run_coela_luckypipe_unit_tests.cc)
    target_link_libraries(test_coela_luckypipe 
                        coela_luckypipe 
                        ${UnitTestpp_LIBRARIES})
    install(TARGETS test_coela_luckypipe DESTINATION ${TEST_FOLDER})

add_dependencies(check test_coela_luckypipe)
file(MAKE_DIRECTORY ${TEST_FOLDER})
add_test(NAME runtest_coela_luckypipe
        WORKING_DIRECTORY ${TEST_FOLDER}
        COMMAND  test_coela_luckypipe)

