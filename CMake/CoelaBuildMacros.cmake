#-----------------------------------------------------------------------------------
# coela_add_library():
# -Adds library target
# -Installs shared/static build to ${LIB_FOLDER}, 
# -Adds '_debug' suffix if CMAKE_BUILD_TYPE is debug.

macro(coela_add_library _name)
  add_library(${_name} ${ARGN})
  install(TARGETS ${_name}
    ARCHIVE DESTINATION ${LIB_FOLDER}
    LIBRARY DESTINATION ${LIB_FOLDER}
  )
  
#  set_target_properties(${_name} PROPERTIES PREFIX ${LIB_PREFIX})
  set_target_properties(${_name} PROPERTIES OUTPUT_NAME_DEBUG ${_name}_debug)
endmacro(coela_add_library)
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# coela_add_executable():
# -Adds executable target
# -Installs target to ${BIN_FOLDER}, 
# -Adds '_debug' suffix if CMAKE_BUILD_TYPE is debug.

macro(coela_add_executable _name)
  add_executable(${_name} ${ARGN})
  install(TARGETS ${_name} DESTINATION ${BIN_FOLDER})
  set_target_properties(${_name} PROPERTIES OUTPUT_NAME_DEBUG ${_name}_debug)
endmacro(coela_add_executable)


#-----------------------------------------------------------------------------------
# General unit-testing related macros
# (Makes use of UnitTest++ framework)

enable_testing()

# Workaround the cmake bug that 'make test' won't by default *build the test* first:
# http://www.cmake.org/Wiki/CMakeEmulateMakeCheck
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} )

#Place unittests in a separate binary folder to general usage programs:
set(TEST_BIN_FOLDER ${CMAKE_BINARY_DIR}/unittests )

# coela_add_test_exec()
# -Add test executable
# -Change binary dir from default
# -Add to 'make check' dependency list
# -Run unit-test in project-specific working directory,
#  to separate any output files from other projects.
# -Pass the path to TEST_RESOURCE_FOLDER (if defined) as arg to unittest runner,
#  in case it's useful.
macro(coela_add_test_exec _name)
  add_executable(${_name} ${ARGN})
  install(TARGETS ${_name} DESTINATION ${TEST_BIN_FOLDER})
  add_dependencies(check ${_name})    

  file(MAKE_DIRECTORY ${TEST_BIN_FOLDER}/${PROJECT_NAME} )
  add_test(NAME ${_name}
          WORKING_DIRECTORY ${TEST_BIN_FOLDER}/${PROJECT_NAME}
          COMMAND $<TARGET_FILE:${_name}> ${TEST_RESOURCE_FOLDER}
          )
endmacro(coela_add_test_exec)

#-----------------------------------------------------------------------------------
