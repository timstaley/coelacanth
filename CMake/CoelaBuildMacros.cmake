
#------------------------------------------------------------------------------
# coela_print_compiler_flags()
# -Print a concatenated list of all compiler flags 
#  (can check this for yourself, try `make VERBOSE=1`, but this makes it more 
#  obvious).
macro(coela_print_compiler_flags)
  if (CMAKE_BUILD_TYPE STREQUAL "debug")
  message(STATUS "CXX_FLAGS=${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  elseif (CMAKE_BUILD_TYPE STREQUAL "release")
  message(STATUS "CXX_FLAGS=${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  endif()
endmacro(coela_print_compiler_flags)

#------------------------------------------------------------------------------
# coela_determine_build_type()
# -sets CMAKE_BUILD_TYPE depending upon the current build directory location,
#  relative to the top-level CMakeLists file.
#  i.e. are we in TOPDIR/build/debug or TOPDIR/build/release?
#  NB Only runs if CMAKE_BUILD_TYPE is not already set.
#  If the build-dir path does not match the presets, default to debug-build.
macro(coela_determine_build_type)
  if (NOT CMAKE_BUILD_TYPE )
    if ( "${CMAKE_CURRENT_SOURCE_DIR}/build/debug" STREQUAL "${PROJECT_BINARY_DIR}" )
      set( CMAKE_BUILD_TYPE debug)
      message(STATUS "Building in debug mode (debug build subfolder used)")
    elseif ( "${CMAKE_CURRENT_SOURCE_DIR}/build/release" STREQUAL "${PROJECT_BINARY_DIR}" )
      set( CMAKE_BUILD_TYPE release)
      message(STATUS "Building in release mode (release build subfolder used)")
    else()
      set( CMAKE_BUILD_TYPE debug)
      message(STATUS "Building in debug mode (default)")
    endif()

  else()
    message(STATUS "Building in ${CMAKE_BUILD_TYPE} mode (user specified)")
  endif(NOT CMAKE_BUILD_TYPE )
endmacro(coela_determine_build_type)


#------------------------------------------------------------------------------
# coela_add_library(exec_name ${lib_sources}):
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
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# coela_add_executable(exec_name ${exec_sources}):
# -Adds executable target
# -Installs target to ${BIN_FOLDER}, 
# -Adds '_debug' suffix if CMAKE_BUILD_TYPE is debug.

macro(coela_add_executable _name)
  add_executable(${_name} ${ARGN})
  install(TARGETS ${_name} DESTINATION ${BIN_FOLDER})
  set_target_properties(${_name} PROPERTIES OUTPUT_NAME_DEBUG ${_name}_debug)
endmacro(coela_add_executable)


#------------------------------------------------------------------------------
# General unit-testing related macros
# (Makes use of UnitTest++ framework)

enable_testing()

# Workaround the cmake bug that 'make test' won't by default *build the test* first:
# http://www.cmake.org/Wiki/CMakeEmulateMakeCheck
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} )

#Place unittests in a separate binary folder to general usage programs:
set(TEST_BIN_FOLDER ${CMAKE_BINARY_DIR}/unittests )

# coela_add_test_exec(exec_name ${exec_sources})
# -Add test executable
# -Change binary dir from default
# -Add _debug suffix if building in debug mode.
# -Add to 'make check' dependency list
# -Run unit-test in project-specific working directory,
#  to separate any output files from other projects.
# -Pass the path to TEST_RESOURCE_FOLDER (if defined) as arg to unittest runner,
#  in case it's useful.
macro(coela_add_test_exec _name)
  add_executable(${_name} ${ARGN})
  install(TARGETS ${_name} DESTINATION ${TEST_BIN_FOLDER})
  set_target_properties(${_name} PROPERTIES OUTPUT_NAME_DEBUG ${_name}_debug)
  add_dependencies(check ${_name})    

  file(MAKE_DIRECTORY ${TEST_BIN_FOLDER}/${PROJECT_NAME} )
  add_test(NAME ${_name}
          WORKING_DIRECTORY ${TEST_BIN_FOLDER}/${PROJECT_NAME}
          COMMAND $<TARGET_FILE:${_name}> ${TEST_RESOURCE_FOLDER}
          )
endmacro(coela_add_test_exec)

#------------------------------------------------------------------------------
