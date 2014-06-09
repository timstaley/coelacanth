cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../CMake)

set(BIN_FOLDER bin)
set(LIB_FOLDER lib)
set(INCLUDE_FOLDER include/coela)
set(TEST_FOLDER ${CMAKE_BINARY_DIR}/unittests)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
enable_testing()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
include_directories( ../../ )
add_subdirectory(../../coela_utility/src bin_utility)
add_subdirectory(../../coela_random/src bin_random)
add_subdirectory(../../coela_core/src bin_core)
add_subdirectory(../src bin_analysis)

#set(CMAKE_VERBOSE_MAKEFILE TRUE) #For debug purposes

