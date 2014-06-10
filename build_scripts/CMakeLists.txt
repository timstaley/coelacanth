cmake_minimum_required(VERSION 2.8)

set(COELA_TOP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
list(APPEND CMAKE_MODULE_PATH ${COELA_TOP_DIR}/CMake)
include(${COELA_TOP_DIR}/CMake/CoelaBuildMacros.cmake)

set(BIN_FOLDER bin)
set(LIB_FOLDER lib)
set(INCLUDE_FOLDER include/coela)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
include_directories( ${COELA_TOP_DIR} )
include_directories( ${COELA_TOP_DIR}/getoptpp )

add_subdirectory(${COELA_TOP_DIR}/coela_utility/src ${BIN_FOLDER}/coela_utility)
add_subdirectory(${COELA_TOP_DIR}/coela_random/src ${BIN_FOLDER}/coela_random)
add_subdirectory(${COELA_TOP_DIR}/coela_core/src ${BIN_FOLDER}/coela_core)
add_subdirectory(${COELA_TOP_DIR}/coela_analysis/src ${BIN_FOLDER}/coela_analysis)
add_subdirectory(${COELA_TOP_DIR}/coela_luckypipe/src ${BIN_FOLDER}/coela_luckypipe)
add_subdirectory(${COELA_TOP_DIR}/getoptpp ${BIN_FOLDER}/getoptpp)
add_subdirectory(${COELA_TOP_DIR}/luckycam_pipeline ${BIN_FOLDER}/luckycam_pipeline)

#set(CMAKE_VERBOSE_MAKEFILE TRUE) #For debug purposes

