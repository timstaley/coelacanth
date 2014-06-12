# Preamble
cmake_minimum_required(VERSION 2.8)
set(COELA_TOP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
list(APPEND CMAKE_MODULE_PATH ${COELA_TOP_DIR}/CMake)

#Load macros - these are all prefixed 'coela_'.
include(${COELA_TOP_DIR}/CMake/CoelaBuildMacros.cmake)

#Guess at debug or release build based on current dir.
coela_determine_build_type()

# Subfolders of install destination:
set(BIN_FOLDER bin)
set(LIB_FOLDER lib)
set(INCLUDE_FOLDER include/coela)

# Find the headers
include_directories( ${COELA_TOP_DIR} )
include_directories( ${COELA_TOP_DIR}/getoptpp )

# Carefully tuned compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# NB -O3 -DNDEBUG are the cmake release-flag defaults, but not necessarily 
# what the user wants. Ditch them.
if (CMAKE_CXX_FLAGS_RELEASE STREQUAL "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -march=native -O3")

#ipa-cp-clone causes template specification problems when combined with O3.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-ipa-cp-clone")  

# Slightly faster math, esp. floor() & ceil()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mfpmath=sse -ffast-math -fno-math-errno ") 

# See all the gory details of your compiler flags 
#(helpful when figuring out details of what e.g. O3 actually does)
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fverbose-asm -save-temps ") 

#Print the result
coela_print_compiler_flags()

# Load the code
add_subdirectory(${COELA_TOP_DIR}/coela_utility/src ${BIN_FOLDER}/coela_utility)
add_subdirectory(${COELA_TOP_DIR}/coela_random/src ${BIN_FOLDER}/coela_random)
add_subdirectory(${COELA_TOP_DIR}/coela_core/src ${BIN_FOLDER}/coela_core)
add_subdirectory(${COELA_TOP_DIR}/coela_analysis/src ${BIN_FOLDER}/coela_analysis)
add_subdirectory(${COELA_TOP_DIR}/coela_luckypipe/src ${BIN_FOLDER}/coela_luckypipe)
add_subdirectory(${COELA_TOP_DIR}/getoptpp ${BIN_FOLDER}/getoptpp)
add_subdirectory(${COELA_TOP_DIR}/luckycam_pipeline ${BIN_FOLDER}/luckycam_pipeline)
add_subdirectory(${COELA_TOP_DIR}/examples ${BIN_FOLDER}/examples)

#set(CMAKE_VERBOSE_MAKEFILE TRUE) #For debug purposes


