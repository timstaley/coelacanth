# http://statmath.wu.ac.at/software/RngStreams/

if(NOT RNGStreams_FOUND)

#  message(STATUS "Home dir is $ENV{HOME}")

  find_path(RNGStreams_INCLUDE_DIR RngStream.h
    HINTS $ENV{HOME}/local  $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES include)

  find_library(RNGStreams_LIBRARY rngstreams    
    HINTS $ENV{HOME}/local $ENV{LIBRARY_PATH}
    PATH_SUFFIXES lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(RNGStreams DEFAULT_MSG
    RNGStreams_LIBRARY RNGStreams_INCLUDE_DIR)

  set(RNGStreams_INCLUDE_DIRS ${RNGStreams_INCLUDE_DIR})
  set(RNGStreams_LIBRARIES ${RNGStreams_LIBRARY})

endif(NOT RNGStreams_FOUND)
