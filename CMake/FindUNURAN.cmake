# http://statmath.wu.ac.at/unuran/

if(NOT UNURAN_FOUND)

#  message(STATUS "Home dir is $ENV{HOME}")

  find_path(UNURAN_INCLUDE_DIR unuran.h
    HINTS $ENV{HOME}/local  $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES include)

  find_library(UNURAN_LIBRARY unuran    
    HINTS $ENV{HOME}/local $ENV{LIBRARY_PATH}
    PATH_SUFFIXES lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(UNURAN DEFAULT_MSG
    UNURAN_LIBRARY UNURAN_INCLUDE_DIR)

  set(UNURAN_INCLUDE_DIRS ${UNURAN_INCLUDE_DIR})
  set(UNURAN_LIBRARIES ${UNURAN_LIBRARY})

endif(NOT UNURAN_FOUND)
