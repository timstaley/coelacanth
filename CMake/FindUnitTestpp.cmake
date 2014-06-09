# http://unittest-cpp.sourceforge.net/
if(NOT UnitTestpp_FOUND)

#  message(STATUS "Home dir is $ENV{HOME}")

  find_path(UnitTestpp_INCLUDE_DIR UnitTest++.h
    HINTS $ENV{HOME}/local $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES include/UnitTest++ UnitTest++)

  find_library(UnitTestpp_LIBRARY UnitTest++    
    HINTS $ENV{HOME}/local $ENV{LIBRARY_PATH}
    PATH_SUFFIXES lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(UnitTestpp DEFAULT_MSG
    UnitTestpp_LIBRARY UnitTestpp_INCLUDE_DIR)

  set(UnitTestpp_INCLUDE_DIRS ${UnitTestpp_INCLUDE_DIR})
  set(UnitTestpp_LIBRARIES ${UnitTestpp_LIBRARY})

endif(NOT UnitTestpp_FOUND)
