# https://www.threadingbuildingblocks.org/

if(NOT TBB_FOUND)

  if(TBB_SRC_DIR)
    file(GLOB _TBB_SRC_DEBUG_LIB ${TBB_SRC_DIR}/build/*/libtbb_debug.so)
    get_filename_component(_TBB_SRC_DEBUG_LIB_DIR ${_TBB_SRC_DEBUG_LIB} PATH)

    file(GLOB _TBB_SRC_REL_LIB ${TBB_SRC_DIR}/build/*/libtbb.so)
    get_filename_component(_TBB_SRC_REL_LIB_DIR ${_TBB_SRC_REL_LIB} PATH)
  endif(TBB_SRC_DIR)

  find_path(TBB_INCLUDE_DIR tbb/tbb.h
    HINTS $ENV{HOME}/local ${TBB_SRC_DIR}  $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES include tbb)

  find_library(TBB_DEBUG_LIBRARY tbb_debug    
    HINTS $ENV{HOME}/local ${_TBB_SRC_DEBUG_LIB_DIR} $ENV{LIBRARY_PATH}
    PATH_SUFFIXES lib)

  find_library(TBB_LIBRARY tbb
    HINTS $ENV{HOME}/local ${_TBB_SRC_REL_LIB_DIR} $ENV{LIBRARY_PATH}
    PATH_SUFFIXES lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(TBB DEFAULT_MSG
    TBB_LIBRARY TBB_INCLUDE_DIR)

  find_package_handle_standard_args(TBB_DEBUG DEFAULT_MSG
    TBB_DEBUG_LIBRARY TBB_INCLUDE_DIR)

  set(TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIR})
  set(TBB_LIBRARIES ${TBB_LIBRARY})    
  if (UNIX)
    if (NOT APPLE) #i.e. Linux:
      set(TBB_LIBRARIES ${TBB_LIBRARY} rt)    
    endif(NOT APPLE)  
  endif(UNIX)

  set(TBB_DEBUG_LIBRARIES ${TBB_DEBUG_LIBRARY})

endif(NOT TBB_FOUND)

