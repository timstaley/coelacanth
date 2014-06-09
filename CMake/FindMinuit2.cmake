# http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/

if(NOT Minuit2_FOUND)

  find_path(Minuit2_INCLUDE_DIR Minuit2/Minuit2Minimizer.h
    HINTS $ENV{HOME}/local $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES include)

  find_library(Minuit2_LIBRARY Minuit2    
    HINTS $ENV{HOME}/local $ENV{LIBRARY_PATH}
    PATH_SUFFIXES lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Minuit2 DEFAULT_MSG
    Minuit2_LIBRARY Minuit2_INCLUDE_DIR)

  set(Minuit2_INCLUDE_DIRS ${Minuit2_INCLUDE_DIR})
  set(Minuit2_LIBRARIES ${Minuit2_LIBRARY})

endif(NOT Minuit2_FOUND)
