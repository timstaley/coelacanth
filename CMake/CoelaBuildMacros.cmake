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
endmacro(add_component_library)
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



