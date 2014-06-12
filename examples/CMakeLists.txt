project(example_programs) 

set(pipeline_deps 
    getoptpp
    coela_luckypipe)
#(And all that those imply.)

set(COELA_BIN_PREFIX coela_)
macro(coela_add_example_exec _name)
  coela_add_executable(${_name} ${ARGN})
  set_target_properties(${_name} PROPERTIES PREFIX ${COELA_BIN_PREFIX})
  target_link_libraries(${_name} ${pipeline_deps})
endmacro(coela_add_example_exec)

coela_add_example_exec(drizzle_combine
                       drizzle_combine.cc)
                        
coela_add_example_exec(getopt_usage_example
                       getopt_usage_example.cc)

