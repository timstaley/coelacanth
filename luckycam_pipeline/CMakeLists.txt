project(luckycam_pipeline) 

set(pipeline_deps 
    getoptpp
    coela_luckypipe)
#(And all that those imply.)

set(LUCKY_BIN_PREFIX lucky_)
macro(coela_add_luckypipe_exec _name)
  coela_add_executable(${_name} ${ARGN})
  set_target_properties(${_name} PROPERTIES PREFIX ${LUCKY_BIN_PREFIX})
  target_link_libraries(${_name} ${pipeline_deps})
endmacro(coela_add_luckypipe_exec)

coela_add_luckypipe_exec(dark_frame_bias_calibrator 
                        dark_frame_bias_calibrator.cc)
coela_add_luckypipe_exec(snapshot 
                        snapshot.cc)
coela_add_luckypipe_exec(clean_and_register 
                        clean_and_register.cc)
coela_add_luckypipe_exec(drizzle 
                        drizzle.cc
                        )


