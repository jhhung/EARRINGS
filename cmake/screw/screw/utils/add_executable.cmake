if(NOT DEFINED screw_add_executable_included)
set(screw_add_executable_included)

macro(screw_add_executable target_name )
    add_executable(${target_name} ${ARGN})
    list(APPEND __screw_exe_targets ${target_name})
endmacro()

else()
endif()