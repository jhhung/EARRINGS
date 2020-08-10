if(NOT DEFINED screw_add_library_included)
set(screw_add_library_included)

macro(screw_add_library target_name )
    add_library(${target_name} ${ARGN})
    list(APPEND __screw_lib_targets ${target_name})
endmacro()

else()
endif()