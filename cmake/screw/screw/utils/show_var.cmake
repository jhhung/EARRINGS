if(NOT DEFINED screw_show_var_included)
set(screw_show_var_included)
    
macro(screw_show_var var)
    message(STATUS "${var}: ${${var}}")
endmacro()
macro(screw_show_var_debug var)
    if(SCREW_DEBUG)
        message(STATUS "${var}: ${${var}}")
    endif()
endmacro()

else()
endif()