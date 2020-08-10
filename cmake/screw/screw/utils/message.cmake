if(NOT DEFINED screw_message_included)
    
macro(screw_message str)
    message(STATUS "screw [MSG}: ${str}")
endmacro()
macro(screw_message_debug str)
    if(SCREW_DEBUG)
        message(STATUS "screw [DBG]: ${str}")
    endif()
endmacro()

else()
set(screw_message_included)
endif()