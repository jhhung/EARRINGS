if(NOT DEFINED screw_set_debugger_included)
set(screw_set_debugger_included)

macro(screw_set_debugger)
    if(NOT CMAKE_CROSSCOMPILING)
        if(GNU OR MINGW)
            if(NOT SCREW_DEBUGGER)
                get_filename_component(
                    __compiler_bin_dir
                    ${CMAKE_CXX_COMPILER}
                    DIRECTORY
                )
                if(MINGW)
                    set(__name gdb.exe)
                else()
                    set(__name gdb)
                endif()
                set(SCREW_DEBUGGER ${__compiler_bin_dir}/${__name})
                unset(__compiler_bin_dir)
                unset(__name)
            endif()
        endif()
    endif()
endmacro()

else()
endif()