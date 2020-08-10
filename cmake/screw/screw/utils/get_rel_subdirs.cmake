if(NOT DEFINED screw_get_rel_subdirs_included)
set(screw_get_rel_subdirs_included)


macro(screw_get_rel_subdirs result curdir)
    file(GLOB children RELATIVE ${curdir} ${curdir}/*)
    set(dirlist "")
    foreach(child ${children})
        if(IS_DIRECTORY ${curdir}/${child})
            list(APPEND dirlist ${child})
        endif()
    endforeach()
    set(${result} ${dirlist})
    unset(children)
    unset(dirlist)
endmacro()


else()
endif()