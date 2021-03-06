if(NOT DEFINED screw_get_qt_win_toolset_included)
set(screw_get_qt_win_toolset_included)

include(${SCREW_DIR}/utils/get_bits.cmake)
include(${SCREW_DIR}/utils/show_var.cmake)

macro( screw_get_qt_win_toolset result)
    screw_get_bits(bits)
    if(MSVC)
        if(MSVC_VERSION EQUAL 1500)
            set(ide_ver 2008)
        elseif(MSVC_VERSION EQUAL 1600)
            set(ide_ver 2010)
        elseif(MSVC_VERSION EQUAL 1700)
            set(ide_ver 2012)
        elseif(MSVC_VERSION EQUAL 1800)
            set(ide_ver 2013)
        elseif(MSVC_VERSION EQUAL 1900)
            set(ide_ver 2015)
        elseif((MSVC_VERSION EQUAL 1910) OR (MSVC_VERSION GREATER 1910))
            set(ide_ver 2017)
        endif()
    endif()
    set(${result} "msvc${ide_ver}_${bits}")
    # message(STATUS "win_toolset: ${${result}}")
endmacro()

else()
endif()