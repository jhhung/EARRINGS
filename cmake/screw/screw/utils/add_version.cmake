if(NOT DEFINED screw_version_included)
set(screw_version_included)

macro(screw_add_version major minor micro)
    set(PROJECT_MAJOR_VERSION ${major})
    set(PROJECT_MINOR_VERSION ${minor})
    set(PROJECT_MICRO_VERSION ${micro})
    set(PROJECT_VERSION 
        ${major}.${minor}.${micro}
    )
    set(${PROJECT_NAME}_MAJOR_VERSION ${major})
    set(${PROJECT_NAME}_MINOR_VERSION ${minor})
    set(${PROJECT_NAME}_MICRO_VERSION ${micro})
    set(${PROJECT_NAME}_VERSION 
        ${major}.${minor}.${micro}
    )
    set(optional_args ${ARGN})
    screw_show_var(optional_args)
    list(GET optional_args 0 version_header_path)
    screw_show_var(version_header_path)
    if(version_header_path)
        configure_file(
            ${SCREW_DIR}/utils/add_version/version.h.in 
            ${version_header_path} 
            @ONLY
        )
    endif()
    unset(optional_args)
    unset(version_header_path)
endmacro(screw_add_version)


else()
endif()