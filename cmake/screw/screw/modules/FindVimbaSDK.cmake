# SET(CMAKE_FIND_LIBRARY_PREFIXES "")
# SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(VC_LIB_DIR Win64)
else()
  set(VC_LIB_DIR Win32)
endif()

if(WIN32)
    set(VIMBA_DIR_FEATURE "Program Files/Allied Vision/Vimba_2.1/")
    screw_get_win_disk_paths(dps)
    screw_show_var(dps)
    find_path(
        VIMBA_DIR 
        ${VIMBA_DIR_FEATURE}
        HINTS ${dps}
    )
    screw_show_var(VIMBA_DIR)
    set(VIMBA_DIR "${VIMBA_DIR}/${VIMBA_DIR_FEATURE}")
    unset(VIMBA_DIR_FEATURE)


    find_library( 
        VIMBA_C_LIBRARIES NAMES VimbaC 
        HINTS "${VIMBA_DIR}/VimbaC/Lib/${VC_LIB_DIR}"
    )

    find_library( 
        VIMBA_CPP_LIBRARIES NAMES VimbaCPP
        HINTS "${VIMBA_DIR}/VimbaCPP/Lib/${VC_LIB_DIR}"
    )

    find_path( 
        VIMBA_C_INCLUDE_DIRS 
        NAMES VimbaC/Include/VimbaC.h
        HINTS "${VIMBA_DIR}"
    )
    find_path( 
        VIMBA_CPP_INCLUDE_DIRS 
        NAMES VimbaCPP/Include/VimbaCPP.h 
        HINTS "${VIMBA_DIR}"
    )
    message(STATUS "VIMBA_DIR:              ${VIMBA_DIR}" )
    message(STATUS "VIMBA_LIB_DIR:          ${VIMBA_DIR}/VimbaC/Lib/${VC_LIB_DIR}")
    message(STATUS "VIMBA_C_LIBRARIES:      ${VIMBA_C_LIBRARIES}" )
    message(STATUS "VIMBA_CPP_LIBRARIES:    ${VIMBA_CPP_LIBRARIES}" )


    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args(
        VimbaSDK  DEFAULT_MSG 
        VIMBA_C_LIBRARIES
        VIMBA_CPP_LIBRARIES 
        VIMBA_C_INCLUDE_DIRS
        VIMBA_CPP_INCLUDE_DIRS
    )

    mark_as_advanced(VIMBA_C_INCLUDE_DIRS VIMBA_C_LIBRARIES VIMBA_CPP_INCLUDE_DIRS VIMBA_CPP_LIBRARIES)

    set(VIMBA_INCLUDE_DIRS ${VIMBA_C_INCLUDE_DIRS} ${VIMBA_CPP_INCLUDE_DIRS} "${VIMBA_DIR}/VimbaImageTransform/Include")
    set(VIMBA_LIBRARIES ${VIMBA_C_LIBRARIES} ${VIMBA_CPP_LIBRARIES} "${VIMBA_DIR}/VimbaImageTransform/Lib/${VC_LIB_DIR}/VimbaImageTransform.lib")
    find_path(
        VIMBA_C_RT_DIR NAMES VimbaC.dll
        HINTS "${VIMBA_DIR}/VimbaC/Bin/${VC_LIB_DIR}"
    )
    find_path(
        VIMBA_CPP_RT_DIR NAMES VimbaCPP.dll
        HINTS "${VIMBA_DIR}/VimbaCPP/Bin/${VC_LIB_DIR}"
    )
    set(VIMBA_RT_DIRS ${VIMBA_CPP_RT_DIR} ${VIMBA_C_RT_DIR} "${VIMBA_DIR}/VimbaImageTransform/Bin/${VC_LIB_DIR}")
    message(STATUS "VIMBA_RT_DIRS: ${VIMBA_RT_DIRS}")
else()
    message(FATAL_ERROR "find package not support this platform")
endif()