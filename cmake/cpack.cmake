if( ENABLE_CPACK )
    if( WIN )
        set(CPACK_GENERATOR "NSIS")
    endif()
    set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
    set(CPACK_PACKAGE_CONTACT "cpack.contact")
    set(CPACK_PACKAGE_VENDOR  "Centrillion.TW")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ProjectTemplate - example")
    set(CPACK_PACKAGE_VERSION ${${PROJECT_NAME}_VERSION})
    set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-v${${PROJECT_NAME}_VERSION}")
    set(CPACK_COMPONENTS_ALL
        ConfigFiles
        TargetFiles
        Library
        Headers
        Runtime
    )
    if( MSVC )
        if (CMAKE_CL_64)
            set (arch_name "x64")
        else ()
            set (arch_name "x86")
        endif ()
    endif()
    if ( "${CPACK_GENERATOR}" STREQUAL "NSIS" )
        set( CPACK_NSIS_DISPLAY_NAME "${PROJECT_NAME} v${${PROJECT_NAME}_VERSION}")
        set( CPACK_NSIS_CONTACT "cpack.nsis.contact")
        set (CPACK_NSIS_COMPRESSOR "/SOLID lzma")
    endif()
    include(CPack)
    # cpack_add_component_group( Development
    #     DISPLAY_NAME "Summit software development kit"
    #     EXPENDED)
    # cpack_add_component ( SDK
    #     DISPLAY_NAME "Summit program, headers and library"
    #     INSTALL_TYPES FullInstall
    #     GROUP Development
    # )
    # 
    # cpack_add_install_type( FullInstall
    #     DISPLAY_NAME "Full install, including header and library"
    # )
endif()