if(NOT DEFINED screw_generate_launch_included)
set(screw_generate_launch_included)

macro(screw_generate_launch)
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/.vscode/launch.json")
        get_property(__screw_launch_tasks GLOBAL PROPERTY __screw_launch_tasks)
        if(SCREW_DEBUG)
            screw_show_var(__screw_launch_tasks)
        endif()
        set(
            template 
            ${SCREW_DIR}/utils/vscode/launch/launch.json.in 
        )
        if(SCREW_DEBUG)
            screw_show_var(template)
        endif()
        configure_file(
            ${template}
            ${CMAKE_SOURCE_DIR}/.vscode/launch.json
            @ONLY
        )
        unset(template)
    endif()
endmacro()


else()
endif()