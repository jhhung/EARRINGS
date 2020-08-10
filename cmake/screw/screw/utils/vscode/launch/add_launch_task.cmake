if(NOT DEFINED screw_add_launch_task_included)
set(screw_add_launch_task_included)


macro(screw_add_launch_task target)
    if(NOT SCREW_VS_LAUNCH_TASK_TMP)
        set(SCREW_VS_LAUNCH_TASK_TMP ${CMAKE_BINARY_DIR}/screw_vs_task_tmp.json)
    endif()
    set(__screw_launch_task_name ${target})
    configure_file(
        ${SCREW_DIR}/utils/vscode/launch/task_gdb.in
        ${SCREW_VS_LAUNCH_TASK_TMP}
    )
    file(READ 
        ${SCREW_VS_LAUNCH_TASK_TMP} 
        __screw_launch_task
    )
    get_property(__screw_launch_tasks GLOBAL PROPERTY __screw_launch_tasks)
    if(__screw_launch_tasks)
        set(__screw_launch_tasks "${__screw_launch_tasks},${__screw_launch_task}")
    else()
        set(__screw_launch_tasks "${__screw_launch_task}")
    endif()
    set_property(GLOBAL PROPERTY __screw_launch_tasks "${__screw_launch_tasks}")
endmacro(screw_add_launch_task )


else()
endif()