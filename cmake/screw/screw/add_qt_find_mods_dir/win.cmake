include(${CMAKE_CURRENT_LIST_DIR}/versions.cmake)
include(${SCREW_DIR}/utils/get_win_disk_label.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/get_qt_win_toolset.cmake)
include(${SCREW_DIR}/utils/show_var.cmake)
screw_get_win_disk_label(disk_labels)
foreach(dl ${disk_labels})
    list( APPEND qt_search_paths 
        "${dl}:/"
    )
    screw_show_var(qt_search_paths)
endforeach()
find_file(
    qt_root_dir Qt
    HINTS ${qt_search_paths}
)
screw_show_var(qt_root_dir)
screw_get_qt_win_toolset(winrt)
macro(try_find_qt5_dir qt_subdir_feature )
    screw_show_var(qt_subdir_feature)
    find_path(
        qt_sdk_dir ${qt_subdir_feature}
        HINTS ${qt_root_dir}
    )
    set(qt_sdk_dir ${qt_sdk_dir}/${qt_subdir_feature})
    screw_show_var(qt_sdk_dir)
    if( EXISTS "${qt_sdk_dir}")
        set(Qt5_DIR ${qt_sdk_dir}/lib/cmake/Qt5)
        set(Qt5_RT_DIRS ${qt_sdk_dir}/bin)
        unset(qt_sdk_dir)
        break()
    endif()
    unset(qt_sdk_dir)
endmacro()
foreach(ver ${screw_qt_search_versions})
    try_find_qt5_dir("${ver}/${winrt}")
endforeach()
if( NOT EXISTS "${Qt5_DIR}" )
    foreach(ver ${screw_qt_search_versions})
        try_find_qt5_dir("Qt${ver}/${ver}/${winrt}")
    endforeach()
endif()
unset(winrt)
unset(qt_temp_subdir)
unset(disk_labels)
unset(qt_search_paths)
unset(screw_qt_search_versions)
unset(qt_root_dir)




