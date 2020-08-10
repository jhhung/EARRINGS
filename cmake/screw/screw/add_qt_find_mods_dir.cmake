# This module defined:
# Qt5_DIR : the cmake path of qt5
# Qt5_RT_DIR : the binrary path of qt5

if(NOT DEFINED screw_add_qt_find_mods_dir_included)
set(screw_add_qt_find_mods_dir_included)

if(WIN32)
    include(${CMAKE_CURRENT_LIST_DIR}/add_qt_find_mods_dir/win.cmake)
elseif(APPLE)
elseif(UNIX)
else()
    message(WARNING "unknown platform for QT path detection.")
endif()

else()
endif()