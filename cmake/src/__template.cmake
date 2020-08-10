screw_path_to_target_name(__screw_rel_src_file __screw_target)
if(THIS_IS_EXE_TARGET)
    screw_add_executable(${__screw_target} ${__screw_src_file})
    screw_add_launch_task(${__screw_target})
else()
    screw_add_library(${__screw_target} ${__screw_src_file})
endif()
set(THIS_IS_EXE_TARGET OFF)
target_include_directories(${__screw_target} PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)