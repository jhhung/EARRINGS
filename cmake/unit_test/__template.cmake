screw_path_to_target_name(__screw_rel_src_file __screw_target)
screw_add_executable(${__screw_target} ${__screw_src_file})
target_include_directories(${__screw_target} PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
screw_add_launch_task(${__screw_target})
# Open this comment if you want to use gtest
# target_link_libraries(
#     ${__screw_target} PRIVATE
#     GTest::gtest 
# )
# Open this comment if you need gtest main
# target_link_libraries(
#     ${__screw_target} PRIVATE
#     GTest::main
# )
add_test(
    NAME ${__screw_target}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_INSTALL_PREFIX}/bin/${__screw_target}
)