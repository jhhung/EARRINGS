cmake_minimum_required(VERSION 2.8)
project(Screw)

set(
    project_config
    "${CMAKE_CURRENT_BINARY_DIR}/generated/${PROJECT_NAME}Config.cmake"
)

set(config_install_dir "lib/cmake/${PROJECT_NAME}")
# set(MODULE_INSTALL_DIR "")

include(CMakePackageConfigHelpers)

# Use:
# * PROJECT_NAME
# * MODULE_INSTALL_DIR
configure_package_config_file(
    src/Config.cmake.in
    ${project_config}
    INSTALL_DESTINATION "${config_install_dir}"
)

install(
    DIRECTORY screw 
    DESTINATION ${config_install_dir}
)

install(FILES ${project_config} DESTINATION ${config_install_dir})
