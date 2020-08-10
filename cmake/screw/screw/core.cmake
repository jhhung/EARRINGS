set(SCREW_DIR ${CMAKE_CURRENT_LIST_DIR})
file(GLOB_RECURSE screw_mods "${SCREW_DIR}/utils/*.cmake")
foreach(mod ${screw_mods})
    message(STATUS "load screw module: ${mod}")
    include("${mod}")
endforeach(mod)

