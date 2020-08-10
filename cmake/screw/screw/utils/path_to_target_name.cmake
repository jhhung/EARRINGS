
if(NOT DEFINED screw_path_to_target_name_included)
set(screw_path_to_target_name_included)

macro(screw_path_to_target_name path_str out)
    string( REPLACE "/" "-" ${out} ${${path_str}})
    string( REPLACE ".cpp" "" ${out} ${${out}})
    # screw_show_var(${out})
endmacro()

endif()