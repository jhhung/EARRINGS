if(NOT DEFINED screw_generate_c_cpp_properties_included)
set(screw_generate_c_cpp_properties_included)

macro(screw_generate_c_cpp_properties template)
    set(C_CPP_PROP "${CMAKE_SOURCE_DIR}/.vscode/c_cpp_properties.json")
    if(NOT EXISTS ${C_CPP_PROP})
        configure_file(
            ${template}
            ${C_CPP_PROP}
            @ONLY
        )
    endif()
endmacro(screw_generate_c_cpp_properties)


else()
endif()
