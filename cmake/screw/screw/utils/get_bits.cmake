if(NOT DEFINED screw_get_bits_included)
set(screw_get_bits_included)

macro(screw_get_bits b)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(${b} 64)
    elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
        set(${b} 32)
    endif()
endmacro()

else()
endif()