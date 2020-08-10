if(NOT DEFINED screw_get_win_disk_label_included)
set(screw_get_win_disk_label_included)

macro(screw_get_win_disk_label var)
    set(${var}
        C D E F G H I J
    )
endmacro()
macro(screw_get_win_disk_paths var)
    set(${var}
        "C:/"
        "D:/" 
        "E:/"
        "F:/" 
        "G:/"
        "H:/"
        "I:/"
        "J:/"
    )
endmacro()

endif()