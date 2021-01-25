set(THIS_IS_EXE_TARGET ON)
screw_extend_template()

add_library(skewer
    ${CMAKE_SOURCE_DIR}/src/skewer/fastq.cpp
    ${CMAKE_SOURCE_DIR}/src/skewer/main.cpp
    ${CMAKE_SOURCE_DIR}/src/skewer/matrix.cpp
    ${CMAKE_SOURCE_DIR}/src/skewer/parameter.cpp
)

target_link_libraries(${__screw_target} PUBLIC
    Boost::filesystem
    Boost::system
    Boost::program_options
    Boost::serialization
    Boost::iostreams
    pthread
    stdc++fs
    skewer
    ZLIB::zlib
)

