set(THIS_IS_EXE_TARGET ON)
screw_extend_template()

target_link_libraries(${__screw_target} PUBLIC
    Boost::filesystem
    Boost::system
    Boost::program_options
    Boost::serialization
    Boost::iostreams
    pthread
    stdc++fs
    ZLIB::zlib
    biovoltron
)

