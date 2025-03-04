# gtest
if(BUILD_TESTS AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1 )
    find_package(GTest CONFIG REQUIRED)
endif()

# boost 
find_package(Boost CONFIG COMPONENTS 
     thread 
     system 
     filesystem 
     graph
     program_options
     serialization
     iostreams
     # more boost module goes here
     REQUIRED
)

include_directories(/usr/local/include/libsimdpp-2.1)

find_package(range-v3 CONFIG REQUIRED)

find_package(ZLIB REQUIRED)
set(ZLIB_LIBRARIES "/usr/lib/x86_64-linux-gnu/libz.so")

set(BIOVOLTRON_TESTS OFF)
add_subdirectory(${CMAKE_SOURCE_DIR}/submodules/Biovoltron)
include_directories(${CMAKE_SOURCE_DIR}/submodules/Biovoltron/submodules/htslib)
