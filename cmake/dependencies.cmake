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

find_package(ZLIB REQUIRED)

set(BIOVOLTRON_TESTS OFF)
add_subdirectory(${CMAKE_SOURCE_DIR}/submodules/Biovoltron)
include_directories(${CMAKE_SOURCE_DIR}/submodules/Biovoltron/submodules/htslib)
include_directories(${CMAKE_SOURCE_DIR}/submodules/Biovoltron/submodules/libsimdpp)
include_directories(${CMAKE_SOURCE_DIR}/submodules/Biovoltron/submodules/range-v3)
