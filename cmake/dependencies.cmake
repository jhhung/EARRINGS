# gtest
if(BUILD_TESTS)
    hunter_add_package(GTest)
    find_package(GTest CONFIG REQUIRED)
endif()

# boost 
hunter_add_package(Boost COMPONENTS 
     thread 
     system 
     filesystem
     graph
     program_options
     serialization
     # more boost module goes here
)
find_package(Boost CONFIG COMPONENTS 
     thread 
     system 
     filesystem 
     graph
     program_options
     serialization
     # more boost module goes here
     REQUIRED
)

hunter_add_package(libsimdpp)
find_package(libsimdpp CONFIG REQUIRED)

hunter_add_package(range-v3)
find_package(range-v3 CONFIG REQUIRED)

hunter_add_package(ZLIB)
find_package(ZLIB CONFIG REQUIRED)

include(${SCREW_DIR}/hunter_root.cmake)
