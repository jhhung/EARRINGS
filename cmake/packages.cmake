if(MSVC)
    set(GTest_CMAKE_ARGS CMAKE_ARGS 
        CMAKE_CXX_FLAGS=/D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
    )
endif()


# hunter_config(
#     Boost
#     VERSION "1.75.0"
# )
hunter_config(
    range-v3
    VERSION "0.3.0-p"
)
hunter_config(
    GTest
    VERSION ${HUNTER_GTest_VERSION}
    ${GTest_CMAKE_ARGS}
)

