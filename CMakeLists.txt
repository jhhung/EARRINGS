cmake_minimum_required(VERSION 3.10.0)
include(cmake/HunterGate.cmake)
HunterGate(
    # URL "https://github.com/CHChang810716/hunter.git"
    # SHA1 24c3c526a85ffcc9d34da7d0bed556afd751102f
    # URL "https://github.com/PhoebeWangintw/hunter.git"
    # SHA1 67fadbe7388254c8626b2130c0afb0a5f1817e47
    URL "https://github.com/joyeli/hunter.git"
    SHA1 161d924b5467924d5f960464084f0e0a68a82f4c
    FILEPATH ${CMAKE_CURRENT_LIST_DIR}/cmake/packages.cmake
)
project(EARRINGS LANGUAGES C CXX)

include(cmake/screw/screw/core.cmake)

include(cmake/version.cmake)
include(cmake/options.cmake)
include(cmake/dependencies.cmake)
include(cmake/compiler.cmake)
include(cmake/src.cmake)
include(cmake/unit_test.cmake)
include(cmake/ide.cmake)
include(cmake/install.cmake)
include(cmake/cpack.cmake)
