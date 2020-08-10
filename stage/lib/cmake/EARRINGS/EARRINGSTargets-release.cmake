#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "EARRINGS::EARRINGS" for configuration "Release"
set_property(TARGET EARRINGS::EARRINGS APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EARRINGS::EARRINGS PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/EARRINGS"
  )

list(APPEND _IMPORT_CHECK_TARGETS EARRINGS::EARRINGS )
list(APPEND _IMPORT_CHECK_FILES_FOR_EARRINGS::EARRINGS "${_IMPORT_PREFIX}/bin/EARRINGS" )

# Import target "EARRINGS::skewer-fastq" for configuration "Release"
set_property(TARGET EARRINGS::skewer-fastq APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EARRINGS::skewer-fastq PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libskewer-fastq.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS EARRINGS::skewer-fastq )
list(APPEND _IMPORT_CHECK_FILES_FOR_EARRINGS::skewer-fastq "${_IMPORT_PREFIX}/lib/libskewer-fastq.a" )

# Import target "EARRINGS::skewer-main" for configuration "Release"
set_property(TARGET EARRINGS::skewer-main APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EARRINGS::skewer-main PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libskewer-main.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS EARRINGS::skewer-main )
list(APPEND _IMPORT_CHECK_FILES_FOR_EARRINGS::skewer-main "${_IMPORT_PREFIX}/lib/libskewer-main.a" )

# Import target "EARRINGS::skewer-matrix" for configuration "Release"
set_property(TARGET EARRINGS::skewer-matrix APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EARRINGS::skewer-matrix PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libskewer-matrix.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS EARRINGS::skewer-matrix )
list(APPEND _IMPORT_CHECK_FILES_FOR_EARRINGS::skewer-matrix "${_IMPORT_PREFIX}/lib/libskewer-matrix.a" )

# Import target "EARRINGS::skewer-parameter" for configuration "Release"
set_property(TARGET EARRINGS::skewer-parameter APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EARRINGS::skewer-parameter PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libskewer-parameter.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS EARRINGS::skewer-parameter )
list(APPEND _IMPORT_CHECK_FILES_FOR_EARRINGS::skewer-parameter "${_IMPORT_PREFIX}/lib/libskewer-parameter.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
