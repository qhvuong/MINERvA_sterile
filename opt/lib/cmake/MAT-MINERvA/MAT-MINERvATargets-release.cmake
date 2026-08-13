#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "MAT-MINERvA" for configuration "Release"
set_property(TARGET MAT-MINERvA APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(MAT-MINERvA PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libMAT-MINERvA.so"
  IMPORTED_SONAME_RELEASE "libMAT-MINERvA.so"
  )

list(APPEND _cmake_import_check_targets MAT-MINERvA )
list(APPEND _cmake_import_check_files_for_MAT-MINERvA "${_IMPORT_PREFIX}/lib/libMAT-MINERvA.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
