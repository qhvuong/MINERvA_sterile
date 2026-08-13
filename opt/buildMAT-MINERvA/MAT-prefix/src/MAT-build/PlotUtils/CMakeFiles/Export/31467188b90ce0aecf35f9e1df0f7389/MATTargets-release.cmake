#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "MAT" for configuration "Release"
set_property(TARGET MAT APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(MAT PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libMAT.so"
  IMPORTED_SONAME_RELEASE "libMAT.so"
  )

list(APPEND _cmake_import_check_targets MAT )
list(APPEND _cmake_import_check_files_for_MAT "${_IMPORT_PREFIX}/lib/libMAT.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
