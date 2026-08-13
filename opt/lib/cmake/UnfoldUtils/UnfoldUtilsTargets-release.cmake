#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "UnfoldUtils" for configuration "Release"
set_property(TARGET UnfoldUtils APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(UnfoldUtils PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libUnfoldUtils.so"
  IMPORTED_SONAME_RELEASE "libUnfoldUtils.so"
  )

list(APPEND _cmake_import_check_targets UnfoldUtils )
list(APPEND _cmake_import_check_files_for_UnfoldUtils "${_IMPORT_PREFIX}/lib/libUnfoldUtils.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
