#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "jsoncpp_static" for configuration "Release"
set_property(TARGET jsoncpp_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(jsoncpp_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/jsoncpp.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS jsoncpp_static )
list(APPEND _IMPORT_CHECK_FILES_FOR_jsoncpp_static "${_IMPORT_PREFIX}/lib/jsoncpp.lib" )

# Import target "jsoncpp_object" for configuration "Release"
set_property(TARGET jsoncpp_object APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(jsoncpp_object PROPERTIES
  IMPORTED_COMMON_LANGUAGE_RUNTIME_RELEASE ""
  IMPORTED_OBJECTS_RELEASE "${_IMPORT_PREFIX}/lib/objects-Release/jsoncpp_object/json_reader.obj;${_IMPORT_PREFIX}/lib/objects-Release/jsoncpp_object/json_value.obj;${_IMPORT_PREFIX}/lib/objects-Release/jsoncpp_object/json_writer.obj"
  )

list(APPEND _IMPORT_CHECK_TARGETS jsoncpp_object )
list(APPEND _IMPORT_CHECK_FILES_FOR_jsoncpp_object "${_IMPORT_PREFIX}/lib/objects-Release/jsoncpp_object/json_reader.obj;${_IMPORT_PREFIX}/lib/objects-Release/jsoncpp_object/json_value.obj;${_IMPORT_PREFIX}/lib/objects-Release/jsoncpp_object/json_writer.obj" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
