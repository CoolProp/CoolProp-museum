cmake_minimum_required (VERSION 2.6)
project (CoolProp)

# The version number.
set (CoolProp_VERSION_MAJOR 5)
set (CoolProp_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/CoolPropConfig.h.in"
  "${PROJECT_BINARY_DIR}/CoolPropConfig.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

# add the executable
add_executable(CoolProp CoolProp.cpp)
