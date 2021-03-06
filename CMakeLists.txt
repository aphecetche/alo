cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

### CMP0025   Compiler id for Apple Clang is now AppleClang.
### CMP0042   MACOSX_RPATH is enabled by default.

FOREACH (p
    CMP0025 # CMake 3.0
    CMP0042 # CMake 3.0
    CMP0028 # CMake 3.0
    )
  IF (POLICY ${p})
    cmake_policy(SET ${p} NEW)
  ENDIF ()
endforeach ()

project(alo)

set(CMAKE_MODULE_PATH "${alo_SOURCE_DIR}/cmake")

find_package(RapidJSON)
if(ROOTSYS)
find_package(ROOT 6.08)
  endif()
find_package(ALIROOT)

include(CMakeALICE)

set(CMAKE_CXX_STANDARD 14)

if(RAPIDJSON_FOUND)
  add_subdirectory(jsonmap)
endif()

if(ROOT_FOUND AND AliRoot_FOUND)
add_subdirectory(dhlt)
add_subdirectory(r23)
add_subdirectory(binmap)
endif()

find_package(Boost COMPONENTS REQUIRED unit_test_framework)
enable_testing()

