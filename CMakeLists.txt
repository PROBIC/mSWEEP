cmake_minimum_required(VERSION 2.8.12)
project(mSWEEP)

## Determine build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
endif()
if(CMAKE_BUILD_TYPE MATCHES Release)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xHost -fp-model fast=2 -no-prec-div -fimf-use-svml=true -ipo -m64")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -xHost -fp-model fast=2 -no-prec-div -fimf-use-svml=true -ipo -m64")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -O3 -xHost -fp-model fast=2 -no-prec-div -fimf-use-svml=true -ipo -m64")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -flto -funroll-loops -m64")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ffast-math -flto -funroll-loops -m64")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ffast-math -flto -funroll-loops -m64")
  endif()
  if(CMAKE_CXX_COMPILER STREQUAL "GNU")
    set(CMAKE_AR  "gcc-ar")
    set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_C_ARCHIVE_FINISH   true)
  endif()
endif()

## Check dependencies
find_package(Threads REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
  set(MSWEEP_OPENMP_SUPPORT 1)
else()
  set(MSWEEP_OPENMP_SUPPORT 0)
endif()
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(ZLIB)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE _BUILD_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message( STATUS "mSWEEP version: ${_BUILD_VERSION}")
else()
  message(STATUS "GIT not found, version numbering disabled")
  set(_BUILD_VERSION "unknown")
endif()
string(TIMESTAMP _BUILD_TIMESTAMP)

## Generate a version.h file containing build version and timestamp
configure_file(include/version.h.in ${CMAKE_BINARY_DIR}/include/version.h @ONLY)
configure_file(include/openmp_config.hpp.in ${CMAKE_BINARY_DIR}/include/openmp_config.hpp @ONLY)

if(NOT ZLIB_FOUND OR CMAKE_BUILD_ZLIB)
  configure_file(CMakeLists.txt.in ${CMAKE_BINARY_DIR}/external/zlib-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/zlib-download )
  if(result)
    message(FATAL_ERROR "CMake step for zlib failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/zlib-download )
  if(result)
    message(FATAL_ERROR "Build step for zlib failed: ${result}")
  endif()
  set(ZLIB ${CMAKE_BINARY_DIR}/lib/libz.a)
else()
  set(ZLIB ZLIB::ZLIB)
endif()

add_subdirectory(${CMAKE_SOURCE_DIR}/external/telescope
                 ${CMAKE_BINARY_DIR}/external/telescope/build)
set_target_properties(telescope PROPERTIES EXCLUDE_FROM_ALL 1)

include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR}/external/zstr ${CMAKE_SOURCE_DIR}/external ${CMAKE_SOURCE_DIR}/include/tools)
file(GLOB SOURCES "src/*.cpp")

add_executable(mSWEEP ${SOURCES})
add_executable(matchfasta ${CMAKE_SOURCE_DIR}/src/tools/matchfasta.cpp)

## Handle compiling zlib from source
if(CMAKE_BUILD_ZLIB)
  add_subdirectory(${CMAKE_BINARY_DIR}/external/zlib
                   ${CMAKE_BINARY_DIR}/external/zlib)
  set_target_properties(example PROPERTIES EXCLUDE_FROM_ALL 1)
  set_target_properties(minigzip PROPERTIES EXCLUDE_FROM_ALL 1)
  set_target_properties(zlib PROPERTIES EXCLUDE_FROM_ALL 1)
  add_dependencies(mSWEEP zlibstatic)
endif()

# Link libraries
target_link_libraries(mSWEEP Threads::Threads ${ZLIB} libtelescope)
if (OPENMP_FOUND)
  target_link_libraries(mSWEEP OpenMP::OpenMP_CXX)
endif()
