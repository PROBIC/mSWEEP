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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -funroll-loops")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ffast-math -funroll-loops")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ffast-math -funroll-loops")
  endif()
endif()

## Check dependencies
find_package(OpenMP)
if (OPENMP_FOUND)
  set(MSWEEP_OPENMP_SUPPORT 1)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_PARALLEL")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -D_GLIBCXX_PARALLEL")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -D_GLIBCXX_PARALLEL")
else()
  set(MSWEEP_OPENMP_SUPPORT 0)
endif()
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

add_executable(mSWEEP
${CMAKE_SOURCE_DIR}/src/BootstrapSample.cpp
${CMAKE_SOURCE_DIR}/src/Reference.cpp
${CMAKE_SOURCE_DIR}/src/Sample.cpp
${CMAKE_SOURCE_DIR}/src/likelihood.cpp
${CMAKE_SOURCE_DIR}/src/main.cpp
${CMAKE_SOURCE_DIR}/src/matrix.cpp
${CMAKE_SOURCE_DIR}/src/parse_arguments.cpp
${CMAKE_SOURCE_DIR}/src/process_reads.cpp
${CMAKE_SOURCE_DIR}/src/rcg.cpp
${CMAKE_SOURCE_DIR}/src/read_bitfield.cpp)

## Check supported compression types
find_package(BZip2)
if (BZIP2_FOUND)
  include_directories(${BZIP2_INCLUDE_DIRS})
  target_link_libraries(mSWEEP ${BZIP2_LIBRARIES})
endif()    
find_package(LibLZMA)
if (LIBLZMA_FOUND)
  include_directories(${LIBLZMA_INCLUDE_DIRS})
  target_link_libraries(mSWEEP ${LIBLZMA_LIBRARIES})
endif()
find_package(ZLIB)
if (ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(mSWEEP ${ZLIB_LIBRARIES})
endif()

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

## Check dependencies and download them if not given
## telescope
if (DEFINED CMAKE_TELESCOPE_LIBRARY AND DEFINED CMAKE_TELESCOPE_HEADERS)
  find_library(TELESCOPE NAMES telescope HINTS ${CMAKE_TELESCOPE_LIBRARY})
  target_link_libraries(mSWEEP ${TELESCOPE})
  include_directories("${CMAKE_TELESCOPE_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-telescope.txt.in ${CMAKE_BINARY_DIR}/external/telescope-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/telescope-download )
  if(result)
    message(FATAL_ERROR "CMake step for telescope failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/telescope-download )
  if(result)
    message(FATAL_ERROR "Build step for telescope failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/telescope
                   ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/build)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/telescope/include)
  set_target_properties(telescope PROPERTIES EXCLUDE_FROM_ALL 1)
  target_link_libraries(mSWEEP libtelescope)
endif()

## bxzstr
if (DEFINED CMAKE_BXZSTR_LIBRARY AND DEFINED CMAKE_BXZSTR_HEADERS)
  find_library(BXZSTR NAMES telescope HINTS ${CMAKE_BXZSTR_LIBRARY})
  include_directories("${CMAKE_BXZSTR_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-bxzstr.txt.in ${CMAKE_BINARY_DIR}/external/bxzstr-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "CMake step for bxzstr failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "Build step for bxzstr failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr
                   ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr/build)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr/include)
endif()

include_directories(${CMAKE_SOURCE_DIR}/include
${CMAKE_SOURCE_DIR}/external ${CMAKE_SOURCE_DIR}/include/tools
${CMAKE_SOURCE_DIR}/external/cxxio
${CMAKE_CURRENT_BINARY_DIR}/include)

## Handle compiling zlib from source
if(CMAKE_BUILD_ZLIB)
  add_subdirectory(${CMAKE_BINARY_DIR}/external/zlib
                   ${CMAKE_BINARY_DIR}/external/zlib)
  set_target_properties(example PROPERTIES EXCLUDE_FROM_ALL 1)
  set_target_properties(minigzip PROPERTIES EXCLUDE_FROM_ALL 1)
  set_target_properties(zlib PROPERTIES EXCLUDE_FROM_ALL 1)
  add_dependencies(mSWEEP zlibstatic)
endif()

add_library(msweeptools ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/matchfasta.cpp)
add_executable(matchfasta ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main.cpp)

# Link libraries
target_link_libraries(mSWEEP ${ZLIB} msweeptools)
target_link_libraries(matchfasta msweeptools)
if (OPENMP_FOUND)
  target_link_libraries(mSWEEP OpenMP::OpenMP_CXX)
  target_link_libraries(matchfasta OpenMP::OpenMP_CXX)
  target_link_libraries(msweeptools OpenMP::OpenMP_CXX)
endif()
