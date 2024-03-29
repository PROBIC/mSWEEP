cmake_minimum_required(VERSION 2.8.12)
project(mSWEEP)

## Determine build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g -Wall -Wextra -Wpedantic")
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

if(CMAKE_BUILD_WITH_FLTO)
  cmake_policy(SET CMP0069 NEW)
  set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION 1)
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
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

### MPI
if (CMAKE_ENABLE_MPI_SUPPORT)
  find_package(MPI REQUIRED)
  set(MSWEEP_MPI_SUPPORT 1)
  include_directories(MPI_C_INCLUDE_DIRS)
  if (CMAKE_MPI_MAX_PROCESSES)
    set(MSWEEP_MPI_MAX_PROCESSES ${CMAKE_MPI_MAX_PROCESSES})
  else()
    set(MSWEEP_MPI_MAX_PROCESSES 1024)
  endif()
else()
  set(MSWEEP_MPI_SUPPORT 0)
endif()

set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/mSWEEP_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/mSWEEP_version.h @ONLY)
## Configure OpenMP if it supported on the system.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/mSWEEP_openmp_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/include/mSWEEP_openmp_config.hpp @ONLY)
## Configure MPI if it's supported on the system.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/mSWEEP_mpi_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/include/mSWEEP_mpi_config.hpp @ONLY)

add_executable(mSWEEP ${CMAKE_CURRENT_SOURCE_DIR}/src/mSWEEP.cpp)

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

## bxzstr
if (DEFINED CMAKE_BXZSTR_HEADERS)
  message(STATUS "bxzstr headers provided in: ${CMAKE_BXZSTR_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-bxzstr.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "CMake step for bxzstr failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "Build step for bxzstr failed: ${result}")
  endif()
  set(CMAKE_BXZSTR_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr/include)
endif()
include_directories(${CMAKE_BXZSTR_HEADERS})

## Check dependencies and download them if not given
## cxxio
if (DEFINED CMAKE_CXXIO_HEADERS)
  message(STATUS "cxxio headers provided in: ${CMAKE_CXXIO_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-cxxio.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/cxxio-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/cxxio-download )
  if(result)
    message(FATAL_ERROR "CMake step for cxxio failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/cxxio-download )
  if(result)
    message(FATAL_ERROR "Build step for cxxio failed: ${result}")
  endif()
  set(CMAKE_CXXIO_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/cxxio/include)
endif()
include_directories("${CMAKE_CXXIO_HEADERS}")

## cxxargs
if (DEFINED CMAKE_CXXARGS_HEADERS)
  message(STATUS "cxxargs headers provided in: ${CMAKE_CXXARGS_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-cxxargs.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs-download )
  if(result)
    message(FATAL_ERROR "CMake step for cxxargs failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs-download )
  if(result)
    message(FATAL_ERROR "Build step for cxxargs failed: ${result}")
  endif()
  set(CMAKE_CXXARGS_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs/include)
endif()
include_directories("${CMAKE_CXXARGS_HEADERS}")

## alignment-writer
if (DEFINED CMAKE_ALIGNMENT_WRITER_HEADERS)
  message(STATUS "alignment-writer headers provided in: ${CMAKE_ALIGNMENT_WRITER_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-alignment-writer.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer-download )
  if(result)
    message(FATAL_ERROR "CMake step for alignment-writer failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer-download )
  if(result)
    message(FATAL_ERROR "Build step for alignment-writer failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer
                   ${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer/build)
  set_target_properties(alignment-writer PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_ALIGNMENT_WRITER_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer/include)
endif()
include_directories(${CMAKE_ALIGNMENT_WRITER_HEADERS})

## telescope
if (DEFINED CMAKE_TELESCOPE_LIBRARY AND DEFINED CMAKE_TELESCOPE_HEADERS)
  find_library(TELESCOPE NAMES telescope HINTS ${CMAKE_TELESCOPE_LIBRARY})
  target_link_libraries(mSWEEP ${TELESCOPE})
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-telescope.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/telescope-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/telescope-download )
  if(result)
    message(FATAL_ERROR "CMake step for telescope failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/telescope-download )
  if(result)
    message(FATAL_ERROR "Build step for telescope failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/telescope
                   ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/build)
  set_target_properties(telescope PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_TELESCOPE_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/include)
  set(CMAKE_BITMAGIC_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/external/BitMagic-7.12.3/src)
  get_property(CMAKE_TELESCOPE_LIBRARY TARGET telescope PROPERTY LOCATION)
  target_link_libraries(mSWEEP libtelescope libalignmentwriter)
endif()
include_directories(${CMAKE_TELESCOPE_HEADERS} ${CMAKE_BITMAGIC_HEADERS})

## seamat
if (DEFINED CMAKE_SEAMAT_HEADERS)
  message(STATUS "seamat headers provided in: ${CMAKE_SEAMAT_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-seamat.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/seamat-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/seamat-download )
  if(result)
    message(FATAL_ERROR "CMake step for seamat failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/seamat-download )
  if(result)
    message(FATAL_ERROR "Build step for seamat failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/seamat
                   ${CMAKE_CURRENT_BINARY_DIR}/external/seamat/build)
  target_link_libraries(mSWEEP rcgomp)
  set(CMAKE_SEAMAT_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/seamat/include)
endif()
include_directories(${CMAKE_SEAMAT_HEADERS})

## rcgpar
if (DEFINED CMAKE_RCGPAR_LIBRARY AND DEFINED CMAKE_RCGPAR_HEADERS)
  find_library(RCGPAR NAMES rcgpar HINTS ${CMAKE_RCGPAR_LIBRARY})
  target_link_libraries(mSWEEP ${RCGPAR})
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-rcgpar.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/rcgpar-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/rcgpar-download )
  if(result)
    message(FATAL_ERROR "CMake step for rcgpar failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/rcgpar-download )
  if(result)
    message(FATAL_ERROR "Build step for rcgpar failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/rcgpar
                   ${CMAKE_CURRENT_BINARY_DIR}/external/rcgpar/build)
  target_link_libraries(mSWEEP rcgomp)
  set(CMAKE_RCGPAR_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/rcgpar/include)
endif()
include_directories(${CMAKE_RCGPAR_HEADERS})

## BitMagic
set(CMAKE_BITMAGIC_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/external/BitMagic-7.12.3/src)
include_directories(${CMAKE_BITMAGIC_HEADERS})

## mGEMS
if (DEFINED CMAKE_MGEMS_LIBRARY AND DEFINED CMAKE_MGEMS_HEADERS)
  find_library(MGEMS NAMES mgems HINTS ${CMAKE_MGEMS_LIBRARY})
  target_link_libraries(mSWEEP ${MGEMS})
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-mGEMS.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/mGEMS-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/mGEMS-download )
  if(result)
    message(FATAL_ERROR "CMake step for libmgems failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/mGEMS-download )
  if(result)
    message(FATAL_ERROR "Build step for libmgems failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/mGEMS
                   ${CMAKE_CURRENT_BINARY_DIR}/external/mGEMS/build)
  set_target_properties(mGEMS PROPERTIES EXCLUDE_FROM_ALL 1)
  target_link_libraries(mSWEEP libmgems)
  set(CMAKE_MGEMS_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/mGEMS/include)
endif()
include_directories(${CMAKE_MGEMS_HEADERS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
${CMAKE_CURRENT_SOURCE_DIR}/external ${CMAKE_CURRENT_SOURCE_DIR}/include/tools
${CMAKE_CURRENT_SOURCE_DIR}/external/cxxio
${CMAKE_CURRENT_BINARY_DIR}/include)

if(MPI_FOUND)
  add_dependencies(mSWEEP rcgmpi bigmpi)
  target_link_libraries(mSWEEP rcgmpi ${LIBRARY_OUTPUT_PATH}/libbigmpi.a)
endif()

add_library(libmsweep
${CMAKE_CURRENT_SOURCE_DIR}/src/Sample.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/PlainSample.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/BootstrapSample.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/Reference.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/Grouping.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/OutfileDesignator.cpp)

set_target_properties(libmsweep PROPERTIES OUTPUT_NAME msweep)
target_link_libraries(mSWEEP libmsweep)

# Link libraries
if (OPENMP_FOUND)
  target_link_libraries(mSWEEP OpenMP::OpenMP_CXX)
endif()
