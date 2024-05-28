cmake_minimum_required(VERSION 3.11)
project(mSWEEP)
include(FetchContent)

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

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${CURRENT_BINARY_DIR}/cmake")

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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/mSWEEP_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/mSWEEP_version.h @ONLY)
## Configure OpenMP if it supported on the system.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/mSWEEP_openmp_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/include/mSWEEP_openmp_config.hpp @ONLY)
## Configure MPI if it's supported on the system.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/mSWEEP_mpi_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/include/mSWEEP_mpi_config.hpp @ONLY)

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
  FetchContent_Declare(bxzstr
    GIT_REPOSITORY    https://github.com/tmaklin/bxzstr.git
    GIT_TAG           v1.1.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/bxzstr"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS        -D ZSTD_FOUND=0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(bxzstr)
  set(CMAKE_BXZSTR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/bxzstr/include)
endif()
include_directories(${CMAKE_BXZSTR_HEADERS})

## Check dependencies and download them if not given
## cxxio
if (DEFINED CMAKE_CXXIO_HEADERS)
  message(STATUS "cxxio headers provided in: ${CMAKE_CXXIO_HEADERS}")
else()
  FetchContent_Declare(cxxio
    GIT_REPOSITORY    https://github.com/tmaklin/cxxio.git
    GIT_TAG           v0.1.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/cxxio"
    BUILD_IN_SOURCE   0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(cxxio)
  set(CMAKE_CXXIO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxio/include)
endif()
include_directories("${CMAKE_CXXIO_HEADERS}")

## cxxargs
if (DEFINED CMAKE_CXXARGS_HEADERS)
  message(STATUS "cxxargs headers provided in: ${CMAKE_CXXARGS_HEADERS}")
else()
  FetchContent_Declare(cxxargs
    GIT_REPOSITORY    https://github.com/tmaklin/cxxargs.git
    GIT_TAG           v1.1.4
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/cxxargs"
    BUILD_IN_SOURCE   0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(cxxargs)
  set(CMAKE_CXXARGS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxargs/include)
endif()
include_directories("${CMAKE_CXXARGS_HEADERS}")

## alignment-writer
if (DEFINED CMAKE_ALIGNMENT_WRITER_HEADERS AND DEFINED CMAKE_ALIGNMENT_WRITER_LIBRARY)
  message(STATUS "alignment-writer headers provided in: ${CMAKE_ALIGNMENT_WRITER_HEADERS}")
  message(STATUS "alignment-writer library provided in: ${CMAKE_ALIGNMENT_WRITER_LIBRARY}")
else()
  FetchContent_Declare(alignment-writer
    GIT_REPOSITORY    https://github.com/tmaklin/alignment-writer.git
    GIT_TAG           v0.5.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/alignment-writer"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS      -D CMAKE_BXZSTR_HEADERS=${CMAKE_BXZSTR_HEADERS}
		    -D CMAKE_CXXARGS_HEADERS=${CMAKE_CXXARGS_HEADERS}
		    -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		    -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
		    -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
		    -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
		    -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(alignment-writer)
  add_dependencies(mSWEEP libalignmentwriter)
  set_target_properties(alignment-writer PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_ALIGNMENT_WRITER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/alignment-writer/include)
  set(CMAKE_ALIGNMENT_WRITER_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/lib/libalignment-writer.a)
endif()
include_directories(${CMAKE_ALIGNMENT_WRITER_HEADERS})
target_link_libraries(mSWEEP ${CMAKE_ALIGNMENT_WRITER_LIBRARY})

## telescope
if (DEFINED CMAKE_TELESCOPE_LIBRARY AND DEFINED CMAKE_TELESCOPE_HEADERS)
  message(STATUS "telescope headers provided in: ${CMAKE_TELESCOPE_HEADERS}")
  message(STATUS "telescope library provided in: ${CMAKE_TELESCOPE_LIBRARY}")
else()
  FetchContent_Declare(telescope
    GIT_REPOSITORY    https://github.com/tmaklin/telescope.git
    GIT_TAG           v0.7.0-prerelease
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/telescope"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/telescope"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS      -D CMAKE_BXZSTR_HEADERS=${CMAKE_BXZSTR_HEADERS}
		    -D CMAKE_CXXARGS_HEADERS=${CMAKE_CXXARGS_HEADERS}
		    -D CMAKE_CXXIO_HEADERS=${CMAKE_CXXIO_HEADERS}
		    -D CMAKE_ALIGNMENT_WRITER_HEADERS=${CMAKE_ALIGNMENT_WRITER_HEADERS}
		    -D CMAKE_ALIGNMENT_WRITER_LIBRARY=${CMAKE_ALIGNMENT_WRITER_LIBRARY}
		    -D CMAKE_BITMAGIC_HEADERS=${CMAKE_CURRENT_SOURCE_DIR}/external/telescope/external/BitMagic-7.12.3/src
		    -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		    -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
		    -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
		    -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
		    -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(telescope)
  add_dependencies(telescope libalignmentwriter)
  add_dependencies(mSWEEP telescope)
  set_target_properties(telescope PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_TELESCOPE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/telescope/include)
  set(CMAKE_TELESCOPE_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/lib/libtelescope.a)
  set(CMAKE_BITMAGIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/telescope/external/BitMagic-7.12.3/src)
endif()
include_directories(${CMAKE_TELESCOPE_HEADERS} ${CMAKE_BITMAGIC_HEADERS})
target_link_libraries(mSWEEP ${CMAKE_TELESCOPE_LIBRARY})

## seamat
if (DEFINED CMAKE_SEAMAT_HEADERS)
  message(STATUS "seamat headers provided in: ${CMAKE_SEAMAT_HEADERS}")
else()
  FetchContent_Declare(seamat
    GIT_REPOSITORY    https://github.com/tmaklin/seamat.git
    GIT_TAG           v0.2.1
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/seamat"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS        -D CMAKE_BUILD_TESTS=0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(seamat)
  set(CMAKE_SEAMAT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/seamat/include)
endif()
include_directories(${CMAKE_SEAMAT_HEADERS})

## rcgpar
if (DEFINED CMAKE_RCGPAR_LIBRARY AND DEFINED CMAKE_RCGPAR_HEADERS AND DEFINED CMAKE_RCGUTILS_LIBRARY)
  message(STATUS "rcgpar headers provided in: ${CMAKE_RCGPAR_HEADERS}")
  message(STATUS "rcgpar library provided in: ${CMAKE_RCGPAR_LIBRARY}")
  message(STATUS "rcgutils library provided in: ${CMAKE_RCGUTILS_LIBRARY}")
else()
  FetchContent_Declare(rcgpar
    GIT_REPOSITORY    https://github.com/tmaklin/rcgpar.git
    GIT_TAG           v1.1.2
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/rcgpar"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/rcgpar"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS      -D CMAKE_ENABLE_MPI_SUPPORT=${MSWEEP_MPI_SUPPORT}
                    -D CMAKE_SEAMAT_HEADERS=${CMAKE_SEAMAT_HEADERS}
		    -D CMAKE_BITMAGIC_HEADERS=${CMAKE_BITMAGIC_HEADERS}
		    -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		    -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
		    -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
		    -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
		    -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(rcgpar)
  add_dependencies(mSWEEP rcgomp)
  set(CMAKE_RCGPAR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/rcgpar/include)
  set(CMAKE_RCGPAR_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/librcgomp.a")
  set(CMAKE_RCGUTILS_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/librcgutils.a")
endif()
target_link_libraries(mSWEEP ${CMAKE_RCGPAR_LIBRARY} ${CMAKE_RCGUTILS_LIBRARY})
include_directories(${CMAKE_RCGPAR_HEADERS})

## mGEMS
if (DEFINED CMAKE_MGEMS_LIBRARY AND DEFINED CMAKE_MGEMS_HEADERS)
  message(STATUS "mGEMS headers provided in: ${CMAKE_MGEMS_HEADERS}")
  message(STATUS "mGEMS library provided in: ${CMAKE_MGEMS_LIBRARY}")
else()
  FetchContent_Declare(mGEMS
    GIT_REPOSITORY    https://github.com/PROBIC/mGEMS.git
    GIT_TAG           v1.3.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/mGEMS"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/mGEMS"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS      -D CMAKE_BXZSTR_HEADERS=${CMAKE_BXZSTR_HEADERS}
		    -D CMAKE_CXXARGS_HEADERS=${CMAKE_CXXARGS_HEADERS}
		    -D CMAKE_CXXIO_HEADERS=${CMAKE_CXXIO_HEADERS}
		    -D CMAKE_ALIGNMENT_WRITER_HEADERS=${CMAKE_ALIGNMENT_WRITER_HEADERS}
		    -D CMAKE_BITMAGIC_HEADERS=${CMAKE_CURRENT_BINARY_DIR}/external/telescope/external/BitMagic-7.12.3/src
		    -D CMAKE_SEAMAT_HEADERS=${CMAKE_SEAMAT_HEADERS}
		    -D CMAKE_TELESCOPE_HEADERS=${CMAKE_TELESCOPE_HEADERS}
		    -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		    -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
		    -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
		    -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
		    -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(mGEMS)
  add_dependencies(mGEMS telescope libalignmentwriter)
  add_dependencies(mSWEEP libmgems)
  set_target_properties(mGEMS PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_MGEMS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/mGEMS/include)
  set(CMAKE_MGEMS_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/external/mGEMS/lib/libmgems.a)
endif()
target_link_libraries(mSWEEP ${CMAKE_MGEMS_LIBRARY} ${CMAKE_TELESCOPE_LIBRARY} ${CMAKE_ALIGNMENT_WRITER_LIBRARY})
include_directories(${CMAKE_MGEMS_HEADERS})

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external ${CMAKE_CURRENT_SOURCE_DIR}/include/tools
  ${CMAKE_CURRENT_BINARY_DIR}/include
  )

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
