cmake_minimum_required (VERSION 2.8.12)
project (mSWEEP)

## Set C++11 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 11)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++11)
endif()

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
configure_file(include/version.h.in ${CMAKE_SOURCE_DIR}/include/version.h @ONLY)

include_directories(include)
file(GLOB SOURCES "src/*.cpp")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(mSWEEP ${SOURCES})

# Link libraries
target_link_libraries(mSWEEP pthread z)

# Configure compiler flags, release enables -O2 optimization.
set(CMAKE_BUILD_TYPE Release)
