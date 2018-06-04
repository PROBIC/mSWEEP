cmake_minimum_required (VERSION 2.8.12)
project (mSWEEP)

## Set C++11 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 11)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++11)
endif()

include_directories(include)
file(GLOB SOURCES "src/*.cpp")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(mSWEEP ${SOURCES})

# Link libraries
target_link_libraries(mSWEEP pthread)

# Configure compiler flags, release enables -O2 optimization.
set(CMAKE_BUILD_TYPE Release)
