cmake_minimum_required(VERSION 3.1)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(eSNP-Sets)

set(CMAKE_CXX_FLAGS "-g -Wall")
option(march "Use --march." OFF)
if(march)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
# Use C++11
set(CMAKE_CXX_STANDARD 11)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(${PROJECT_NAME} "main.cpp"
    "commander.cpp"
    "misc.cpp"
    "gzstream.cpp"
    "functions.h")
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( ${PROJECT_NAME} PRIVATE ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )




