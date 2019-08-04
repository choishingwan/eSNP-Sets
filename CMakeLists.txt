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
set(CMAKE_CXX_EXTENSIONS OFF)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
################################
#  Googletest configuration
################################
# Download and unpack googletest at configure time
configure_file(${CMAKE_SOURCE_DIR}/test/CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build")

# The gtest/gmock targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
                        "${gmock_SOURCE_DIR}/include")
endif()

include_directories(${CMAKE_SOURCE_DIR}/inc)
add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/src/commander.cpp
    ${CMAKE_SOURCE_DIR}/src/misc.cpp
    ${CMAKE_SOURCE_DIR}/src/gzstream.cpp
    ${CMAKE_SOURCE_DIR}/inc/functions.h)

include_directories(${CMAKE_SOURCE_DIR}/test/inc)
add_executable(UnitTest
    ${CMAKE_SOURCE_DIR}/test/main.cpp
    ${CMAKE_SOURCE_DIR}/src/commander.cpp
    ${CMAKE_SOURCE_DIR}/src/misc.cpp
    ${CMAKE_SOURCE_DIR}/src/gzstream.cpp
    ${CMAKE_SOURCE_DIR}/inc/functions.h
    ${CMAKE_SOURCE_DIR}/test/src/function_test.cpp
    ${CMAKE_SOURCE_DIR}/test/inc/global.h)


find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( ${PROJECT_NAME} PRIVATE ${ZLIB_LIBRARIES} )
    target_link_libraries( UnitTest PRIVATE ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
target_link_libraries(UnitTest PRIVATE gtest gtest_main)
add_test(NAME unitTest COMMAND UnitTest "${CMAKE_CURRENT_LIST_DIR}/test/data/")



