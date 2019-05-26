cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_COMPILER "mpicxx")
project(shell-element)

set(CMAKE_CXX_STANDARD 14)
set(LIBMESH_INCLUDE_DIR	 "" CACHE PATH "LOCATION OF LIBMESH")

option(BUILD_SHARED_LIBS "build using shared libraries" ON)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(LIBMESH_SEARCH_PATHS
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local
    /usr
    /sw # Fink
    /opt/local # DarwinPorts
    /opt/csw # Blastwave
    /opt
    ${LIBMESH_DIR})

find_library(libmesh_LIBRARY
    NAMES mesh_opt
    HINTS ${LIBMESH_INCLUDE_DIR}
    PATH_SUFFIXES lib
    PATH ${LIBMESH_SEARCH_PATHS})

enable_testing()
add_subdirectory(source)
add_subdirectory(tests)
