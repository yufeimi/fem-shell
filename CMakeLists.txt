cmake_minimum_required(VERSION 2.8.12)
project(shell-element)

set(CMAKE_CXX_STANDARD 14)
set(LIBMESH_INCLUDE_DIR	 "" CACHE PATH "LOCATION OF LIBMESH")

option(BUILD_SHARED_LIBS "build using shared libraries" ON)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

add_subdirectory(src/)
