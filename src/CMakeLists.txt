set(TARGET_SRC fem-shell.cpp)

add_library(shell-element ${TARGET_SRC})
target_include_directories(shell-element PRIVATE "${CMAKE_SOURCE_DIR}/fem-shell/src/")
target_include_directories(shell-element PRIVATE "${LIBMESH_INCLUDE_DIR}")

file(GLOB public_headers "${CMAKE_SOURCE_DIR}/src/*.h")
set_target_properties(shell-element PROPERTIES PROJECT_PUBLIC_HEADER "${public_headers}")

install(TARGETS shell-element
	LIBRARY DESTINATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
	PUBLIC_HEADER DESTINATION "${CMAKE_BINARY_DIR}/include/shell-element")
