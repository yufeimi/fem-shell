set(TARGET_SRC fem-shell.cpp)

add_library(shell-element ${TARGET_SRC})
target_include_directories(shell-element PUBLIC "${CMAKE_SOURCE_DIR}/source/")
target_include_directories(shell-element PUBLIC "${LIBMESH_INCLUDE_DIR}")
target_link_libraries(shell-element mesh_opt)

file(GLOB public_headers "${CMAKE_SOURCE_DIR}/source/*.h")
set_target_properties(shell-element PROPERTIES PROJECT_PUBLIC_HEADER "${public_headers}")

install(TARGETS shell-element
	LIBRARY DESTINATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
	PUBLIC_HEADER DESTINATION "${CMAKE_BINARY_DIR}/include/shell-element")
