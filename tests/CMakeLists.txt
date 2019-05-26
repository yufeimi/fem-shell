set(shell_tests test_A_uv_t
                test_B_uv_q
                test_C_w_tA16
                test_D_w_q_uni16
                test_E_uvw_t)
              #  test_F_032_ss_uni
              #  test_G_MPI_64_q)

foreach(test ${shell_tests})
    set(meshfile ${CMAKE_SOURCE_DIR}/example-meshes/${test}.xda)
    set(forcefile ${CMAKE_SOURCE_DIR}/example-meshes/${test}_f)
    set(output ${CMAKE_CURRENT_BINARY_DIR}/${test})
    file(MAKE_DIRECTORY ${output})
    add_executable(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.cpp)
    target_include_directories(${test} PRIVATE "${CMAKE_SOURCE_DIR}/source")
    target_include_directories(${test} PRIVATE "${LIBMESH_INCLUDE_DIR}")
    target_link_libraries(${test} shell-element stdc++fs ${libmesh_LIBRARY})
    add_test(NAME ${test} COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test}
    -nu 0.25 -e 30000 -t 1.0 -mesh ${meshfile} -out ${output}/solid WORKING_DIRECTORY ${output})
endforeach()
