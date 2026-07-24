# Test suite. Enabled by PCHS_TESTS (default ON). Wired into CTest.

add_executable(pchs_tests
  run_tests.cpp
  test_invariants.cpp
  test_components.cpp
  test_native_geom.cpp
  test_native_mesh.cpp
)
target_link_libraries(pchs_tests PRIVATE
  progressive_convex_hull_simplification
  igl::core
  igl_copyleft::cgal
)
target_compile_definitions(pchs_tests PRIVATE
  PCHS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
)
target_compile_options(pchs_tests PRIVATE -Wno-format)

# Native convex hull (qhull) — compiled into the test binary when enabled.
if(PCHS_QHULL)
  target_sources(pchs_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../native_hull.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_native_hull.cpp
  )
  target_link_libraries(pchs_tests PRIVATE pchs_qhull_r)
  target_compile_definitions(pchs_tests PRIVATE PCHS_QHULL)
endif()

# One CTest entry per test tier so failures are localized.
add_test(NAME unit.invariants COMMAND pchs_tests invariants_)
add_test(NAME unit.components COMMAND pchs_tests comp_)
add_test(NAME unit.native_geom COMMAND pchs_tests native_geom)
add_test(NAME unit.native_mesh COMMAND pchs_tests native_mesh)
if(PCHS_QHULL)
  add_test(NAME unit.native_hull COMMAND pchs_tests native_hull)
endif()

# Golden regression on the icosahedron (byte-exact for the CGAL backend).
add_test(NAME golden.icosahedron
  COMMAND ${CMAKE_COMMAND} -E env bash
          ${CMAKE_CURRENT_SOURCE_DIR}/run_golden.sh
          $<TARGET_FILE:pchs>
          ${CMAKE_CURRENT_SOURCE_DIR}/golden
          ${CMAKE_CURRENT_BINARY_DIR}/golden-work
)
