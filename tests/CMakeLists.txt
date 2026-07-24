# Test suite. Enabled by PCHS_TESTS (default ON). Wired into CTest.

add_executable(pchs_tests
  run_tests.cpp
  test_invariants.cpp
  test_components.cpp
  test_native_geom.cpp
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

# One CTest entry per test tier so failures are localized.
add_test(NAME unit.invariants COMMAND pchs_tests invariants_)
add_test(NAME unit.components COMMAND pchs_tests comp_)
add_test(NAME unit.native_geom COMMAND pchs_tests native_)

# Golden regression on the icosahedron (byte-exact for the CGAL backend).
add_test(NAME golden.icosahedron
  COMMAND ${CMAKE_COMMAND} -E env bash
          ${CMAKE_CURRENT_SOURCE_DIR}/run_golden.sh
          $<TARGET_FILE:pchs>
          ${CMAKE_CURRENT_SOURCE_DIR}/golden
          ${CMAKE_CURRENT_BINARY_DIR}/golden-work
)
