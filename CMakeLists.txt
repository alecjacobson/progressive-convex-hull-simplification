# libigl's dependency recipes require CMake >= 3.24 (e.g. its Boost recipe).
cmake_minimum_required(VERSION 3.24)
project(progressive-convex-hull-simplification)

# Some fetched dependencies still declare very old cmake_minimum_required()
# values that CMake 4.x rejects. Setting the policy floor here lets those
# subprojects configure without every user having to pass
# -DCMAKE_POLICY_VERSION_MINIMUM=3.5 on the command line. Harmless on CMake < 4.
if(NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
endif()

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(PCHS_INTERACTIVE "Build the interactive polyscope viewer (requires polyscope + embree)" ON)
option(PCHS_PYTHON_BINDINGS "Build Python bindings via nanobind" OFF)
option(PCHS_TESTS "Build the test suite (CTest)" ON)
option(PCHS_NATIVE_REGRESSION "Add a CTest that builds the native backend and A/B-compares vs CGAL (nested build; needs network)" OFF)
set(PCHS_BACKEND "CGAL" CACHE STRING "Mesh backend: CGAL or NATIVE")
# The native backend needs qhull; force it on when selected.
if(PCHS_BACKEND STREQUAL "NATIVE")
  set(PCHS_QHULL ON CACHE BOOL "" FORCE)
endif()
option(PCHS_QHULL "Fetch qhull and build the native backend's convex hull" ON)

include(dependencies)

# Wrapper exposing Shewchuk robust orient3d to native_geom.h. Anything that links
# it gets PCHS_ROBUST_PREDICATES, so nat::orient3d uses the exact predicate.
if(PCHS_QHULL)
  add_library(pchs_predicates_wrap STATIC native_predicates.cpp)
  target_include_directories(pchs_predicates_wrap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(pchs_predicates_wrap PUBLIC pchs_predicates)
  target_compile_definitions(pchs_predicates_wrap PUBLIC PCHS_ROBUST_PREDICATES)
endif()

#
set(PCHS_LIB_SOURCES
  chebyshev_center.cpp
  convex_hull_simplification.cpp
)
if(PCHS_BACKEND STREQUAL "NATIVE")
  list(APPEND PCHS_LIB_SOURCES native_hull.cpp)
else()
  # CGAL-only translation unit (Euler ops, ear-clip path).
  list(APPEND PCHS_LIB_SOURCES vertex_erasure.cpp)
endif()

# The native backend is CGAL-free: link only igl::core + SDLP + qhull.
if(PCHS_BACKEND STREQUAL "NATIVE")
  set(PCHS_MESH_LIBS igl::core SDLP pchs_qhull_r)
else()
  set(PCHS_MESH_LIBS igl::core igl_copyleft::cgal SDLP)
endif()

add_library(progressive_convex_hull_simplification STATIC ${PCHS_LIB_SOURCES})
target_include_directories(progressive_convex_hull_simplification PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(progressive_convex_hull_simplification PUBLIC ${PCHS_MESH_LIBS})
if(PCHS_BACKEND STREQUAL "NATIVE")
  target_compile_definitions(progressive_convex_hull_simplification PUBLIC PCHS_BACKEND_NATIVE)
  if(TARGET pchs_predicates_wrap)
    target_link_libraries(progressive_convex_hull_simplification PUBLIC pchs_predicates_wrap)
  endif()
endif()

add_executable(pchs main.cpp)
target_link_libraries(pchs PUBLIC
  ${PCHS_MESH_LIBS}
  progressive_convex_hull_simplification
)

if(PCHS_INTERACTIVE)
  target_compile_definitions(pchs PRIVATE PCHS_INTERACTIVE)
  target_link_libraries(pchs PUBLIC
    polyscope
    igl::embree
  )
endif()

target_compile_options(pchs PRIVATE
  -Wno-format
  $<$<CONFIG:RelWithDebInfo>:-g>
  $<$<CONFIG:RelWithDebInfo>:-fno-omit-frame-pointer>
  $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
)

if(PCHS_PYTHON_BINDINGS)
  nanobind_add_module(pchs_python_module STABLE_ABI pchs_bindings.cpp)
  set_target_properties(pchs_python_module PROPERTIES OUTPUT_NAME "pchs")
  target_link_libraries(pchs_python_module PRIVATE
    progressive_convex_hull_simplification
  )
  target_compile_options(pchs_python_module PRIVATE -Wno-format)
endif()

if(PCHS_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

find_program(DSYMUTIL_PROGRAM dsymutil)

if(APPLE)
  # For Xcode: ask Xcode itself to emit dSYMs.
  set_target_properties(pchs PROPERTIES
    XCODE_ATTRIBUTE_GCC_GENERATE_DEBUGGING_SYMBOLS[variant=Debug] YES
    XCODE_ATTRIBUTE_GCC_GENERATE_DEBUGGING_SYMBOLS[variant=RelWithDebInfo] YES
    XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT[variant=Debug] "dwarf-with-dsym"
    XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT[variant=RelWithDebInfo] "dwarf-with-dsym"
  )

  # For non-Xcode generators: generate the dSYM explicitly.
  # Still annoying. Seems I need to rm pchs.dSYM everytime
  if(DSYMUTIL_PROGRAM AND NOT CMAKE_GENERATOR MATCHES "Xcode")
    add_custom_command(TARGET pchs POST_BUILD
      COMMAND ${DSYMUTIL_PROGRAM}
              "$<TARGET_FILE:pchs>"
              -o
              "$<TARGET_FILE:pchs>.dSYM"
      BYPRODUCTS "pchs.dSYM"
      VERBATIM
      COMMENT "Generating dSYM bundle for pchs"
    )
  endif()
endif()
