cmake_minimum_required(VERSION 3.16)
project(progressive-convex-hull-simplification)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(PCHS_INTERACTIVE "Build the interactive polyscope viewer (requires polyscope + embree)" ON)
option(PCHS_PYTHON_BINDINGS "Build Python bindings via nanobind" OFF)
option(PCHS_TESTS "Build the test suite (CTest)" ON)

include(dependencies)

#
add_library(progressive_convex_hull_simplification STATIC
  chebyshev_center.cpp
  primal_change.cpp
  vertex_erasure.cpp
  convex_hull_simplification.cpp
)
target_include_directories(progressive_convex_hull_simplification PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(progressive_convex_hull_simplification PUBLIC
  igl::core
  igl_copyleft::cgal
  SDLP
)

add_executable(pchs main.cpp)
target_link_libraries(pchs PUBLIC
  igl::core
  igl_copyleft::cgal
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
