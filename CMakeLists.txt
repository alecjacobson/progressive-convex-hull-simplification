cmake_minimum_required(VERSION 3.16)
project(progressive-convex-hull-simplification)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG 87105143e839217150f062abf6b7a66728bbebb6
)
FetchContent_MakeAvailable(libigl)

include(FetchContent)
FetchContent_Declare(
    SDLP
    GIT_REPOSITORY https://github.com/alecjacobson/SDLP.git
    GIT_TAG        origin/main
)
FetchContent_Populate(SDLP)
# Make an interface target so you can link it normally
add_library(SDLP INTERFACE)
target_include_directories(SDLP INTERFACE ${sdlp_SOURCE_DIR}/include)


# Enable CGAL copyleft
#igl_include(imgui)
igl_include(copyleft cgal)

# Polyscope
include(FetchContent)
FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        v2.3.0
)
FetchContent_MakeAvailable(polyscope)

#
add_library(progressive_convex_hull_simplification STATIC
  chebyshev_center.cpp
  primal_volume_subtended.cpp
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
  polyscope
  igl_copyleft::cgal
  progressive_convex_hull_simplification
)

target_compile_options(pchs PRIVATE
  -Wno-format
  $<$<CONFIG:RelWithDebInfo>:-g>
  $<$<CONFIG:RelWithDebInfo>:-fno-omit-frame-pointer>
  $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
)

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
