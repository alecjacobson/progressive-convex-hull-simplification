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

# Fetch Geometry Central
FetchContent_Declare(
  GeometryCentral
  GIT_REPOSITORY https://github.com/nmwsharp/geometry-central
  GIT_TAG        v1.0.0
)
FetchContent_MakeAvailable(GeometryCentral)
add_library(GeometryCentral::geometry-central ALIAS geometry-central)



# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)
igl_include(copyleft cgal)

#
add_library(progressive_convex_hull_simplification STATIC 
  measure.cpp
  dual_volume.cpp
  convex_hull.cpp
  boundary_matrix.cpp
  cells.cpp
  half_edge_matrix.cpp
  half_edge_next.cpp
  quad_edge.cpp
  #manifold_surface_mesh_to_polygons.cpp
  convex_polytope_from_halfspaces.cpp
  face_halfspaces.cpp
  chebyshev_center.cpp
  volume_from_halfspaces.cpp
  vector_area_and_barycenters.cpp
  robust_dual_hull.cpp
  primal_volume_subtended.cpp
)
target_include_directories(progressive_convex_hull_simplification PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(progressive_convex_hull_simplification PUBLIC 
  igl::core
  igl_copyleft::cgal
  SDLP
  #GeometryCentral::geometry-central
)

# Add your project files
add_executable(basic main.cpp)

# Link igl (and the glfw module) to your project
target_link_libraries(basic PUBLIC 
  igl::glfw
  igl_copyleft::cgal
  SDLP
  progressive_convex_hull_simplification
  )

# Add your project files
add_executable(dynamic_hull dynamic_hull.cpp)

# Link igl (and the glfw module) to your project
target_link_libraries(dynamic_hull PUBLIC 
  igl::glfw
  igl_copyleft::cgal
  igl::imgui
  SDLP
  progressive_convex_hull_simplification
  )

add_executable(manifold manifold.cpp)
target_link_libraries(manifold PUBLIC
  igl::core
  igl_copyleft::cgal
  progressive_convex_hull_simplification
)

target_compile_options(manifold PRIVATE
  -Wno-format
  $<$<CONFIG:RelWithDebInfo>:-g>
  $<$<CONFIG:RelWithDebInfo>:-fno-omit-frame-pointer>
  $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
)

find_program(DSYMUTIL_PROGRAM dsymutil)

if(APPLE)
  # For Xcode: ask Xcode itself to emit dSYMs.
  set_target_properties(manifold PROPERTIES
    XCODE_ATTRIBUTE_GCC_GENERATE_DEBUGGING_SYMBOLS[variant=Debug] YES
    XCODE_ATTRIBUTE_GCC_GENERATE_DEBUGGING_SYMBOLS[variant=RelWithDebInfo] YES
    XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT[variant=Debug] "dwarf-with-dsym"
    XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT[variant=RelWithDebInfo] "dwarf-with-dsym"
  )

  # For non-Xcode generators: generate the dSYM explicitly.
  # Still annoying. Seems I need to rm manifold.dSYM everytime
  if(DSYMUTIL_PROGRAM AND NOT CMAKE_GENERATOR MATCHES "Xcode")
    add_custom_command(TARGET manifold POST_BUILD
      COMMAND ${DSYMUTIL_PROGRAM}
              "$<TARGET_FILE:manifold>"
              -o
              "$<TARGET_FILE:manifold>.dSYM"
      BYPRODUCTS "manifold.dSYM"
      VERBATIM
      COMMENT "Generating dSYM bundle for manifold"
    )
  endif()
endif()
