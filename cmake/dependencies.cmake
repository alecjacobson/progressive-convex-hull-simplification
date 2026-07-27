include(FetchContent)

# Libigl
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG 87105143e839217150f062abf6b7a66728bbebb6
)
FetchContent_MakeAvailable(libigl)

# SDLP
FetchContent_Declare(
    SDLP
    GIT_REPOSITORY https://github.com/alecjacobson/SDLP.git
    GIT_TAG        origin/main
)
FetchContent_Populate(SDLP)
add_library(SDLP INTERFACE)
target_include_directories(SDLP INTERFACE ${sdlp_SOURCE_DIR}/include)

# Enable CGAL copyleft — needed by the CGAL backend and by the test suite
# (the tests validate the native backend against CGAL). A native-only build with
# tests off is fully CGAL-free.
if(NOT PCHS_BACKEND STREQUAL "NATIVE" OR PCHS_TESTS)
  igl_include(copyleft cgal)
endif()

if(PCHS_INTERACTIVE)
  igl_include(embree)

  FetchContent_Declare(
      polyscope
      GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
      GIT_TAG        v2.3.0
  )
  FetchContent_MakeAvailable(polyscope)
endif()

if(PCHS_QHULL)
  # qhull for the native backend's convex hull. We only want the reentrant
  # library; qhull's own CMakeLists pulls in apps + CTest smoketests, so we
  # populate the source but build just src/libqhull_r/ ourselves (SOURCE_SUBDIR
  # points at a dir with no CMakeLists so MakeAvailable skips add_subdirectory).
  FetchContent_Declare(
      qhull
      GIT_REPOSITORY https://github.com/qhull/qhull.git
      GIT_TAG        2020.2
      SOURCE_SUBDIR  src/libqhull_r
  )
  FetchContent_MakeAvailable(qhull)
  file(GLOB PCHS_QHULL_R_SRC ${qhull_SOURCE_DIR}/src/libqhull_r/*.c)
  add_library(pchs_qhull_r STATIC ${PCHS_QHULL_R_SRC})
  target_include_directories(pchs_qhull_r PUBLIC ${qhull_SOURCE_DIR}/src)
  set_target_properties(pchs_qhull_r PROPERTIES POSITION_INDEPENDENT_CODE ON)

  # Shewchuk robust predicates (adaptive orient3d) for the native backend.
  # Precomputed constants (constants.c) -> no runtime exactinit needed. We build
  # just the two source files (SOURCE_SUBDIR points at a header-only dir so
  # MakeAvailable skips the project's own CMakeLists/tests).
  FetchContent_Declare(
      predicates
      GIT_REPOSITORY https://github.com/danshapero/predicates.git
      GIT_TAG        master
      SOURCE_SUBDIR  include
  )
  FetchContent_MakeAvailable(predicates)
  add_library(pchs_predicates STATIC
      ${predicates_SOURCE_DIR}/src/predicates.c
      ${predicates_SOURCE_DIR}/src/constants.c)
  target_include_directories(pchs_predicates PUBLIC ${predicates_SOURCE_DIR}/include)
  set_target_properties(pchs_predicates PROPERTIES POSITION_INDEPENDENT_CODE ON)
  # Shewchuk's exactness relies on error-free transformations; FMA contraction
  # (default on many compilers, esp. arm64/Apple Silicon) would break them.
  target_compile_options(pchs_predicates PRIVATE
    $<$<C_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>)
endif()

if(PCHS_PYTHON_BINDINGS)
  find_package(Python 3.8
    REQUIRED COMPONENTS Interpreter Development.Module
    OPTIONAL_COMPONENTS Development.SABIModule)

  FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG        v2.4.0
  )
  FetchContent_MakeAvailable(nanobind)
endif()
