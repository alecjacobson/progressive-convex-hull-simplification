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

# Enable CGAL copyleft
igl_include(copyleft cgal)

if(PCHS_INTERACTIVE)
  igl_include(embree)

  FetchContent_Declare(
      polyscope
      GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
      GIT_TAG        v2.3.0
  )
  FetchContent_MakeAvailable(polyscope)
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
