# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

Two build directories exist: `build/` (release) and `build-debug/` (debug). Use CMake with Ninja or Make.

```bash
# Configure (first time or after CMakeLists.txt changes)
cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo

# Build the executable
cmake --build build-debug --target pchs
cmake --build build --target pchs
```

The executable is `pchs`. Run it with an optional mesh argument (PLY, OBJ, etc.):

```bash
./build/pchs Actaeon.ply
./build/pchs  # defaults to icosahedron
```

`MAX_DEGREE_FOR_FLIPS` can be set via environment variable to control which triangulation method is used for high-degree vertices.

## Dependencies (fetched automatically by CMake via FetchContent)

- **libigl**: geometry processing (mesh I/O, linear programming, min-heap, etc.)
- **CGAL** (via libigl copyleft): exact/inexact kernels, convex hull, polyhedron operations
- **SDLP** (`alecjacobson/SDLP`): small-dimension linear programming used in `chebyshev_center`
- **Eigen**: matrix math

## Architecture

The project implements progressive convex hull simplification — iteratively removing vertices from the *dual* of a convex hull to produce a sequence of simpler convex approximations. The algorithm lives entirely in `main.cpp` plus two helper translation units.

### Primal / Dual duality

The core idea: given an input mesh, compute its convex hull (`primal`), then construct the *dual* convex polyhedron (`dual`) via a polarity transform centered at the Chebyshev center of the primal halfspaces. Simplifying the dual (removing dual vertices) corresponds to adding halfspace constraints to the primal — i.e., cutting off primal vertices to produce a tighter convex approximation.

### Key data types

- `EPolyhedron` — `CGAL::Polyhedron_3<EK, Polyhedron_items_with_id_3>` with exact kernel (`EPECK`). Used for the primal convex hull and Chebyshev center computation.
- `IPolyhedron` — same but with inexact kernel (`EPICK`). Used for the dual during simplification (faster arithmetic).
- Vertex `.id()` fields are used as stable indices throughout (must be kept consistent across operations).

### Simplification loop (`main.cpp:1148–1216`)

1. Compute cost of erasing each dual vertex using `measure_vertex_erasure` (cost = primal volume subtended by the vertex's one-ring).
2. Maintain a lazy-deletion `igl::min_heap` priority queue keyed on `(cost, vertex_id, visit_count)`. Stale entries are discarded by checking the visit counter.
3. Pop the minimum-cost vertex, call `erase_vertex_and_clip_ears` to remove it from the dual and re-triangulate the resulting hole, then recompute costs for affected neighbors.

### One-ring triangulation

When a dual vertex is to be removed, its one-ring (the boundary polygon left after removal) must be triangulated convexly so that the dual remains a valid convex polyhedron:

- **`one_ring_triangulation_convex_via_flips`** (used when `degree ≤ MAX_DEGREE_FOR_FLIPS`): copy the one-ring, then flip non-convex edges via `flip_until_all_interior_edges_are_convex` (stack-based, not the O(N³) restart variant).
- **`one_ring_triangulation_convex_via_convex_hull`** (fallback for high-degree vertices): compute the convex hull of the one-ring neighbors and remove faces on the wrong side of the removed vertex.

The triangulation order is encoded as an *ear-clipping path* (`clip_ears`) so it can be replayed on the actual dual (`erase_vertex_and_clip_ears`) without re-triangulating.

### Helper files

- **`chebyshev_center.{h,cpp}`**: finds the largest inscribed ball center (Chebyshev center) of a set of halfplanes by solving a linear program via SDLP. Input is an `N×4` matrix of `[nx, ny, nz, b]` plane coefficients. Template is explicitly instantiated for the two call-site types.
- **`primal_volume_subtended.{hpp,cpp}`**: computes the signed primal volume and dual volume subtended by a vertex's one-ring (used as the simplification cost). Has both an Eigen matrix overload and a CGAL Polyhedron overload.

### Output

After simplification, the simplified dual is converted back to a primal polygon-mesh representation (`dual_to_primal_mesh`) and printed in MATLAB/libigl format (`pV`, `pPI`, `pPC`) to stdout for inspection.
