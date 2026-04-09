# Progressive Convex Hull Simplification

Research implementation accompanying `convex-hull-simplification.tex`. The code is a working but unpolished prototype.

## What it does

Given a 3D mesh, produces a sequence of progressively simpler convex hulls that are guaranteed to strictly contain the input at every step (conservative / exterior simplification).

The key idea is to work in the *dual* of the convex hull. Removing a vertex from the dual corresponds to removing a halfspace from the primal hull's H-representation — i.e., relaxing one face constraint and slightly enlarging the hull. Vertices are removed greedily in order of the primal volume added (analogous to quadric-error simplification for triangle meshes).

## Build

Requires CMake ≥ 3.16 and a C++17 compiler. Dependencies (libigl, CGAL, SDLP, Eigen, GLFW) are fetched automatically.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --target pchs
```

## Usage

```
./build/pchs [options] [mesh.ply]
```

Defaults to an icosahedron if no mesh is given.

| Option | Description |
|---|---|
| `--target N` | Simplify to N dual vertices / halfspaces (default: 18) |
| `--interactive` | Open a viewer (see below) |
| `--stats` | Print per-phase timing after simplification |

In non-interactive mode the simplified hull is printed to stdout as a polygon mesh in libigl format (`pV`, `pPI`, `pPC`) suitable for pasting into MATLAB/Octave.

### Interactive viewer

```bash
./build/pchs --interactive --target 18 Actaeon.ply
```

The viewer shows the input mesh (gray) and the current simplified hull (orange, 75% transparent).

| Key | Action |
|---|---|
| `space` | Remove one dual vertex (one greedy step) |
| `a` | Simplify to `--target` in one go |
| `r` | Reset to the full convex hull |

## Algorithm sketch

1. Compute the convex hull of the input mesh (primal, exact kernel via CGAL).
2. Find the Chebyshev center of the primal halfspaces via a small linear program (SDLP).
3. Map each non-(nearly-)degenerate primal face to a dual vertex via the polarity transform about the Chebyshev center.
4. Compute the convex hull of the dual points (using exact predicates).
5. For each dual vertex, measure the cost of removing it: the primal volume that would be added to the hull (`primal_volume_subtended`).
6. Maintain a lazy-deletion min-heap. Repeatedly pop the cheapest vertex, erase it from the dual (retriangulating the resulting hole convexly), and update neighbors' costs.
7. Convert the simplified dual back to a primal polygon mesh.

## Library API

The core logic is in `ConvexHullSimplification` (`convex_hull_simplification.h`):

```cpp
ConvexHullSimplification chs(V, F);   // build primal hull, dual, init queue
chs.simplify_to(18);                  // greedily simplify to 18 dual vertices
auto [pV, pPI, pPC] = chs.get_primal_mesh();  // get current hull as polygon mesh
```

For step-by-step control:

```cpp
while(chs.step()) { ... }             // one greedy removal per call
chs.num_dual_vertices();              // current vertex count
chs.stats();                          // timing breakdown (t_primal_hull, t_dual_hull,
                                      //   t_queue_init, t_last_simplify)
```

`MAX_DEGREE_FOR_FLIPS` (env var, default 100) controls the vertex-degree threshold above which the convex-hull-based one-ring triangulation is used instead of the flip-based method.
