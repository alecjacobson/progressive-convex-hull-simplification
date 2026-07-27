# Progressive Convex Hull Simplification

![Progressive convex hull simplification of the Actaeon mesh](pchs-actaeon.gif)

Implementation of the algorithm from the paper

https://arxiv.org/abs/2604.14468

## What it does

Given a 3D mesh, produces a sequence of progressively simpler convex hulls that are guaranteed to strictly contain the input at every step (conservative / exterior simplification).

The key idea is to work in the *dual* of the convex hull. Removing a vertex from the dual corresponds to removing a halfspace from the primal hull's H-representation — i.e., eliminating one face constraint and slightly enlarging the hull. Dual vertices are removed greedily in order of a per-vertex cost (primal volume added, surface area change, or mean width change), analogous to quadric-error simplification for triangle meshes.

## Build

Requires CMake ≥ 3.16 and a C++17 compiler. Dependencies are fetched automatically via FetchContent.

```bash
# Full build (includes interactive viewer; polyscope + embree fetched automatically)
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target pchs

# Headless build (no polyscope or embree dependency)
cmake -B build-headless -DCMAKE_BUILD_TYPE=Release -DPCHS_INTERACTIVE=OFF
cmake --build build-headless --target pchs

# Python bindings (nanobind fetched automatically; pass your Python executable)
cmake -B build-py -DCMAKE_BUILD_TYPE=Release -DPCHS_INTERACTIVE=OFF -DPCHS_PYTHON_BINDINGS=ON \
      -DPython_EXECUTABLE=$(python3 -c "import sys; print(sys.executable)")
cmake --build build-py --target pchs_python_module
```

### Mesh backend (CGAL or native)

The geometry/mesh backend is swappable at configure time via `-DPCHS_BACKEND`:

- `CGAL` (default) — `CGAL::Polyhedron_3` + CGAL predicates + `CGAL::convex_hull_3`.
- `NATIVE` — a dependency-light backend using only the C++ standard library, a
  small half-edge mesh, [qhull](http://www.qhull.org) for the convex hull, and
  [Shewchuk's robust predicates](https://www.cs.cmu.edu/~quake/robust.html) for
  `orient3d`. The native build fetches no CGAL and links no gmp/mpfr:

  ```bash
  cmake -B build-native -DCMAKE_BUILD_TYPE=Release -DPCHS_INTERACTIVE=OFF \
        -DPCHS_TESTS=OFF -DPCHS_BACKEND=NATIVE
  cmake --build build-native --target pchs
  ```

Both backends produce the same hulls; they differ only in floating-point tie-breaks
that can change the greedy removal order late in a run. See
`docs/backend-abstraction-plan.md` for the design and `tests/native_regression.sh`
(also wired as an opt-in CTest via `-DPCHS_NATIVE_REGRESSION=ON`) for the A/B check.

## Usage

```
./build/pchs [options] [mesh.ply]
```

Defaults to an icosahedron if no mesh is given.

| Option | Description |
|---|---|
| `--target N` | Simplify to N dual vertices / halfspaces (default: 18); if N < 0, subtracted from the initial count (e.g. `-1` removes one face) |
| `--cost-function F` | Cost metric: `volume`, `area`, or `mean-width` (default: `volume`) |
| `--primal-output file.ply` | Write simplified hull as a polygon PLY (default: `<stem>-primal.ply`) |
| `--dual-output file` | Write simplified dual as a triangle mesh in any format igl supports (default: `<stem>-dual.ply`) |
| `--primal-initial file.ply` | Write the initial (unsimplified) primal hull |
| `--dual-initial file` | Write the initial (unsimplified) dual hull |
| `--costs [file.dmat]` | Write per-vertex removal costs; NaN for surviving vertices (default: `<stem>-costs.dmat`) |
| `--popped-ids [file.dmat]` | Write vertex removal order (default: `<stem>-popped-ids.dmat`) |
| `--stats` | Print per-phase timing after simplification |
| `--interactive` | Open a viewer (see below; requires `PCHS_INTERACTIVE=ON` build) |
| `--help`, `-h` | Print usage |

Default output filenames are derived from the input stem, e.g. `Actaeon.ply` → `Actaeon-primal.ply`, `Actaeon-dual.ply`, `Actaeon-costs.dmat`, etc.

### Interactive viewer

```bash
./build/pchs --interactive --target 18 Actaeon.ply
```

Opens a [polyscope](https://polyscope.run) window showing the input mesh (ambient occlusion computed in the background) and the current simplified hull (transparent, colored by face normals or graph coloring). Controls are in the left panel under **Hull Simplification**:

- **Target** / **Animate to target** / **Step** / **Go to target** — drive simplification one step at a time or animate to the target count.
- **Cost function** — switch between `Volume`, `Area`, and `Mean Width` cost metrics; switching resets to the full hull and re-simplifies to the current vertex count.
- **Hull alpha** — adjust hull transparency.
- **Coloring** — graph coloring or normal pseudocolor.
- **Reset** — rebuild from the full hull.

## Algorithm sketch

1. Compute the convex hull of the input primal mesh (CGAL).
2. Find the Chebyshev center of the primal halfspaces via a small linear program (SDLP).
3. Map each non-(nearly-)degenerate primal face to a dual vertex via the polarity transform about the Chebyshev center.
4. Compute the convex hull of the dual points (CGAL).
5. For each dual vertex, measure the cost of removing it and re-triangulating the hole convexly: the primal volume added, the change in primal surface area, or the change in mean width.
6. Maintain a lazy-deletion min-heap. Repeatedly pop the cheapest vertex, erase it from the dual (retriangulating the resulting hole convexly), and update neighbors' costs.
7. Convert the simplified dual back to a primal polygon mesh.

## Library API

The core logic is in `ConvexHullSimplification` (`convex_hull_simplification.h`):

```cpp
ConvexHullSimplification chs(V, F);   // build primal hull, dual, init queue
chs.simplify_to(18);                  // greedily simplify to 18 dual vertices
auto [pV, pPI, pPC] = chs.get_primal_mesh();  // get current hull as polygon mesh
auto [dV, dF]       = chs.get_dual_mesh();    // get current dual as triangle mesh
```

For step-by-step control:

```cpp
while(chs.step()) { ... }             // one greedy removal per call
chs.num_dual_vertices();              // current vertex count
chs.popped_dual_vertex_costs();       // VectorXd, NaN for surviving vertices
chs.popped_dual_vertex_ids();         // removal order
chs.stats();                          // timing breakdown (t_primal_hull, t_dual_hull,
                                      //   t_queue_init, t_last_simplify)
```

The cost function defaults to `CostFunction::PRIMAL_VOLUME`; alternatives are `CostFunction::PRIMAL_AREA` and `CostFunction::PRIMAL_MEAN_WIDTH`:

```cpp
ConvexHullSimplification chs(V, F, /*max_degree_for_flips=*/100, CostFunction::PRIMAL_AREA);
```

`MAX_DEGREE_FOR_FLIPS` (env var, default 100) controls the vertex-degree threshold above which the convex-hull-based one-ring triangulation is used instead of the flip-based method.


## Python API

After building with `-DPCHS_PYTHON_BINDINGS=ON`, add the build directory to `PYTHONPATH` and import `pchs`:

```python
import pchs
import igl

V, F = igl.read_triangle_mesh("mesh.obj")

chs = pchs.ConvexHullSimplification(V, F,
    max_degree_for_flips=100,
    cost_function=pchs.CostFunction.volume)  # or pchs.CostFunction.area

chs.simplify_to(18)

pV, pPI, pPC = chs.get_primal_mesh()   # polygon mesh
dV, dF       = chs.get_dual_mesh()     # triangle mesh

costs = chs.popped_dual_vertex_costs() # VectorXd, NaN for surviving vertices
ids   = chs.popped_dual_vertex_ids()   # removal order

s = chs.stats()
print(s.t_primal_hull, s.t_dual_hull, s.t_queue_init, s.t_last_simplify)

# Or use the one-shot wrapper:
pV, pPI, pPC = pchs.simplify_convex_hull(V, F, 18)
```

## LLM use

I wrote the original prototype of this algorithm in matlab without much input
from llms or coding agents (probably copilot tab-complete helped a little bit).
Then I wrote most of first draft of the paper. My thought was that Claude Code
etc. would have no trouble taking the paper source and the matlab prototype and
spitting out a C++ implementation. A few failed attempts and a lot of spent
tokens later, I gave up doing it completely automatically. The llms were
simultaneously immensely helpful parsing through CGAL's otherwise forboding
documentation and templating and also terribly confused by how to build anything
more complicated than simple loops over elements using CGAL's half-edge
datastructure. Using webchat llms effectively as a search engine on the CGAL
documentation, I wrote the initial core of the `vertex_erasure`,
`primal_change`, and `convex_hull_simplification` functions more or less
manually (with tab-complete copilot stuff). When a very minimal version of this
was working, I let Claude Code at it with CMake build junk, adding visualization, flags, python
bindings, and most of this readme file. I did not use coding agents of llms to
write any of the paper text (beyond looking for typos or double-checking
derivations).

I added the mean-width cost metric after the first round of reviews. Again, I
thought that because the area and volume metrics were already in place this
would be an easy fully automatic task for Claude Code. It burnt an entire
session quota and output nothing. So, I ended up writing that manually (with
copilot tab-complete), and then Claude code cleaned up the flags and
visualization code to expose the new cost metric.
