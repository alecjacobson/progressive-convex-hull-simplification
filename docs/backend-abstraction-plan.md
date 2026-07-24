# Backend Abstraction Plan

Goal: make the geometry/mesh backend of Progressive Convex Hull Simplification
swappable between

- **Backend A — CGAL** (current, working): `CGAL::Polyhedron_3` half-edge mesh,
  CGAL Epick kernel predicates/constructions, `CGAL::convex_hull_3`.
- **Backend B — native**: an `std`-only half-edge surface mesh, Shewchuk robust
  predicates (`orient3d`), and qhull for convex hull.

The switch must be compile-time, both backends must pass the **same** test suite,
and Backend A must remain byte-for-byte unchanged while we build Backend B.

---

## 1. What the algorithm actually needs from a backend

The core algorithm (see `README.md` "Algorithm sketch") lives in templated
headers plus two explicitly-instantiated `.cpp` files. Everything is templated on
a `Polyhedron` type and calls CGAL free functions. There is **one geometric
kernel in play**: `EK == IK == Epick` (the exact-kernel path is commented out in
`convex_hull_simplification.h`), so the backend is effectively *double precision
throughout*. That simplifies Backend B enormously — no exact/inexact kernel
juggling, and `CGAL::Cartesian_converter` collapses to identity.

The backend surface area, extracted from the source, is three layers:

### 1a. Geometry kernel (`geom`)
Types: `FT = double`, `Point3`, `Vector3`.

| Used as (current) | Meaning | Backend B |
|---|---|---|
| `CGAL::orientation(a,b,c,d)` | **orient3d** sign of tetra | Shewchuk `orient3d` |
| `CGAL::orientation(a,b,c,p)` returning `Oriented_side` | side of plane | sign of orient3d |
| `CGAL::cross_product(u,v)` | 3d cross | trivial |
| `CGAL::scalar_product(u,v)` / `u*v` | dot | trivial |
| `CGAL::centroid(a,b,c)` | triangle centroid | trivial |
| `CGAL::squared_area(a,b,c)` | ¼‖(b−a)×(c−a)‖² | trivial |
| `CGAL::squared_distance(a,b)` | ‖a−b‖² | trivial |
| `CGAL::sqrt(x)` | √ | `std::sqrt` |
| `p - q` (pt−pt→vec), `p - ORIGIN`, `ORIGIN + v`, `v/s`, `s*v`, `v+v` | affine ops | operator overloads |
| `CGAL::ORIGIN`, `CGAL::NULL_VECTOR` | constants | constants |
| `CGAL::to_double(x)` | → double | identity |

Only **orient3d** needs a robust predicate; everything else is ordinary double
arithmetic. orient3d is the single predicate that governs mesh convexity
decisions (`edge_side`, flip loop, one-ring containment asserts, `primal_change`
origin test), so robustness there is what matters.

### 1b. Half-edge mesh datastructure (`mesh`)
Element API the templates rely on (must be matched by both backends):

- Iteration: `vertices_begin/end`, `facets_begin/end`, `halfedges_begin/end`,
  `edges_begin/end` (+ const), `size_of_{vertices,facets,halfedges}()`.
- Vertex: `->halfedge()`, `->id()` (mutable), `->point()` (mutable), `->degree()`.
- Halfedge: `->vertex()`, `->next()`, `->prev()`, `->opposite()`, `->facet()`,
  `->id()` (mutable), `->is_border()`, `->is_border_edge()`.
- Facet: `->halfedge()`, `->id()` (mutable).
- Global predicates: `is_pure_triangle()`, `is_closed()`, `reserve(v,h,f)`.
- Nested typedefs used by templates: `Point_3`, `Traits::FT`, `Traits::Point_3`,
  `Vertex_handle`, `Vertex_const_handle`, `Halfedge_handle`,
  `Halfedge_const_handle`, `Facet_(const_)iterator`, etc.

**Euler / mutating operations** (the hard, must-be-correct part):

| Op | Where used | Semantics |
|---|---|---|
| `erase_center_vertex(h)` | `erase_vertex_and_clip_ears` | delete `h`'s target vertex + spokes, merge star into one facet; return a halfedge on that facet |
| `split_facet(h1,h2)` | ear replay | add diagonal `h1→h2` inside a facet, split in two; return new halfedge |
| `erase_facet(h)` | `clip_ears`, hull one-ring cleanup | turn `h`'s facet into a hole (border) |
| `make_triangle(p0,p1,p2)` | fan build | new isolated triangle, return a halfedge |
| `add_vertex_and_facet_to_border(h,g)` | fan build | grow a triangle strip along a border |
| incremental builder (`delegate`) | `extract_copy_of_one_ring` | build fan from polygon soup |
| `CGAL::Euler::flip_edge(e)` | convexifying flips | flip an interior edge |

### 1c. Global algorithms
- `convex_hull_3(begin,end, mesh)` → qhull in Backend B (build half-edge mesh from
  qhull's facets, ensure outward orientation, triangulate).
- `bounding_box(begin,end)` → trivial min/max.
- `polyhedron_chebyshev_center` already reduces to a double `N×4` LP solved by
  **SDLP** — backend-independent, no change.

---

## 2. The seam

Introduce two backend-neutral headers, selected at compile time
(`PCHS_BACKEND=CGAL|NATIVE`):

- `geometry.h` — `namespace geom` with `Point3`, `Vector3`, `FT`, and the
  free functions above. Backend A: thin `inline` forwarders to `CGAL::`.
  Backend B: implementations for the native point/vector types (orient3d via
  Shewchuk).
- `mesh_backend.h` — `using mesh::Polyhedron = <IPolyhedron | NativeMesh>` plus the
  hull/bbox entry points.

The templated algorithm headers (`dual_hull.h`, `convex_triangulation.h`,
`vertex_erasure.*`, `primal_change.*`) are edited to call `geom::…` instead of
`CGAL::…` and to be instantiated on `mesh::Polyhedron`. Because both backends
present the *same member API*, the templates themselves barely change.

`ConvexHullSimplification` keeps its structure; only its concrete
`IPolyhedron`/`EK::Point_3` typedefs and its two `convex_hull_3`/`bounding_box`
call sites move behind `mesh_backend.h`.

---

## 3. Test strategy (build this FIRST)

Three tiers, all wired into CTest. Tiers 1–2 are **backend-independent oracles**
and golden files; tier 3 pins internal components. The suite is written and made
green against Backend A, then reused verbatim to validate Backend B.

### Tier 1 — Global invariant oracles (no golden files; must hold for ANY backend, at every step)
- Dual stays `is_pure_triangle()` and `is_closed()`.
- Euler characteristic V−E+F = 2 on the dual throughout.
- `num_dual_vertices()` decreases by exactly 1 per successful `step()`.
- **Conservative containment**: every primal face plane keeps *all* input
  vertices on the interior side (the hull strictly contains the input) at every
  target.
- Cost monotonicity/sign: reported per-removal costs are ≥ 0; cumulative volume
  is non-decreasing (volume metric).

### Tier 2 — Global golden regression (Backend A: byte-exact; Backend B: tolerance)
- Icosahedron → target 8 for `volume`/`area`/`mean-width`: `popped-ids`,
  `costs.dmat`, `*-primal.ply` match `tests/golden/`.
- Actaeon → target 50/500: determinism + tier-1 invariants + final
  volume/area/mean-width within tolerance of golden scalars.

### Tier 3 — Local unit tests (pin each component; reused across backends)
1. `geom` predicates — orient3d signs on known tetrahedra (incl. coplanar);
   cross/dot/centroid/squared_area/squared_distance vs hand values.
2. Chebyshev center — cube & regular tetrahedron have analytic center+inradius.
3. Dual round-trip — icosahedron primal→dual→primal recovers face planes.
4. Half-edge Euler ops — after each of `erase_center_vertex`, `split_facet`,
   `erase_facet`, `flip_edge`, `make_triangle`,
   `add_vertex_and_facet_to_border`: connectivity invariants + expected
   V/E/F deltas on a small hand-built mesh.
5. One-ring convex triangulation — `via_flips` and `via_convex_hull` both leave
   the removed apex on the negative side of every new face, and produce the same
   triangle set.
6. `clip_ears` / `erase_vertex_and_clip_ears` — recorded path replays to the same
   connectivity that measuring produced.
7. `primal_change` costs — volume/area for one icosahedron removal vs a
   finite-difference oracle (full-hull measure before/after); mean-width vs the
   `ConvexHullSimplification::mean_width()` difference.
8. `mean_width()` — cube/icosahedron vs analytic value.

---

## 4. Roadmap (phased, each phase ends green)

- **Phase 0 — Baseline & tests** ✅ *done*: headless build fixed, deterministic
  golden refs captured (`tests/golden/`), 3-tier test harness (`tests/`) green.
- **Phase A — Seam, no behavior change** ✅ *done*: `geometry.h` (`geom::`) +
  `mesh_backend.h` (`mesh::`) added; the algorithm routes all geometry through
  `geom::` and all datastructure/hull calls through `mesh::`. Default
  `PCHS_BACKEND_CGAL`. Golden outputs byte-identical, all tests green.
- **Phase B — Native backend, isolated** *(in progress)*: implement, unit-tested
  against CGAL:
  - ✅ `native_geom.h` (`nat::`) — double Point3/Vector3 + operators + predicates.
    Landed with a determinant `orient3d` (unit-tested vs CGAL on exact-integer
    inputs, `tests/test_native_geom.cpp`). Note: `nat::orient3d` matches CGAL's
    sign convention (the textbook formula is negated — see the test).
    **Remaining: swap in Shewchuk `orient3d` for near-degenerate robustness.**
  - ✅ `native_mesh.h` (`nat::Mesh`) — std-only, index-based half-edge mesh:
    paired opposites (`opposite(i)==i^1`), skip-deleted handle/iterators with
    CGAL-like navigation, explicit border half-edges, `build()` from a triangle
    soup, `flip_edge()`, `ring_vertices()`, `garbage_collect()`, and — the key
    design win — `retriangulate_star(v, tris)`: erase a vertex's fan and fill
    the hole directly from a triangle list (no ear-clip path replay). Unit +
    integration tested (`tests/test_native_mesh.cpp`), including sequential
    retriangulation of an icosahedron down to a tetrahedron.
  - ✅ `native_hull.{h,cpp}` — `nat::convex_hull_3(points) -> nat::Mesh` via
    qhull's reentrant C API (Qt-triangulated, outward, compacted). Validated vs
    CGAL (vertex count, volume, containment) on cube/icosahedron/random cloud
    (`tests/test_native_hull.cpp`). Built behind `PCHS_QHULL` (default ON),
    fetching only `src/libqhull_r/` so qhull's apps/tests stay out.
  - **Remaining: Shewchuk `orient3d` (robustness hardening, optional for clean
    inputs); `bounding_box` (trivial).**
- **Phase C — Integration & validation**: `PCHS_BACKEND=NATIVE` compiles the
  algorithm on the native mesh; run tiers 1–3 + global goldens (tolerance) on the
  native backend; A/B-compare native vs CGAL on icosahedron + Actaeon.
  Note on the retriangulation splice: rather than reproduce CGAL's
  `erase_center_vertex` + `split_facet`-replay path on `nat::Mesh`, the native
  backend uses `retriangulate_star(v, tris)` directly — the candidate one-ring
  triangulation (as a triangle list over the ring) is spliced into the main mesh
  in one call. So the native `measure_vertex_erasure`/apply step is a thin
  backend-conditional path, not a faithful port of the CGAL ear-clip machinery.

## 5. Risks / watch-list
- **Euler op semantics** must match CGAL exactly (return-halfedge conventions);
  `erase_vertex_and_clip_ears` walks `->next()` from the returned halfedge, and
  `clip_ears` relies on `->opposite()->prev()->opposite()` border walks. Pin with
  tier-3 #4/#6 before trusting integration.
- **Hull face ordering** differs between qhull and CGAL, so Backend B goldens are
  tolerance-based, not byte-exact; tier-1 invariants + scalar comparisons carry
  the correctness burden there.
- **Robustness**: qhull may merge/triangulate coplanar facets differently; the
  `primal_squared_area_tol` degenerate-face filter and `remove_duplicate_vertices`
  dedup must be reproduced.
- `add_vertex_and_facet_to_border` and `erase_center_vertex` are the two CGAL ops
  with the least obvious contracts — implement + test these earliest in Phase B.
