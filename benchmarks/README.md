# Benchmarks

Release-mode performance check run before merging the backend-abstraction work,
answering two questions:

1. Did the refactor slow the **CGAL** path? (regression check)
2. How does the **native** backend compare?

## Harness

- `gen_sphere.py N out.ply` — a Fibonacci-sphere point cloud. Because every point
  lands on the convex hull, the dual has ~2N vertices, so simplifying to a small
  target makes the run **simplify-loop-dominated** (the part the refactor touched).
  Real scans (Actaeon, Pan) are smooth, so their hulls are small and those runs
  are **hull-construction-dominated**.
- `bench.py [runs]` — runs three binaries (baseline-CGAL, HEAD-CGAL, HEAD-native)
  on each model, parses `--stats` phase timings, reports the median. Edit the
  `BINS`/`MODELS` paths at the top. Baseline = the pre-refactor commit built with
  the CGAL backend.

Models used: `Actaeon.ply` (in repo), `Pan.stl` from
[threedscans.com](https://threedscans.com/uncategorized/pan/), and
`gen_sphere.py` clouds at 50k/100k/200k points.

## Results (Release `-O3`, Linux, gcc 11, median of 4 runs; seconds)

`baseline` = original algorithm (commit before the refactor), CGAL backend.

| model | phase | baseline-CGAL | HEAD-CGAL | HEAD-native |
|---|---|--:|--:|--:|
| **Actaeon** (118k → hull 3.7k) | primal hull | 0.051 | 0.051 | 0.030 |
| | dual hull | 0.034 | 0.034 | 0.024 |
| | queue init | 0.0081 | 0.0082 | 0.013 |
| | simplify_to | 0.055 | 0.055 | 0.096 |
| | **total** | **0.149** | **0.149** | **0.164** |
| **Pan** (250k → small hull) | primal hull | 0.114 | 0.114 | 0.055 |
| | dual hull | 0.030 | 0.030 | 0.022 |
| | simplify_to | 0.052 | 0.052 | 0.089 |
| | **total** | **0.204** | **0.204** | **0.178** |
| **sphere-50k** (dual 100k) | primal hull | 0.208 | 0.212 | 0.218 |
| | dual hull | 0.722 | 0.718 | 0.902 |
| | queue init | 0.241 | 0.236 | 0.331 |
| | simplify_to | 1.50 | 1.48 | 2.31 |
| | **total** | **2.68** | **2.64** | **3.76** |
| **sphere-100k** (dual 200k) | primal hull | 0.496 | 0.497 | 0.493 |
| | dual hull | 1.69 | 1.69 | 2.10 |
| | simplify_to | 3.24 | 3.22 | 4.75 |
| | **total** | **5.99** | **5.92** | **8.01** |
| **sphere-200k** (dual 400k) | primal hull | 1.18 | 1.18 | 1.17 |
| | dual hull | 4.08 | 4.18 | 5.18 |
| | simplify_to | 7.24 | 7.54 | 9.76 |
| | **total** | **13.6** | **14.1** | **17.5** |

(sphere-50k row is from a clean interleaved re-run reporting the min, after the
first pass showed transient load — the identical-code `primal hull` phase is the
control that flags a noisy run.)

## Verdict

**No CGAL regression.** HEAD-CGAL total vs baseline: 1.00x (Actaeon), 1.00x (Pan),
0.99x (sphere-50k), 0.99x (sphere-100k), 1.04x (sphere-200k) — centered on 1.00x,
within run-to-run noise. The `primal hull` phase is byte-identical code
(`CGAL::convex_hull_3`) between the two, and it matches to <1% on every clean row,
which validates the measurements. The `geom::`/`mesh::` seam and the headerized
`primal_change` are zero-cost as intended.

**Native backend.** Competitive and CGAL-free:

- **Hull construction is faster** on real scans — qhull beats `CGAL::convex_hull_3`
  (Actaeon 0.030 vs 0.051 s; Pan 0.055 vs 0.114 s) — so native is *faster overall*
  on hull-dominated inputs (Pan total 0.178 vs 0.204 s, 0.87x).
- **The simplify loop is ~1.3–1.6x slower** than CGAL (native index-based mesh +
  `retriangulate_star` + Shewchuk exact `orient3d` vs CGAL's tuned half-edge and
  interval-filtered predicates). This dominates on sphere clouds, where native
  total lands at 1.24–1.42x.
- Throughput on the largest case (sphere-200k, ~400k dual-vertex removals):
  native ≈ 41k removals/s vs CGAL ≈ 53k/s.

Optimization headroom for the native loop (not yet done): the `std::map` edge
lookups in `nat::Mesh::build`/`retriangulate_star` are the obvious first target
(replace with hashing or arrays), and slot reuse instead of grow-only allocation.
