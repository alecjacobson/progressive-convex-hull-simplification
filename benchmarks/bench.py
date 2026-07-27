#!/usr/bin/env python3
"""Release-mode benchmark / regression harness for pchs.

Compares three binaries on a set of models by parsing `--stats` phase timings
(median of N runs). Edit BINS/MODELS below to point at your builds and meshes.

  python3 benchmarks/bench.py [runs]        # default 4 runs

BINS: build each in Release, e.g.
  # baseline (some earlier commit), CGAL:
  git worktree add /tmp/base <commit>
  cmake -S /tmp/base -B /tmp/base/build -DCMAKE_BUILD_TYPE=Release -DPCHS_INTERACTIVE=OFF
  cmake --build /tmp/base/build --target pchs
  # HEAD CGAL:
  cmake -B build-headless -DCMAKE_BUILD_TYPE=Release -DPCHS_INTERACTIVE=OFF
  cmake --build build-headless --target pchs
  # HEAD native:
  cmake -B build-native -DCMAKE_BUILD_TYPE=Release -DPCHS_INTERACTIVE=OFF \
        -DPCHS_TESTS=OFF -DPCHS_BACKEND=NATIVE
  cmake --build build-native --target pchs
"""
import subprocess, re, statistics, sys

BINS = {
    "baseline-CGAL": "/tmp/base/build/pchs",
    "HEAD-CGAL":     "build-headless/pchs",
    "HEAD-native":   "build-native/pchs",
}
# (label, mesh path, target). Sphere clouds (from gen_sphere.py) put ~every
# vertex on the hull -> a huge dual -> a simplify-loop-dominated workload.
MODELS = [
    ("Actaeon (118k verts)", "Actaeon.ply", 100),
    ("Pan (250k verts)",     "/tmp/pan_extract/Pan.stl", 100),
    ("sphere-50k",           "/tmp/sphere_50000.ply", 100),
    ("sphere-100k",          "/tmp/sphere_100000.ply", 100),
    ("sphere-200k",          "/tmp/sphere_200000.ply", 100),
]
PHASES = ["primal hull", "dual hull", "queue init", "simplify_to", "total"]


def run(binpath, model, target):
    out = subprocess.run(
        [binpath, "--stats", "--target", str(target), model,
         "--primal-output", "/tmp/bench_p.ply", "--dual-output", "/tmp/bench_d.ply"],
        capture_output=True, text=True, timeout=1200).stdout
    d = {}
    for ph in PHASES:
        m = re.search(rf"^{re.escape(ph)}:\s+([0-9.eE+-]+)\s*s", out, re.M)
        d[ph] = float(m.group(1)) if m else float("nan")
    return d


def median_runs(binpath, model, target, R):
    runs = [run(binpath, model, target) for _ in range(R)]
    return {ph: statistics.median(r[ph] for r in runs) for ph in PHASES}


def main():
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    print(f"Release-mode benchmark, median of {R} runs, seconds\n")
    for name, model, target in MODELS:
        res = {b: median_runs(p, model, target, R) for b, p in BINS.items()}
        print(f"### {name}  (--target {target})")
        print("  phase        " + "".join(f"{b:>16}" for b in BINS))
        for ph in PHASES:
            print(f"  {ph:<12}" + "".join(f"{res[b][ph]:16.4g}" for b in BINS))
        bc, hc, hn = (res[b]["total"] for b in BINS)
        algo = lambda b: res[b]["queue init"] + res[b]["simplify_to"]
        print(f"  total: HEAD-CGAL/baseline {hc/bc:.3f}x   native/HEAD-CGAL {hn/hc:.3f}x")
        print(f"  algorithm-only (queue+simplify): HEAD-CGAL/baseline "
              f"{algo('HEAD-CGAL')/algo('baseline-CGAL'):.3f}x   "
              f"native/CGAL {algo('HEAD-native')/algo('HEAD-CGAL'):.3f}x\n")


if __name__ == "__main__":
    main()
