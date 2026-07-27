#!/usr/bin/env bash
# Native-backend regression: build a CGAL-free native pchs and A/B-compare its
# --stats output against the CGAL pchs. In the early-removal range the two
# backends are geometrically identical (they only diverge later, by greedy-order
# drift from qhull vs CGAL dual points), so we compare there and on Actaeon.
#
# Usage: native_regression.sh <cgal-pchs> <source-dir> <work-dir> [cmake]
set -euo pipefail

CGAL_PCHS="$1"
SRC="$2"
WORK="$3"
CMAKE="${4:-cmake}"

mkdir -p "$WORK"
BUILD="$WORK/build-native"
"$CMAKE" -B "$BUILD" -S "$SRC" \
  -DCMAKE_BUILD_TYPE=Release -DPCHS_INTERACTIVE=OFF -DPCHS_TESTS=OFF \
  -DPCHS_BACKEND=NATIVE -DCMAKE_POLICY_VERSION_MINIMUM=3.5 > "$WORK/cfg.log" 2>&1
"$CMAKE" --build "$BUILD" --target pchs -j > "$WORK/build.log" 2>&1
NATIVE_PCHS="$BUILD/pchs"

fail=0
stat() { "$1" --stats "${@:2}" 2>/dev/null | grep -E "final volume|final mw|final area"; }

compare() {  # <label> <extra pchs args...>
  local label="$1"; shift
  local c n
  c=$(stat "$CGAL_PCHS" "$@")
  n=$(stat "$NATIVE_PCHS" "$@")
  if [ "$c" = "$n" ]; then
    echo "ok: $label"
  else
    echo "MISMATCH: $label"
    echo "  CGAL:   $(echo "$c" | tr '\n' ' ')"
    echo "  NATIVE: $(echo "$n" | tr '\n' ' ')"
    fail=1
  fi
}

# Icosahedron (default mesh), early-removal range where backends are identical.
for t in 18 14 10; do compare "icosahedron --target $t" --target "$t"; done
# Actaeon, if present (matches exactly at this target).
if [ -f "$SRC/Actaeon.ply" ]; then
  compare "Actaeon --target 200" --target 200 "$SRC/Actaeon.ply"
fi

exit $fail
