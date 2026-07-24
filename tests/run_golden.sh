#!/usr/bin/env bash
# Golden regression driver: run pchs on the icosahedron for each cost function
# and byte-compare popped-ids / costs / primal PLY against tests/golden/.
#
# Usage: run_golden.sh <pchs-binary> <golden-dir> <work-dir>
set -euo pipefail

PCHS="$1"
GOLDEN="$2"
WORK="$3"
mkdir -p "$WORK"

fail=0
for cf in volume area mean-width; do
  "$PCHS" --target 8 --cost-function "$cf" \
    --primal-output "$WORK/ico-$cf-primal.ply" \
    --dual-output   "$WORK/ico-$cf-dual.ply" \
    --costs         "$WORK/ico-$cf-costs.dmat" \
    --popped-ids    "$WORK/ico-$cf-popped.dmat" > /dev/null

  for kind in primal.ply costs.dmat popped.dmat; do
    # golden names are ico-<cf>-<primal|costs|popped>.<ext>
    base="ico-$cf-${kind%%.*}.${kind#*.}"
    g="$GOLDEN/$base"
    w="$WORK/$base"
    if ! cmp -s "$g" "$w"; then
      echo "MISMATCH: $base"
      fail=1
    else
      echo "ok: $base"
    fi
  done
done

exit $fail
