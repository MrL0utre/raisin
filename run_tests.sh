#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-only
# Test suite for the RaiSin red-black hypergraph partitioner.
# Usage: ./run_tests.sh [hypergraph] [architecture]

set -euo pipefail

BINARY="${RAISIN_BINARY:-./exe/raisin}"
HG="${1:-hypergraphs/b14.rzn2}"
ARCH="${2:-targets/arch0.arch}"
PASS=0
FAIL=0
RUN_DIR=$(mktemp -d "${TMPDIR:-/tmp}/raisin-tests.XXXXXX")
trap 'rm -rf "$RUN_DIR"' EXIT

run_test() {
    local name="$1"
    shift
    local out rc key

    if out=$(timeout 300 "$BINARY" "$HG" "$@" seed 1 2>&1); then
        rc=0
    else
        rc=$?
    fi
    key=$(echo "$out" | grep -E "pmax;|cost;|clusters;|vertices;|cut;" | head -1 || true)

    if [ "$rc" -eq 0 ]; then
        printf "OK   %-35s %s\n" "$name" "$key"
        PASS=$((PASS + 1))
    else
        printf "FAIL %-35s rc=%d  %s\n" "$name" "$rc" "$(echo "$out" | tail -1)"
        FAIL=$((FAIL + 1))
    fi
}

echo "=== RaiSin test suite ==="
echo "Graph: $HG   Arch: $ARCH"
echo ""

run_test "stats" stats
run_test "cluster_hem" cluster hem size 200 archfile "$ARCH" partfile "$RUN_DIR/t_hem"
run_test "cluster_bsc" cluster bsc size 200 archfile "$ARCH" partfile "$RUN_DIR/t_bsc"
run_test "part_dbfs" part dbfs part_number 4 archfile "$ARCH" partfile "$RUN_DIR/t_dbfs"
run_test "part_ddfs" part ddfs part_number 4 archfile "$ARCH" partfile "$RUN_DIR/t_ddfs"
run_test "part_ccp" part ccp part_number 4 archfile "$ARCH" partfile "$RUN_DIR/t_ccp"
run_test "eval" eval part_number 4 partfile "$RUN_DIR/t_dbfs.sol" archfile "$ARCH"
run_test "multilevel_hem_dbfs_kfm" \
    multilevel cluster hem part dbfs refine kfm \
    part_number 4 archfile "$ARCH" perform 5 partfile "$RUN_DIR/t_ml1"
run_test "multilevel_bsc_ddfs_dkfm" \
    multilevel cluster bsc part ddfs refine dkfm \
    part_number 4 archfile "$ARCH" perform 5 partfile "$RUN_DIR/t_ml2"

echo ""
echo "Results: $PASS passed, $FAIL failed"
if [ "$FAIL" -eq 0 ]; then
    echo "ALL PASS"
else
    echo "SOME FAILURES"
    exit 1
fi
