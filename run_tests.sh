#!/bin/bash
# run_tests.sh — Test suite for Raisin red-black hypergraph partitioner
# Usage: ./run_tests.sh [hypergraph] [arch]
# Default: hypergraphs/b14.rzn2 targets/arch0.arch

set -euo pipefail

BINARY="./exe/raisin"
HG="${1:-hypergraphs/b14.rzn2}"
ARCH="${2:-targets/arch0.arch}"
PASS=0; FAIL=0

run_test() {
    local name="$1"; shift
    local out rc key
    out=$(timeout 300 "$BINARY" "$HG" "$@" 2>&1)
    rc=$?
    key=$(echo "$out" | grep -E "pmax;|cost;|clusters;|vertices;|cut;" | head -1)
    if [ $rc -eq 0 ]; then
        printf "OK   %-35s %s\n" "$name" "$key"
        PASS=$((PASS+1))
    else
        printf "FAIL %-35s rc=%d  %s\n" "$name" "$rc" "$(echo "$out" | tail -1)"
        FAIL=$((FAIL+1))
    fi
}

echo "=== Raisin test suite ==="
echo "Graph: $HG   Arch: $ARCH"
echo ""

run_test "stats"                 stats
run_test "cluster_hem"           cluster hem size 200 archfile "$ARCH" partfile /tmp/t_hem
run_test "cluster_bsc"           cluster bsc size 200 archfile "$ARCH" partfile /tmp/t_bsc
run_test "part_dbfs"             part dbfs part_number 4 archfile "$ARCH" partfile /tmp/t_dbfs
run_test "part_ddfs"             part ddfs part_number 4 archfile "$ARCH" partfile /tmp/t_ddfs
run_test "part_ccp"              part ccp  part_number 4 archfile "$ARCH" partfile /tmp/t_ccp
run_test "eval"                  eval part_number 4 partfile /tmp/t_dbfs.sol archfile "$ARCH"
run_test "multilevel_hem_dbfs_kfm" \
    multilevel cluster hem part dbfs refine kfm \
    part_number 4 archfile "$ARCH" perform 5 partfile /tmp/t_ml1
run_test "multilevel_bsc_ddfs_dkfm" \
    multilevel cluster bsc part ddfs refine dkfm \
    part_number 4 archfile "$ARCH" perform 5 partfile /tmp/t_ml2

echo ""
echo "Results: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ] && echo "ALL PASS ✓" || { echo "SOME FAILURES"; exit 1; }
