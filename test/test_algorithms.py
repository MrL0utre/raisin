#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only

import json
import subprocess
import sys
import tempfile
from pathlib import Path


def run(binary: Path, root: Path, graph: Path, arguments: list[str]) -> dict[str, str]:
    process = subprocess.run(
        [str(binary), str(graph), *arguments],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=240,
        check=False,
    )
    assert process.returncode == 0, process.stderr
    metrics: dict[str, str] = {}
    for line in process.stdout.splitlines():
        if ";" in line:
            key, value = line.split(";", 1)
            metrics[key.strip()] = value.strip()
    return metrics


def objective(metrics: dict[str, str], name: str) -> int:
    matches = [int(value) for key, value in metrics.items() if key.endswith(name)]
    assert len(matches) == 1, (name, metrics)
    return matches[0]


def validate_solution(path: Path, vertices: int, parts: int) -> None:
    assignments = [int(line) for line in path.read_text(encoding="ascii").splitlines()]
    assert len(assignments) == vertices
    assert all(0 <= assignment < parts for assignment in assignments)


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: test_algorithms.py RAISIN REPOSITORY BASELINE", file=sys.stderr)
        return 2

    binary = Path(sys.argv[1]).resolve()
    root = Path(sys.argv[2]).resolve()
    baseline_path = Path(sys.argv[3]).resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    graph = root / baseline["circuit"]
    arch = root / baseline["architecture"]
    seed = str(baseline["seed"])
    part_count = str(baseline["part_count"])

    with tempfile.TemporaryDirectory(prefix="raisin-algorithms-") as temp_dir:
        temp = Path(temp_dir)
        cases = {
            "part_dbfs": [
                "part", "dbfs", "part_number", part_count,
                "archfile", str(arch), "partfile", str(temp / "dbfs"),
                "seed", seed,
            ],
            "part_ddfs": [
                "part", "ddfs", "part_number", part_count,
                "archfile", str(arch), "partfile", str(temp / "ddfs"),
                "seed", seed,
            ],
            "part_ccp": [
                "part", "ccp", "part_number", part_count,
                "archfile", str(arch), "partfile", str(temp / "ccp"),
                "seed", seed,
            ],
        }

        results: dict[str, dict[str, str]] = {}
        for name, arguments in cases.items():
            results[name] = run(binary, root, graph, arguments)
            validate_solution(temp / f"{name.removeprefix('part_')}.sol", 10124, 4)

        refinement_cases = {
            "refine_kfm": [
                "refine", "kfm", "part_file", str(temp / "dbfs.sol"),
                "part_number", part_count, "archfile", str(arch),
                "partfile", str(temp / "kfm"), "perform", "5", "seed", seed,
            ],
            "refine_dkfm": [
                "refine", "dkfm", "part_file", str(temp / "dbfs.sol"),
                "part_number", part_count, "archfile", str(arch),
                "partfile", str(temp / "dkfm"), "perform", "5", "seed", seed,
            ],
        }
        for name, arguments in refinement_cases.items():
            results[name] = run(binary, root, graph, arguments)
            validate_solution(temp / f"{name.removeprefix('refine_')}.sol", 10124, 4)

        multilevel_cases = {
            "multilevel_kfm": [
                "multilevel", "cluster", "hem", "part", "dbfs",
                "refine", "kfm", "part_number", part_count,
                "archfile", str(arch), "partfile", str(temp / "mlkfm"),
                "perform", "5", "seed", seed,
            ],
            "multilevel_dkfm": [
                "multilevel", "cluster", "bsc", "part", "ddfs",
                "refine", "dkfm", "part_number", part_count,
                "archfile", str(arch), "partfile", str(temp / "mldkfm"),
                "perform", "5", "seed", seed,
            ],
            "multilevel_dkfmfast": [
                "multilevel", "cluster", "hem", "part", "ccp",
                "refine", "dkfmfast", "part_number", part_count,
                "archfile", str(arch), "partfile", str(temp / "mlfast"),
                "perform", "5", "seed", seed,
            ],
        }
        multilevel_outputs = {
            "multilevel_kfm": "mlkfm.sol",
            "multilevel_dkfm": "mldkfm.sol",
            "multilevel_dkfmfast": "mlfast.sol",
        }
        for name, arguments in multilevel_cases.items():
            results[name] = run(binary, root, graph, arguments)
            validate_solution(temp / multilevel_outputs[name], 10124, 4)

        for name, expected in baseline["cases"].items():
            actual = results[name]
            assert objective(actual, "pmax") <= expected["pmax"], name
            assert objective(actual, "cut") <= expected["cut"], name
            assert int(actual["communication signal hops"]) <= expected["signal_hops"], name
            assert actual["communication feasible"] == "yes", name
            assert actual["communication repair moves"] == "0", name

    print("Algorithm regression checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
