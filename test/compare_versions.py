#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only

import argparse
import json
import subprocess
import tempfile
from pathlib import Path


def execute(
    binary: Path, root: Path, arguments: list[str], timeout: int
) -> dict[str, object]:
    try:
        process = subprocess.run(
            [str(binary), *arguments],
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        return {
            "returncode": 124,
            "metrics": {},
            "stdout": error.stdout or "",
            "stderr": f"timed out after {timeout} seconds",
        }
    values: dict[str, int | float | str] = {}
    for line in process.stdout.splitlines():
        if ";" not in line:
            continue
        key, raw_value = line.split(";", 1)
        raw_value = raw_value.strip()
        try:
            value: int | float | str = int(raw_value)
        except ValueError:
            try:
                value = float(raw_value)
            except ValueError:
                value = raw_value
        values[key.strip()] = value
    return {
        "returncode": process.returncode,
        "metrics": values,
        "stdout": process.stdout,
        "stderr": process.stderr,
    }


def objective(metrics: dict[str, object], suffix: str) -> object | None:
    for key, value in metrics.items():
        if key == suffix or key.endswith(" " + suffix):
            return value
    return None


def validate_solution(path: Path, vertices: int, parts: int) -> None:
    assignments = [int(line) for line in path.read_text().splitlines()]
    if len(assignments) != vertices:
        raise AssertionError(f"{path}: expected {vertices} assignments")
    if not all(0 <= assignment < parts for assignment in assignments):
        raise AssertionError(f"{path}: assignment outside 0..{parts - 1}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare scientific metrics produced by two RaiSin binaries."
    )
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("repository", type=Path)
    parser.add_argument(
        "--enforce-objectives",
        action="store_true",
        help="fail when candidate pmax/cut/cost is worse than the baseline",
    )
    parser.add_argument("--timeout", type=int, default=180, help="timeout per scenario")
    arguments = parser.parse_args()

    baseline = arguments.baseline.resolve()
    candidate = arguments.candidate.resolve()
    root = arguments.repository.resolve()
    graph = (root / "hypergraphs" / "b14.rzn2").resolve()
    arch = (root / "targets" / "arch0.arch").resolve()

    report: dict[str, object] = {}
    with tempfile.TemporaryDirectory(prefix="raisin-compare-") as temp_dir:
        temp = Path(temp_dir)
        scenarios = {
            "stats": [str(graph), "stats"],
            "cluster_bsc": [
                str(graph), "cluster", "bsc", "size", "200",
                "archfile", str(arch), "partfile", str(temp / "{version}-bsc"),
            ],
            "part_dbfs": [
                str(graph), "part", "dbfs", "part_number", "4",
                "archfile", str(arch), "partfile", str(temp / "{version}-dbfs"),
            ],
            "multilevel": [
                str(graph), "multilevel", "cluster", "bsc", "part", "ddfs",
                "refine", "dkfm", "part_number", "4", "perform", "2",
                "archfile", str(arch), "partfile", str(temp / "{version}-ml"),
            ],
        }

        for name, template in scenarios.items():
            results: dict[str, dict[str, object]] = {}
            for version, binary in (("baseline", baseline), ("candidate", candidate)):
                command = [item.format(version=version) for item in template]
                results[version] = execute(binary, root, command, arguments.timeout)

                if version == "candidate" and results[version]["returncode"] != 0:
                    raise RuntimeError(
                        f"candidate failed in {name} ({results[version]['returncode']}): "
                        f"{str(results[version]['stderr']).strip()}"
                    )

                if results[version]["returncode"] == 0 and name in {"part_dbfs", "multilevel"}:
                    output_base = next(
                        command[index + 1]
                        for index, item in enumerate(command)
                        if item == "partfile"
                    )
                    validate_solution(Path(output_base + ".sol"), 10124, 4)

            baseline_metrics = results["baseline"]["metrics"]
            candidate_metrics = results["candidate"]["metrics"]
            assert isinstance(baseline_metrics, dict)
            assert isinstance(candidate_metrics, dict)

            if name == "stats":
                if results["baseline"]["returncode"] == 0:
                    for key in ("#vertices", "#hyperedges", "#reds", "pmax"):
                        if baseline_metrics.get(key) != candidate_metrics.get(key):
                            raise AssertionError(
                                f"stats mismatch for {key}: "
                                f"{baseline_metrics.get(key)} != {candidate_metrics.get(key)}"
                            )

            deltas: dict[str, object] = {
                "baseline_returncode": results["baseline"]["returncode"],
                "candidate_returncode": results["candidate"]["returncode"],
            }
            if results["baseline"]["returncode"] != 0:
                deltas["baseline_error"] = str(results["baseline"]["stderr"]).strip()
            for metric_name in ("pmax", "cut", "cost", "clusters"):
                old = objective(baseline_metrics, metric_name)
                new = objective(candidate_metrics, metric_name)
                if isinstance(old, (int, float)) and isinstance(new, (int, float)):
                    deltas[metric_name] = {"baseline": old, "candidate": new, "delta": new - old}
                    if arguments.enforce_objectives and metric_name in {"pmax", "cut", "cost"} and new > old:
                        raise AssertionError(f"{name}: {metric_name} regressed ({old} -> {new})")

            report[name] = deltas

    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
