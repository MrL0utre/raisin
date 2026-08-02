#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only

import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path


os.environ.setdefault("PYTHONIOENCODING", "utf-8")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")


def resident_memory_bytes(process_id: int) -> int | None:
    if sys.platform.startswith("linux"):
        try:
            status = Path(f"/proc/{process_id}/status").read_text(encoding="ascii")
        except (FileNotFoundError, PermissionError):
            return None
        for line in status.splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
        return None
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        query = 0x0400
        read = 0x0010
        handle = ctypes.windll.kernel32.OpenProcess(query | read, False, process_id)
        if not handle:
            return None
        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        ok = ctypes.windll.psapi.GetProcessMemoryInfo(
            handle, ctypes.byref(counters), counters.cb
        )
        ctypes.windll.kernel32.CloseHandle(handle)
        return int(counters.WorkingSetSize) if ok else None
    return None


def run(binary: Path, root: Path, graph: Path, arguments: list[str]) -> dict[str, str]:
    with tempfile.TemporaryFile(mode="w+", encoding="utf-8") as stdout_file, (
        tempfile.TemporaryFile(mode="w+", encoding="utf-8")
    ) as stderr_file:
        started = time.perf_counter()
        process = subprocess.Popen(
            [str(binary), str(graph), *arguments],
            cwd=root,
            text=True,
            stdout=stdout_file,
            stderr=stderr_file,
        )
        peak_memory = 0
        deadline = started + 60
        while process.poll() is None:
            sample = resident_memory_bytes(process.pid)
            if sample is not None:
                peak_memory = max(peak_memory, sample)
            if time.perf_counter() >= deadline:
                process.kill()
                process.wait()
                raise subprocess.TimeoutExpired(process.args, 60)
            time.sleep(0.01)
        elapsed = time.perf_counter() - started
        stdout_file.seek(0)
        stderr_file.seek(0)
        stdout = stdout_file.read()
        stderr = stderr_file.read()
    assert process.returncode == 0, stderr
    metrics: dict[str, str] = {}
    for line in stdout.splitlines():
        if ";" in line:
            key, value = line.split(";", 1)
            metrics[key.strip()] = value.strip()
    metrics["__runtime_seconds"] = f"{elapsed:.6f}"
    if peak_memory > 0:
        metrics["__peak_memory_mib"] = f"{peak_memory / (1024 * 1024):.3f}"
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
        clustering_cases = {
            "cluster_hem": [
                "cluster", "hem", "size", "200",
                "archfile", str(arch), "partfile", str(temp / "hem"),
                "seed", seed,
            ],
            "cluster_bsc": [
                "cluster", "bsc", "size", "200",
                "archfile", str(arch), "partfile", str(temp / "bsc"),
                "seed", seed,
            ],
        }
        results: dict[str, dict[str, str]] = {}
        for name, arguments in clustering_cases.items():
            print(f"Running {name}...", flush=True)
            results[name] = run(binary, root, graph, arguments)
            solution = temp / f"{name.removeprefix('cluster_')}.sol"
            assignments = [
                int(line) for line in solution.read_text(encoding="ascii").splitlines()
            ]
            assert len(assignments) == 10124
            assert min(assignments) == 0
            assert max(assignments) + 1 == objective(results[name], "clusters")

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

        for name, arguments in cases.items():
            print(f"Running {name}...", flush=True)
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
            print(f"Running {name}...", flush=True)
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
            print(f"Running {name}...", flush=True)
            results[name] = run(binary, root, graph, arguments)
            validate_solution(temp / multilevel_outputs[name], 10124, 4)

        failures: list[str] = []
        for name, expected in baseline["cases"].items():
            actual = results[name]
            runtime = float(actual["__runtime_seconds"])
            if runtime > expected["runtime_seconds_max"]:
                failures.append(
                    f"{name}: runtime {runtime} exceeds "
                    f"{expected['runtime_seconds_max']} seconds"
                )
            if "__peak_memory_mib" in actual:
                memory = float(actual["__peak_memory_mib"])
                if memory > expected["peak_memory_mib_max"]:
                    failures.append(
                        f"{name}: peak memory {memory} MiB exceeds "
                        f"{expected['peak_memory_mib_max']} MiB"
                    )
            if name.startswith("cluster_"):
                cost = objective(actual, "cost")
                clusters = objective(actual, "clusters")
                if cost > expected["cost"]:
                    failures.append(f"{name}: cost {cost} exceeds {expected['cost']}")
                if clusters != expected["clusters"]:
                    failures.append(
                        f"{name}: clusters {clusters} differs from "
                        f"{expected['clusters']}"
                    )
                continue
            pmax = objective(actual, "pmax")
            cut = objective(actual, "cut")
            signal_hops = int(actual["communication signal hops"])
            balance = float(actual["balance"])
            if pmax > expected["pmax"]:
                failures.append(f"{name}: pmax {pmax} exceeds {expected['pmax']}")
            if cut > expected["cut"]:
                failures.append(f"{name}: cut {cut} exceeds {expected['cut']}")
            if signal_hops > expected["signal_hops"]:
                failures.append(
                    f"{name}: signal hops {signal_hops} exceeds "
                    f"{expected['signal_hops']}"
                )
            if balance > expected["balance_max"]:
                failures.append(
                    f"{name}: balance {balance} exceeds {expected['balance_max']}"
                )
            if actual["communication feasible"] != "yes":
                failures.append(f"{name}: communication is not feasible")
            if actual["communication repair moves"] != "0":
                failures.append(
                    f"{name}: expected no communication repair moves, got "
                    f"{actual['communication repair moves']}"
                )

        if failures:
            raise AssertionError("\n".join(failures))

    print("Algorithm regression checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
