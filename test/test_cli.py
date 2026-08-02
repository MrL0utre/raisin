#!/usr/bin/env python3
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run(binary: Path, root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(binary), *args],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=120,
        check=False,
    )


def metrics(output: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in output.splitlines():
        if ";" in line:
            key, value = line.split(";", 1)
            result[key.strip()] = value.strip()
    return result


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_cli.py RAISIN REPOSITORY", file=sys.stderr)
        return 2

    binary = Path(sys.argv[1]).resolve()
    root = Path(sys.argv[2]).resolve()
    graph = root / "hypergraphs" / "b14.rzn2"
    arch = root / "targets" / "arch0.arch"

    first = run(binary, root, str(graph), "stats", "seed", "7")
    second = run(binary, root, str(graph), "stats", "seed", "7")
    assert first.returncode == 0, first.stderr
    assert first.stdout == second.stdout
    stats = metrics(first.stdout)
    assert stats["seed"] == "7"
    assert stats["#vertices"] == "10124"
    assert stats["#hyperedges"] == "10014"
    assert stats["#reds"] == "245"
    assert stats["pmax"] == "23000"

    invalid = run(
        binary, root, str(graph), "part", "dbfs", "part_number", "256"
    )
    assert invalid.returncode != 0
    assert "Invalid part_number" in invalid.stderr

    unknown = run(binary, root, str(graph), "not-a-mode")
    assert unknown.returncode != 0
    assert "Unknown mode" in unknown.stderr

    help_result = run(binary, root, "help")
    assert help_result.returncode == 0
    assert "dkfmfast" in help_result.stdout
    assert "1..255" in help_result.stdout
    assert "write" not in help_result.stdout

    with tempfile.TemporaryDirectory(prefix="raisin-test-") as temp_dir:
        temp = Path(temp_dir)
        temp_graph = temp / graph.name
        shutil.copyfile(graph, temp_graph)
        args = (
            str(temp_graph),
            "part",
            "dbfs",
            "part_number",
            "4",
            "archfile",
            str(arch),
            "seed",
            "11",
        )
        partition_run = run(binary, root, *args)
        assert partition_run.returncode == 0, partition_run.stderr

        solution = Path(str(temp_graph) + ".part.sol")
        assert solution.is_file()
        assignments = [int(line) for line in solution.read_text().splitlines()]
        assert len(assignments) == 10124
        assert all(0 <= part < 4 for part in assignments)

        original = solution.read_bytes()
        replay = run(binary, root, *args)
        assert replay.returncode == 0, replay.stderr
        assert original == solution.read_bytes()

        truncated = temp / "truncated.sol"
        truncated.write_text("0\n1\n", encoding="ascii")
        invalid_partition = run(
            binary,
            root,
            str(graph),
            "eval",
            "part_number",
            "4",
            "partfile",
            str(truncated),
            "archfile",
            str(arch),
        )
        assert invalid_partition.returncode != 0
        assert "ends before vertex" in invalid_partition.stderr

        out_of_range = temp / "out-of-range.sol"
        out_of_range.write_text(
            "4\n" + "0\n" * 10123,
            encoding="ascii",
        )
        invalid_assignment = run(
            binary,
            root,
            str(graph),
            "eval",
            "part_number",
            "4",
            "partfile",
            str(out_of_range),
            "archfile",
            str(arch),
        )
        assert invalid_assignment.returncode != 0
        assert "outside 0..3" in invalid_assignment.stderr

        invalid_arch = temp / "invalid.arch"
        invalid_arch.write_text("4 1\n100 30 0\n", encoding="ascii")
        malformed_architecture = run(
            binary,
            root,
            str(graph),
            "part",
            "dbfs",
            "part_number",
            "4",
            "archfile",
            str(invalid_arch),
        )
        assert malformed_architecture.returncode != 0
        assert "Invalid architecture connection" in malformed_architecture.stderr

        invalid_graph = temp / "invalid.rzn2"
        invalid_graph.write_text(
            "2 2 1 1 1\n"
            "1 0 2\n"
            "1 1 0 1\n"
            "0 1 0 1\n",
            encoding="ascii",
        )
        malformed_hypergraph = run(binary, root, str(invalid_graph), "stats")
        assert malformed_hypergraph.returncode != 0
        assert "Invalid vertex index" in malformed_hypergraph.stderr

    print("CLI regression checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
