# RaiSin

This repository is a maintained fork of the original
[CEA-LIST/raisin](https://github.com/CEA-LIST/raisin) project. It preserves the
CEA/Inria copyright notices and records the CEA-LIST repository as `upstream`.

RaiSin is a research partitioner for placing synchronous circuits on multi-FPGA
architectures. It represents a circuit as a directed red-black hypergraph:

- red vertices model registers and external I/O ports;
- black vertices model combinational logic;
- directed hyperedges model nets, with one source and one or more sinks;
- vertex weight vectors model target resource consumption;
- target architectures define part capacities and inter-part delays.

The primary objective is to reduce the longest red-to-red path after placement.
RaiSin also considers cut hyperedges and resource balance. The implementation is
based on the algorithms described in Julien Rodriguez's 2024 PhD thesis and the
ICCS 2023/2024 papers listed below.

For placed partitions, RaiSin routes hyperedge signals over the target topology,
reports physical-link loads, and rejects solutions that exceed a declared link
capacity. Architectures may also declare heterogeneous capacities for every
FPGA and resource dimension; these are checked against the circuit's vertex
weight vectors.

## Algorithms

RaiSin provides the following stages:

| Stage | Algorithms |
| --- | --- |
| Clustering | HEM (heavy-edge matching), BSC (binary-search clustering) |
| Initial partitioning | DBFS, DDFS, CCP |
| Refinement | KFM, DKFM |
| Multilevel refinement | KFM, DKFM, DKFMFAST |

The `multilevel` mode combines clustering, initial partitioning and refinement.

## Requirements

- a C11 compiler (GCC or Clang is tested in CI);
- GNU Make, or CMake 3.16 or newer;
- Python 3 for CLI regression tests;
- a POSIX shell and `timeout` for the extended `run_tests.sh` suite.

## Build

With GNU Make:

```sh
make
./exe/raisin help
```

Useful Make targets are `debug`, `sanitize`, `check` and `clean`. Build products
are isolated under `build/` and `exe/`.

With CMake:

```sh
cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build-cmake --parallel
ctest --test-dir build-cmake --output-on-failure
```

To enable AddressSanitizer and UndefinedBehaviorSanitizer with GCC or Clang:

```sh
cmake -S . -B build-sanitize \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DBUILD_TESTING=ON \
  -DRAISIN_ENABLE_SANITIZERS=ON
cmake --build build-sanitize --parallel
ctest --test-dir build-sanitize --output-on-failure
```

Warnings are treated as errors by default. Set
`-DRAISIN_WARNINGS_AS_ERRORS=OFF` only for unsupported toolchains.

## Command line

The general form is:

```text
raisin <graph.rzn2> <mode> <mode arguments> [options]
```

Run `raisin help` for the complete syntax. Common options are:

- `archfile <path>`: target architecture; defaults to `targets/arch0.arch`;
- `partfile <prefix>`: output prefix; RaiSin appends `.sol`;
- `bfactor <n>`: balance factor, default `5`;
- `seed <n>`: pseudo-random seed, default `1`.

RaiSin uses an 8-bit partition identifier by default, so `part_number` must be
between 1 and 255. Builds that define `RAISIN_PART_INT` use integer identifiers
instead.

Examples:

```sh
# Inspect the reference hypergraph.
./exe/raisin hypergraphs/b14.rzn2 stats seed 1

# Build a four-part initial partition.
./exe/raisin hypergraphs/b14.rzn2 part dbfs \
  part_number 4 archfile targets/arch0.arch \
  partfile results/b14-dbfs seed 1

# Refine an existing solution.
./exe/raisin hypergraphs/b14.rzn2 refine dkfm \
  part_file results/b14-dbfs.sol part_number 4 \
  archfile targets/arch0.arch partfile results/b14-dkfm \
  perform 10 tolerance 0 seed 1

# Run the complete multilevel workflow.
./exe/raisin hypergraphs/b14.rzn2 multilevel \
  cluster bsc part ddfs refine dkfmfast part_number 4 \
  archfile targets/arch0.arch partfile results/b14-multilevel \
  perform 5 seed 1

# Evaluate a solution without changing it.
./exe/raisin hypergraphs/b14.rzn2 eval \
  part_number 4 partfile results/b14-multilevel.sol \
  archfile targets/arch0.arch seed 1
```

Create the output directory before using a prefix such as `results/b14-dbfs`.
Without `partfile`, the output is `<graph path>.part.sol`.

## File formats

The repository uses three line-oriented text formats:

- `.rzn2` for directed red-black hypergraphs;
- `.arch` for target architectures;
- `.sol` for vertex-to-part assignments.

Their normative v1.3 grammar and validation rules are documented in
[`docs/formats.md`](docs/formats.md).

## Tests and reproducibility

The test suite validates graph invariants, hand-checkable path and capacity
fixtures, malformed inputs, deterministic replay, every supported algorithm and
the final-v1.2 objective/performance ceilings:

```sh
make check
```

The extended algorithm smoke suite is:

```sh
./run_tests.sh
```

Set `RAISIN_BINARY` to test another executable. The suite writes temporary
solutions outside the repository and removes them on exit.

For scientific comparisons between two builds:

```sh
python3 test/compare_versions.py \
  /path/to/baseline/raisin /path/to/candidate/raisin .
```

Add `--enforce-objectives` to fail when comparable `pmax`, `cut` or `cost`
metrics regress. Always record the graph, architecture, command, seed, compiler
and commit hash with experimental results.

GitHub Actions builds and tests the project with GCC and Clang, then repeats the
suite under ASan and UBSan.

Planned scientific-model, API and scalability work is tracked in
[`ROADMAP.md`](ROADMAP.md).

## Repository layout

```text
include/       public C headers
src/           implementation and CLI
hypergraphs/   reference circuit inputs
targets/       reference target architectures
test/          automated scientific and regression tests
```

The active suites are `test_invariants.c`, `test_arch.c`, `test_capacity.c`,
`test_scientific.c`, `test_cli.py` and `test_algorithms.py`. Historical
interactive programs were removed after their useful coverage was migrated to
automated assertions.

## References

1. Julien Rodriguez, François Galea, François Pellegrini and Lilia Zaourar,
   “A Hypergraph Model and Associated Optimization Strategies for Path
   Length-Driven Netlist Partitioning,” ICCS 2023,
   [doi:10.1007/978-3-031-36024-4_50](https://doi.org/10.1007/978-3-031-36024-4_50).
2. Julien Rodriguez, François Galea, François Pellegrini and Lilia Zaourar,
   “Hypergraph Clustering with Path-Length Awareness,” ICCS 2024,
   [doi:10.1007/978-3-031-63775-9_7](https://doi.org/10.1007/978-3-031-63775-9_7).
3. Julien Rodriguez, *Circuit partitioning for multi-FPGA platforms*, PhD
   thesis, Université de Bordeaux, 2024,
   [HAL tel-04731886](https://theses.hal.science/tel-04731886).

## License

RaiSin is distributed under the GNU General Public License version 3. See
[`LICENSE`](LICENSE) for the complete GPLv3 text. This fork preserves verbatim
the original CEA/Inria copyright, license and warranty notices present in the
inherited source files, without replacing or reinterpreting their wording.
