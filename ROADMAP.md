# RaiSin roadmap after v1.2

Version 1.2 is a consolidation release: it makes the current implementation
buildable, deterministic, testable and defensive without changing the research
model. The following work should be developed on top of that baseline.

## v1.3 — Scientific validation and target constraints

1. **Make target capacities effective.** The architecture loader stores a
   connection-capacity matrix, but the optimization code currently consumes only
   the delay matrix. Define the capacity semantics from the thesis, enforce them
   in partitioning/refinement, and report violations as first-class metrics.
2. **Represent per-part resource capacities.** Match every vertex weight
   dimension with an explicit capacity for every FPGA part. Replace the current
   balance-only proxy with feasibility checks, while retaining balance as an
   optional secondary objective.
3. **Validate topology semantics.** Distinguish a missing connection from a
   zero-delay connection, reject or explicitly support disconnected targets, and
   add non-complete topology fixtures.
4. **Build hand-checkable algorithm fixtures.** Add small golden cases for HEM,
   BSC, DBFS, DDFS, CCP, KFM, DKFM and DKFMFAST, including infeasible cases and
   expected objective values.
5. **Expand scientific regression coverage.** Add representative circuits and
   target families, store structured JSON results, and define accepted ranges for
   `pmax`, cut, cost, balance, runtime and peak memory.
6. **Modernize historical tests.** Convert useful assertions from the five
   exploratory `test_*.c` programs into maintained unit tests, then archive the
   remaining experiment drivers outside the automated-test namespace.

Exit criterion: capacity-feasible solutions and objective metrics are verified
on multiple circuit/topology pairs and compared against a recorded v1.2 baseline.

## v1.4 — Stable library and versioned data model

1. Separate the CLI from a documented `libraisin` API with explicit ownership,
   `const` contracts and structured error values instead of process-wide
   `exit()` macros.
2. Introduce version markers for graph, architecture and solution formats, plus
   backward-compatible readers and round-trip tests.
3. Replace the fixed 16 MiB line buffers with bounded streaming parsers that
   report file, line and field context.
4. Replace global `rand()` state with an RNG object passed through algorithm
   contexts, preserving seeded reproducibility and enabling concurrent runs.
5. Choose and document partition identifier widths (`uint8_t`, `uint16_t` or
   `uint32_t`) from problem size, with compatibility tests above 255 parts.
6. Add API documentation and package metadata suitable for downstream builds.

Exit criterion: applications can link RaiSin as a library, handle errors without
process termination and read explicitly versioned files.

## v1.5 — Scalability and complete thesis model

1. Profile representative workloads before optimizing; publish time and memory
   complexity measurements per phase.
2. Optimize incidence traversal and allocation hotspots with reusable,
   caller-owned workspaces whose lifetime and thread-safety are explicit.
3. Add safe parallelism only after the RNG and algorithm contexts are isolated.
4. Represent sets of directed acyclic hypergraphs sharing red vertices, rather
   than requiring a pre-flattened single hypergraph.
5. Add checkpointing and machine-readable result manifests containing input
   hashes, seed, compiler, build options, commit and metrics.

Exit criterion: the complete multi-DAH model runs reproducibly on larger
benchmarks with measured scaling and resumable experiments.

## Release governance

Before tagging v1.2, the copyright holders must resolve the conflicting license
notices described in `CHANGELOG.md`. Every later release should also publish:

- a completed changelog and migration notes;
- hosted CI and sanitizer results;
- a benchmark comparison against the previous release;
- hashes and provenance for published experimental inputs.
