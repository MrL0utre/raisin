# RaiSin roadmap after v1.3

Version 1.3 completes the scientific-validation and target-constraint work
planned after the v1.2 consolidation baseline. The following API and scalability
work should be developed on top of this release.

## v1.3 — Scientific validation and target constraints (completed 2026-08-02)

1. [x] **Make target capacities effective.** Link demand is routed, measured,
   repaired in partition-producing modes and rejected when infeasible.
2. [x] **Represent per-part resource capacities.** Every vertex-weight dimension
   can be matched to an explicit capacity per FPGA, while legacy balance remains
   available when the matrix is omitted.
3. [x] **Validate topology semantics.** Missing links are distinct from zero-delay
   links, shortest routes are deterministic, and disconnected or duplicate
   topologies are rejected.
4. [x] **Build hand-checkable fixtures and golden algorithm cases.** Small path,
   topology, resource and infeasibility cases have exact expected values; HEM,
   BSC, DBFS, DDFS, CCP, KFM, DKFM and DKFMFAST have deterministic baselines.
5. [x] **Expand scientific regression coverage.** Structured JSON bounds `pmax`,
   cut, cost, balance, runtime and peak memory across the reference circuit and
   the small non-complete-topology fixtures.
6. [x] **Modernize historical tests.** Useful coverage from the five obsolete
   exploratory C programs was migrated to maintained assertions, then the stale
   drivers were removed from the automated-test namespace.

Exit criterion: capacity-feasible solutions and objective metrics are verified
on multiple circuit/topology pairs and compared against a recorded v1.2 baseline.

Evidence is maintained in `test/test_scientific.c`, `test/test_capacity.c`,
`test/test_algorithms.py` and `test/baselines/v1.2-b14-arch0.json`.

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

The fork must retain its CEA-LIST upstream attribution, the complete original
CEA/Inria source notices and GPL version 3 licensing. Every later release should
also publish:

- a completed changelog and migration notes;
- hosted CI and sanitizer results;
- a benchmark comparison against the previous release;
- hashes and provenance for published experimental inputs.
