# Changelog

This project follows semantic versioning for released command-line behavior and
serialized formats. Dates use the ISO `YYYY-MM-DD` form.

## [1.2.0] - Unreleased

### Added

- Portable CMake build alongside the GNU Make build.
- Isolated `build/` and `exe/` build products with generated dependency files.
- Deterministic `seed <n>` CLI option, defaulting to seed `1`.
- Hypergraph invariant tests and deterministic CLI regression tests.
- Cross-version scientific metric comparison harness.
- GCC, Clang, ASan and UBSan GitHub Actions jobs.
- Normative documentation for `.rzn2`, `.arch` and `.sol` files.
- Repository-wide UTF-8 and line-ending policy.
- Canonical GPLv3 license text and upstream fork attribution, while preserving
  the original CEA/Inria source headers verbatim.

### Changed

- C11 compilation is warning-clean with `-Wall -Wextra -Wpedantic -Werror`.
- Scratch memory in hypergraph and vertex-to-hyperedge algorithms is allocated
  per call instead of being shared through static caches.
- Build outputs no longer overwrite source-tree object files.
- CLI help now matches supported modes, algorithms, defaults and limits.
- Default solution names preserve the complete graph path and append
  `.part.sol`.

### Fixed

- Public declarations for hypergraph precomputation and validation helpers.
- Missing string terminators in historical test path allocations.
- Strict numeric CLI parsing, required-option checks and non-zero failure codes.
- Overflow of the default 8-bit `PART` type above 255 partitions.
- Mutation of the precomputed red-vertex mask during criticality computation.
- Unsafe shared scratch buffers that made repeated calls and concurrency fragile.
- Invalid `free()` of pointers returned by `strtok()` in architecture loading.
- BSC merge paths that could use an uninitialized predecessor.
- Non-deterministic BSC results caused by uninitialized feasibility flags.
- Partition readers accepting truncated, malformed or out-of-range assignments.
- Architecture readers dereferencing missing tokens or accepting invalid part IDs.
- Hypergraph readers accepting inconsistent counts, indices and vertex records.
- Readers rejecting a valid final record when the file had no trailing newline.
- Shell tests aborting before they could count and report failures.

### Compatibility notes

- Malformed inputs that were previously accepted or crashed now return failure.
- Explicit help (`help` or `--help`) returns success; missing arguments return
  failure.
- Runs are deterministic by default. Pass a different `seed` to explore another
  pseudo-random sequence.
- The default build supports at most 255 parts. Define `RAISIN_PART_INT` when a
  larger identifier type is required and validate the resulting memory costs.

### Release checklist

- [x] Warning-free GCC build.
- [x] Warning-free Clang build in CI configuration.
- [x] Invariant and CLI regression suites.
- [x] ASan/UBSan CI configuration.
- [x] File-format and reproducibility documentation.
- [x] Confirm GPL version 3 licensing, add the canonical `LICENSE` file and
  preserve the complete original CEA/Inria notices in inherited source files.
- [ ] Run the hosted CI jobs on the final release commit.
- [ ] Record a final v1.1-to-v1.2 comparison report on representative circuits.

## [1.1.0]

Baseline release represented by the `v1.1` Git tag.
