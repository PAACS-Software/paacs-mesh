# Local qualification — 1.0.0-rc.1

Date: 2026-09-29. This report distinguishes executed checks from future CI
configuration. The original v0.10 archive is the comparison baseline; no
original production source file was deleted. See `docs/CHANGE_INVENTORY.json`.

## Environment

- Linux x86-64, glibc 2.39; AMD EPYC 9V74 host, 8 CPU-equivalent container quota,
  8 GiB memory limit. Timings are shared-host observations with no pinned
  frequency or affinity, not performance guarantees.
- GCC 16.2.0 (conda-forge), C++26, CMake 4.4.3, Ninja 1.13.2.
- Release: `-O3 -DNDEBUG -std=c++26 -fno-fast-math -ffp-contract=off`.
  Production library warnings enabled and treated as errors. IPO and native
  architecture tuning were OFF for qualification and timing.
- Boost 1.92 headers supply independent exact-rational test oracles; the
  installed library has no Boost dependency.

## Executed checks

| Check | Result and scope |
|---|---|
| Original source, GCC 13/C++20 | All 4 originally registered CTest probes passed |
| Modern GCC 16/C++26 build | Complete production library, CLI, test library, unit runner, probes and benchmarks built |
| Modern CTest | 12/12 CTest entries passed, including all 475 inherited unit cases; 274.75 s total |
| Allocation-failure injection | 17 injected failures passed across node/cell append, dedup insertion/index rebuild, archive encoding and decoding |
| Archive hostile inputs | All ten supported element types, independent empty-header fixture, all truncations, 2,000 deterministic byte mutations, invalid counts/flags/indices/coordinates, resource limits, transactional destination checks |
| AddressSanitizer + UndefinedBehaviorSanitizer | 4/4 selected checks passed: hardening, allocation failure, pre-OCCT probe, Delaunay probe; leak detection disabled as explained below |
| ThreadSanitizer | Production `Parallel.cpp` plus the dedicated stress executable passed; 50 rounds with 4,096 outer items, nested work, body exceptions, and concurrent compatibility-setting reads/writes |
| Installed SDK | Independent C and C++ consumer projects configured, linked and passed 2/2 tests |
| CLI end-to-end | Tet4, Tet10, cylinder and Hex8 generation, archive reload and supported audits passed; strict numeric rejection checks passed |
| Determinism/compatibility | Tet4 box identity matched between 1 and 4 requested threads and the original CLI; original CLI read the new Tet10 archive |
| Disk-write failure | `/dev/full` was rejected, including stream-close error reporting |
| GCC 13/C++23 compatibility | Configuration plus targeted archive/allocation and parallel stress builds/runs passed; the entire suite was not rerun under GCC 13 |
| Build configuration gates | Older GCC rejected by GCC-16-required preset; compatibility configuration and IPO support probe succeeded |

The historical full-suite checks include nine semantic baseline hashes, exact
predicate differential oracles, topology transactions, constrained recovery,
layer and multiregion workflows. Its old environment-specific strict hashes
remain guarded and were skipped outside their original reference definition;
semantic hashes still ran. The separate Delaunay goldens and CLI comparisons
ran exact hashes on GCC 16.

## Historical test corrections, with baseline reproduction

Enabling every shipped unit-test translation unit uncovered four assertions
that the original four-probe CTest configuration did not run:

1. The lattice test expected process-global predicate counts to rise. The
   existing compact Delaunay implementation intentionally calls
   instrumentation-free filtered predicates. The test now explicitly requests
   local Delaunay operation counts and verifies those counters. Its 125-point
   insertion, adjacency audit, and hull-volume checks are retained; public
   predicate instrumentation is independently checked by the predicate tests.
2. Two protected-feature tests carried stale exact work-count values. A fresh
   build of **unmodified original production sources**, with the same GCC 16,
   optimization and strict FP flags, reproduced the upgraded values:

| Case | Old expected work | Original source observed | Upgraded source observed |
|---|---:|---:|---:|
| Single-pass interface prefill | 1,430,115 | 1,422,984 | 1,422,984 |
| Prepared-constraint prefill | 206,285 | 199,154 | 199,154 |

Only these measured-work reference constants were refreshed. The checks for
181 nodes, 657 tets, unchanged mesh identity, boundary validity and geometric
quality remain. `evidence/baseline-work-counts.log` intentionally contains the
reproduced old assertion failures. They are diagnostic evidence, not failures
of the final candidate. See the final CTest log for the final outcome.

## Measured performance

| Archive operation | Original | Upgraded | Original / upgraded |
|---|---:|---:|---:|
| Encode | 3.315 ms | 1.873 ms | 1.77x |
| Decode | 7.368 ms | 6.138 ms | 1.20x |

| Delaunay input points | Output tets | Median seconds | Million tets/s | Canonical hash |
|---:|---:|---:|---:|---|
| 50,000 | 335,606 | 0.340603 | 0.985 | `0xe763362f526a11b1` |
| 200,000 | 1,346,539 | 1.365306 | 0.986 | `0xfcd63f974312b432` |

Archive timing is a synthetic serialization/reconstruction workload, not a
meshing-speed comparison. It repeats one Tet4 connectivity 200,000 times;
the output is 8,800,168 bytes. Old and new binaries use GCC 16/C++26 with the
same `-O3 -DNDEBUG` and strict FP options. Each invocation warms once and
reports the median of seven measured iterations. Three interleaved invocations
per implementation are retained; the table uses the median of those medians.
Both implementations produced matching byte counts and archive hashes
(`13931809011004280743`).
Decode timing includes the new preflight validation. Peak RSS was not measured.

Delaunay results describe the retained unconstrained serial kernel: one warmup
and seven measured runs per point count, with built-in lattice/random/near-
coplanar goldens verified before timing. Hashes were stable between measured
runs and insertion/audit checks passed. No speedup over the original Delaunay
algorithm is claimed.

## Scope and remaining qualification

- ASan/UBSan initially encountered a LeakSanitizer runtime failure involving
  this environment's `/proc`/ptrace restrictions. Repeating the selected tests
  with `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1` and
  `UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1` passed. Leak freedom is
  **not** established. The full inherited suite was not sanitizer-qualified.
- TSan evidence covers the dedicated parallel stress workload, not every
  meshing pipeline or independent concurrent sessions. The single external
  dispatcher contract remains in force.
- The Hex8 CLI result passed mixed-element quality checks. The older generic
  topology audit defers non-tetrahedral cells; its `audit_passed` field alone
  must not be interpreted as exhaustive Hex8 topology certification.
- No Windows/macOS/Clang build, sustained fuzz campaign, 50-million-cell run,
  final stable ABI, complete C++26 feature coverage, or remote CI run was
  qualified. IPO was configuration-probed, not benchmarked.
- General CAD import/healing, unrestricted anisotropic meshing, and robust
  handling of every difficult domain remain development work. See
  `docs/RELEASE_READINESS.md` for the acceptance roadmap.

## Evidence and reproduction

`evidence/` retains selected raw logs, compiler settings and benchmark output.
`docs/BUILDING.md` documents tool installation and build/test commands.
`docs/BENCHMARK_PROTOCOL.md` defines fair comparative measurements.

The package script creates a source ZIP, verifies its CRCs and every recorded
SHA-256, and includes `SHA256SUMS.txt`. Compiler/runtime binaries and local build
trees are deliberately excluded. The source candidate requires compilation.
