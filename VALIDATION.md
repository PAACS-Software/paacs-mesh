# PAACS Mesh 1.1.0 — integration validation

The [1.2.0-rc.1 candidate validation](docs/RELEASE_READINESS.md) records the
latest results. This page describes the earlier 1.1.0 integration build.

Verified 2026-10-02 UTC on Linux x86_64 with GCC 16.2.0, C++26, CMake 4.4.3,
Ninja 1.13.2, strict floating point, portable CPU flags and warnings as errors.
OCCT is disabled. This release adds SDK/application integration; all existing
meshing, sizing, predicate, core-pool and audit implementation files are unchanged.

## Fresh checks

| Check | Result |
|---|---|
| Release CTest selection | **8/8 passed**, 8.97 s: four SDK/UI suites plus CLI version, hardening, control fields and spatial controls |
| Installed package consumers | **13/13 passed**, 1.05 s, including the new C UI and native background-job examples |
| ASan/UBSan with leak detection | **3/3 SDK suites passed**, 1.41 s; core, adapter and tests built with address/undefined sanitizers; `detect_leaks=1` |
| SDK ThreadSanitizer check | Passed 637 checks; new JobService and its tests instrumented, linked against the ordinary Release core. This is not full-core TSan qualification. |
| Native job semantics | Five byte-identical archive comparisons, copied profile/position lifetimes, immutable ownership, budgets, exceptions, cancellation and active shutdown |
| Cross-service dispatch | 48 jobs, no overlapping engine dispatches; concurrent snapshot polling |
| Publication races | 100 accepted-cancel/terminal-publication races passed |
| C adapter parity | Five direct-native mesh comparisons; 261,532 scalar/page/connectivity/classification checks, including result lifetime after service destruction |
| Python protocol | Ten test methods, including 18 primitive shape/element combinations, fields/profiles, explicit profile disable, malformed/oversized input recovery and cancellation |
| Local HTTP example | Live submission, polling, summary, full 1,296-face cylinder surface preview, entity tags, default hash-off and host/origin/token/export checks passed |
| Version and installation | Source/CMake/SDK/CLI report 1.1.0; installed static library SHA256 matches the tested build |

The browser renderer received a JavaScript syntax check. Actual browser rendering,
mouse interactions and recipe-import event flows were **not automatically exercised**:
no browser binary was available and the test-browser download endpoints were
unavailable. The local HTTP workflow and preview data were exercised separately.
Windows/macOS and frontend deployment are not newly qualified.

The local development workspace retains the full logs, JUnit records, build
metadata and binary hashes under `evidence/1.1.0/`. The initial sandboxed leak check could
not inspect its own threads; the final successful run used the same instrumented
binaries with process-inspection access. Two zero-byte generated objects and a
stale generated shared-library link were rebuilt before the final run. No source
checks or sanitizer findings were suppressed to obtain these results.

## Interactive overhead measurement

`benchmarks/ui_jobs.cpp` compares a direct primitive call with a queued call,
including its summary and result publication. Fifteen measured repetitions per
configuration follow two warmups. Both paths produce the same mesh hash; hash
verification occurs outside the timed regions when automatic hashing is off.

| Target | Nodes / total cells | Direct median, ms | Job median, ms | Automatic result hash |
|---|---:|---:|---:|---|
| 0.4 | 1,651 / 9,072 | 1.607 | 1.486 | Off (default) |
| 0.2 | 1,651 / 9,072 | 1.663 | 1.493 | Off (default) |
| 0.1 | 6,951 / 39,600 | 8.842 | 7.515 | Off (default) |
| 0.4 | 1,651 / 9,072 | 1.470 | 2.494 | On |
| 0.2 | 1,651 / 9,072 | 1.405 | 2.422 | On |
| 0.1 | 6,951 / 39,600 | 8.011 | 11.005 | On |

These small local samples include scheduling/cache effects. The lower job medians
with hashing off are **not an algorithm speedup**. The two coarser targets produce
the same curvature-limited mesh. The useful design result is that strict hashing
can be omitted from the interactive path while preserving every native mesh
check. These measurements do not establish scaling on large industrial CAD or
comparative speed against another mesher. Logical result payloads were 0.525 MiB
and 2.300 MiB, respectively; these are not process RSS measurements.

## Retained qualifications and scope

Original evidence is preserved unchanged. See the [1.0.0 verification record](docs/VALIDATION_1.0.0.md)
and [rc.11 implementation qualification](docs/VALIDATION_RC11.md) for the inherited
38/38 full suites, 475 unit tests, previous sanitizer/allocation campaigns,
geometry corpus and performance work. Those full campaigns were **not rerun in
full for 1.1.0**. Fresh SDK tests and selected existing regressions are listed above.

The new native layer supports five adapter routes and custom native tasks. The
portable C/JSON bridge currently covers primitives and adaptive boxes, with 68
scalar options and field/profile arrays; it is not full transport coverage of
every native operation. Models, geometry controls and layers remain accessible
through C++. [UI integration](docs/UI_INTEGRATION.md) documents the boundaries.

General curved CAD, broader OCCT integration, protected-boundary refinement,
curved CFD layer/core coupling, GPU acceleration and matched external comparisons
remain future work. This release makes no new claim of universal CAD robustness,
downstream solver convergence.

## Reproduce

```sh
cmake --preset gcc16-release -DPAACS_BUILD_UI_BRIDGE=ON
cmake --build --preset gcc16-release --parallel 2
ctest --test-dir build/gcc16-release --output-on-failure \
  -R 'paacs.sdk|paacs.hardening|paacs.sizing.control_fields|paacs.sizing.spatial_controls|paacs.cli.version'
python3 tests/ui_http.py --library build/gcc16-release/libpaacs_ui.so
cmake --build --preset gcc16-release --target paacs-bench-ui-jobs
build/gcc16-release/paacs-bench-ui-jobs
python3 scripts/verify_version.py --final
python3 scripts/package_source.py
```

Use the documented sanitizer presets for an instrumented build. A complete TSan
preset instruments the core as well; the narrower current SDK run must not be
confused with that larger qualification. The HTTP check needs loopback socket
access. The source ZIP excludes build outputs and includes its per-file SHA256
manifest; the packager checks CRCs and file hashes before writing a ZIP checksum.
