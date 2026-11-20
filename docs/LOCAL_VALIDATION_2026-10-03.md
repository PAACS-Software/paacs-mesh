# Local Windows validation and condensation review

This record covers the modified Desktop checkout and the supplied OCCT/ABC
inputs. It began on October 2 and continued October 3, 2026. The older GitHub
checkout was inspected read-only; it was not rebuilt or newly qualified.
Raw records are in `evidence/local-2026-10-02/`.

## What condensation changed

The new package is **not equivalent to the previous CAD implementation**.
In the initial package snapshot, source plus headers shrink only about 4%
(4,318,590 to 4,146,784 bytes), while their file count increases from 167 to 198.
Those counts precede this session's repairs. Exact predicates, core mesh audits
and much transaction/storage infrastructure are preserved. Primitive sizing,
geometry controls and frontend job support have substantial new coverage.

The removed OCCT-backed adapter represented curved surfaces, oriented uses,
XDE assembly occurrences, units and provider bindings. Its replacement accepts
one planar, straight-edge solid and bakes it into a polyhedral model. It
explicitly refuses curved STEP faces and assemblies. Persistent V1 selectors,
model/provenance archives and use-level projection contracts were also removed.
The bounded non-star one-Steiner recovery extension and focused regression
witnesses are absent. Historical old successes do not establish that all of
those capabilities were production-ready, but their implementation and focused
coverage were lost.

See the detailed comparison (`evidence/local-2026-10-02/BASELINE_COMPARISON.md`)
and its JSON inventories/hashes. The old owner amendment retires standalone
T2-C/T5 and legacy per-cell topology routes; those should not be copied back
merely because they exist in the old checkout.

## Baseline and toolchain

| Item | Recorded local value |
|---|---|
| OS | Windows 11 Home x64, 10.0.26200 |
| CPU | Ryzen 9 8945HS, 8 cores / 16 logical processors |
| RAM | 33,601,921,024 bytes, approximately 32 GiB |
| Compiler | MSVC 19.44.35207.1, VS Build Tools 17.14 |
| Build tools | CMake 4.0.2; MSBuild 17.14.8; SDK 10.0.26100.0 |
| Boost | 1.89 headers from the existing vcpkg installation |
| Executed language mode | Explicit C++23 compatibility, `/std:c++latest` |

C++26 configure fails because this CMake/compiler combination does not know
`cxx_std_26`. The project default remains C++26; no C++26 qualification is
claimed from these Windows runs.

The inherited unit harness skips its reference-environment hash battery and
`__int128` differential predicate oracles on this compiler. Other exact,
deterministic and geometry checks still run. This local work adds no new
ASan/UBSan/TSan qualification, and does not rerun the reference Linux toolchain.

The original downloaded `paacs-mesh-1.1.0.zip` was extracted into an isolated
build input. All 1,117 manifest-listed source files matched before and after
testing. Its unmodified build passed. Initial CTest ran 42 tests: 34 passed,
eight failed. Two UI launch timeouts were corrected by adding the built shared
library to `PATH`, and both reruns passed. The six remaining outcomes were:

- inaccurate adaptive-box reported volume;
- mapped-Hex test hashing uninitialized structure padding;
- three allocation-failure executables terminating with `0xc0000409`;
- the inherited unit suite exceeding its original 600-second timeout.

The unit process remained CPU-active and reached 388 case starts. Its last
active case passed independently in 17.773 seconds. The final local invocation
uses an explicit 1,800-second unit timeout; the original timeout result remains
in `unmodified-tests.xml` and `unmodified-summary.json`.
The longer native run completed in 750.30 seconds and revealed two exact
work-count assertions whose totals include standard-library comparator calls.
The untouched executable reproduces both failures: `t2c.fast_single_pass`
reports 1,418,116 against the GNU reference 1,422,984, and `t2c.v4_prepared`
reports 194,286 against 199,154. Both differences are 4,868 work units; all
geometry, strict/semantic hash and quality assertions in those cases pass.
These are inherited Windows baseline failures, not changes introduced by the
local implementation repairs. The reference GNU totals remain recorded and
checked under the qualified reference compiler/library/platform guard.

A scoped diagnostic linked to the untouched test library measures 26,425
comparator calls in each case: coordinate sorting 2,538, insertion sorting
2,551, selected-cell sorting 6,871, global `std::set` ordering 1,757, and partial
edge/face sorting 12,708. These call counts depend on the STL implementation;
they are not floating-point volume error. The downstream construction,
recovery and coordinator work is shared by both routes. Their total work after
subtracting only these five explicitly charged comparator counters is
1,391,691 and 167,861 respectively. The cache route accounts for its additional
work explicitly; no production work counter or budget was changed.

The repaired checks enforce exact accounting decomposition and exact
non-comparator component work on every platform, preserving the geometry
goldens. Exact total-work checks remain for Linux x86_64 GCC 16.2/libstdc++ 16
without debug/parallel STL instrumentation, and for this exact MSVC compiler
build, STL update, iterator-debug mode and strict-FP release configuration.
Other environments report their comparator totals and still enforce the exact
non-comparator checks; no numeric tolerance replaces either assertion. The
new invariants pass against the untouched native library. Their reference
Linux run and the GNU comparator breakdown were not newly executed here, so
the full 4,868-unit cross-library breakdown remains an inference to be checked
by that reference environment. Receipts are `work-counter-original-fast.log`,
`work-counter-original-prepared.log`, `work-counter-components-original.log`
and `work-counter-portable-original.log` in the local evidence directory.

## Repairs and additional checks

- Split installed OCCT targets from the core export. Core-only consumers no
  longer discover OCCT; optional/required component behavior is checked in
  isolated consumer projects. Restored `PAACS::Mesh` and `PAACS::MeshOCCT`
  target aliases, without claiming restoration of the removed APIs/ABI.
- Fixed adaptive volume accumulation with a compensated sum. Exact integer
  determinant accounting proves the original meshes already had volume one;
  the reported value accumulated enough error to fail on Windows, where
  `long double` has double precision. All ten measured reports now equal one;
  the original test tolerance was retained.
- Made the two fixture work-count checks portable while preserving exact
  qualified-toolchain totals and exact work-accounting/non-comparator
  invariants. Their untouched Windows reproductions and reference GNU
  goldens remain separate; this does not claim an old/new speed improvement.
- Moved empty replacement-mesh construction out of `noexcept` publication.
  MSVC's empty hash table may allocate; previously an injected failure at
  geometry allocation ordinal 431 terminated the process. Preparation can
  now return the existing allocation-failure diagnostic before publication.
  The single-use construction storage preconstructs its replacement mesh;
  publication only swaps existing objects. The extended campaign then reached
  Volume-profile HexDominant construction and terminated at allocation ordinal
  1 (trace variant 16): the private `primitive_detail::profiled` helper
  constructed staging meshes inside an unguarded `noexcept` function.
  Removing that private specification lets its sole caller's existing
  allocation/general catch return the error and preserve the destination.
  The public `meshPrimitive` API remains `noexcept`; all potentially throwing
  staging happens before its final destination exchange. Both concrete
  before-fix traces are retained as `allocation-before.stderr.txt` and
  `allocation-after-storage.stderr.txt`, and the existing exhaustive allocation
  campaign covers both paths. Final rerun status is recorded below.
- Replaced a probe's raw structure-byte hash with the existing fieldwise mesh
  hash. Coordinates, cell fields and connectivity already matched; only
  padding differed across repeated runs.
- Removed a false assumption that solid count equals connected mesh component
  count. A shared-face, two-solid witness measures one component, while two
  separated solids publish with two. The wrapper still refuses shared-interface
  publication through its existing exterior-only final boundary audit; that
  separate limitation was retained.
- Restored default face-projection trim checks, finite-value checks and
  unchanged outputs on refusal. Allocation failures propagate distinctly.
- Made archive publication use an exclusively created temporary file beside
  the destination, validate the write and close, then replace the destination
  by rename. Invalid input, allocation failure and write/replace failure
  preserve an existing destination. Tests check successful replacement and
  SDK readback, a missing parent, a directory destination, Windows locked-file
  replacement, and the actual writer's allocation ordinals. Cleanup of owned
  temporary files is attempted on failure. Atomic publication does not promise
  crash-durable storage or preservation of filesystem metadata.
- Explicitly transfer STEP coordinates to millimetres and record source units.
  A 2 x 3 x 4 box round-trips STEP declared in millimetres and metres with
  volume 24 after SDK normalization. Adapter tests also cover scaled dimensions,
  reversed face orientation, shared identities, projection and classification.
- Added opt-in retained allocation/rollback campaigns and registered previously
  unregistered CAD CLI checks. Added structured CLI JSON, bounded corpus
  sampling, repeated processes, timeout/crash detection and Windows peak
  working-set measurement.
- Kept OCCT parser allocations outside the claimed file/topology budgets:
  input bytes, logical mesh output, and audit storage are bounded separately;
  no process-memory limit is enforced. Import has no caller cancellation token.
  The corpus runner kills its direct child at the per-file timeout and retains
  raw logs and publication receipts. Its children inherit Windows error-dialog
  suppression, so DLL/entry-point failures return native crash status instead
  of waiting invisibly on a modal dialog. Eight native malformed-input,
  argument, file/output/audit-budget and quality-refusal probes also preserved
  the bytes of an existing output file.
- Added exact node-plus-cell construction floors for mapped/swept work limits.
  Baseline constrained Tet4 construction has no full work counter and remains
  outside that cap. Logical output and audit budgets are not a process-memory
  cap. The runner measures process memory and enforces a child-process timeout.
- Prevented local executables, objects, DLLs and SDK/runtime artifacts from
  entering source ZIPs. Local evidence is excluded by default; an explicit
  option includes text logs and reproducible witnesses. Both packaging modes
  passed full archived-file checksums, CRC and artifact-exclusion checks.

Global MSBuild vcpkg integration initially copied OCCT 7.9.3 DLLs beside
executables linked against the supplied 8.0.1 SDK. Adjacent DLLs override
`PATH`. We removed 58 build-output copies only after verifying their hashes
against vcpkg, and disabled that deployment in this build. The SDK itself
contains prebuilt Windows x64 libraries and headers, not complete kernel source.
Its combined archive already includes the matching third-party dependencies.
The extraction helper retains approximately 359 MiB of required SDK files in
`build/deps/`, excluding symbols and unrelated payloads.
The supplied SDK version also agrees with the official
[OCCT 8.0.1 release](https://github.com/Open-Cascade-SAS/OCCT/releases/tag/V8.0.1).

## Reproduce the modified Windows build

Run from the repository root in PowerShell. Set `$occtArchive` to the path of
the combined SDK ZIP. The SDK helper is idempotent for the selected destination.
The combined SDK includes the third-party runtime.

```powershell
. ./scripts/setup_occt_sdk.ps1 `
  -Archive $occtArchive `
  -ActivateRuntime
cmake -S . -B build/windows-msvc-release -G 'Visual Studio 17 2022' -A x64 `
  -DPAACS_CXX_STANDARD=23 -DPAACS_BUILD_UI_BRIDGE=ON `
  -DPAACS_BUILD_BENCHMARKS=ON -DPAACS_BUILD_EXTENDED_PROBES=ON `
  -DPAACS_BOOST_INCLUDE_DIR=C:/vcpkg/installed/x64-windows/include `
  -DPAACS_WITH_OCCT=ON `
  -DOpenCASCADE_DIR="$PWD/build/deps/occt-8.0.1/opencascade-8.0.1-vc14-64/cmake" `
  '-DCMAKE_CXX_FLAGS=/MP2 /EHsc' `
  '-DCMAKE_VS_GLOBALS=VcpkgEnabled=false;VcpkgApplocalDeps=false' `
  -DPAACS_UNIT_TEST_TIMEOUT=1800
cmake --build build/windows-msvc-release --config Release --parallel 2
$env:PATH = "$PWD/build/windows-msvc-release/Release;" + $env:PATH
ctest --test-dir build/windows-msvc-release -C Release --parallel 2 --output-on-failure
cmake --install build/windows-msvc-release --config Release --prefix out/windows-install
```

Build/install/component logs preserve the executed commands and results.
The original package's manifest is retained separately; locally changed source
must not be described as byte-identical to the original release.
`evidence/local-2026-10-02/final-repaired-run.ps1` preserves the final invocation,
including the inherited Windows dialog policy for unattended subprocesses.

## Progress against the eight stages

| Stage | Current qualification boundary |
|---|---|
| 1. Baseline | Unmodified Windows build and tests recorded; C++23 compatibility only. Modified full-run result is recorded below. |
| 2. OCCT | Supplied 8.0.1 x64 SDK inspected and linked; real STEP write/read units, dimensions and volume checked. |
| 3. Adapter | Optional planar adapter strengthened. Native curved BRep, tolerance/use-level topology, seams, assemblies and persistent body selection remain restoration work. |
| 4. CAD meshes | Supported planar fixtures and a real ABC sample exercised. Refusals remain explicit; successful planar tests do not qualify general CAD meshing. |
| 5. Sizing | Existing primitive/polyhedral controls tested. Imported curved-cylinder face refinement remains blocked by adapter scope; all three requested imported-CAD demonstrations are not yet qualified. |
| 6. Hardening | Native fault injection, malformed inputs, publication and resource refusals tested. Full process memory limiting and complete Tet4 work accounting remain incomplete. |
| 7. Performance | Corpus runner records per-stage timings, counts, quality and peak working set. Surface and volume generation still share one timer; no old/new scaling or competing-mesher speed claim. |
| 8. Studio | Existing job/bridge consumers exercised. Separate pinned studio repository and complete import/select/section/save workflow remain unimplemented. |

Restore the optional native CAD geometry provider first, using focused old
contracts and witnesses against the current authoritative core. Then qualify
curved STEP and body identity before CAD face sizing and the studio workflow.

## Final run results

The final repaired build and local installation passed. All five isolated
installed component checks passed their expected outcomes, and all 14 installed
consumer examples passed again after the final library rebuild (2.59 seconds),
including actual through-hole CAD import and meshing, C, C++, UI and background
jobs. Installed library hashes match the final build. Receipts are
`package-occt-components.json`, `final-repaired-install.log` and
`installed-consumer-evidence.json`; the earlier consumer run remains separate.

The deterministic ABC screening selected 51 of 10,000 downloaded files:
first 30, six smallest, 12 spread positions and four deliberately added planar
candidates, with overlap removed. Each ran in three fresh processes, with a
60-second timeout. This is a focused support sample, not a random success-rate
estimate. All 153 runs completed:

| Outcome | Files | Runs |
|---|---:|---:|
| Audited mesh and archive | 5 | 15 |
| Curved surfaces unsupported | 29 | 87 |
| Multiple solids unsupported | 16 | 48 |
| Planar surface-route failure | 1 | 3 |
| Crash, timeout or failed-result publication | 0 | 0 |

| ABC file ID | Nodes | All cells, including boundary | Minimum SubTetMeanRatio |
|---|---:|---:|---:|
| 22 | 125 | 160 | 0.756 |
| 638 | 46 | 90 | 0.275 |
| 850 | 32 | 37 | 0.213 |
| 5956 | 64 | 83 | 0.305 |
| 6807 | 70 | 92 | 0.280 |

Every successful archive was independently reloaded through the SDK,
re-audited, re-encoded byte-for-byte and compared with its pinned hash. Hashes
were stable across all repeats. Maximum relative CAD/mesh volume error was
1.21e-15. Maximum classified boundary-node geometric residual was 2.00e-11 mm,
or 3.67e-15 of the bounding-box diagonal. This node/trim check is not a
whole-face Hausdorff certificate.

Successful-run medians were 12.383 ms import, 0.0867 ms preparation, 0.0394 ms
combined generation, 0.0392 ms boundary audit, 0.021 ms quality, 1.321 ms
archive publication, and 14.403 ms CLI total. Fresh-process wall time was
130.878 ms. Peak working set was 24,731,648 bytes for successes and 131,706,880
bytes across the sample. These are small supported planar meshes, not a
throughput/scaling benchmark. The earliest unsupported timings may overlap
brief install/contract checks; successful cases ran later without build/test
overlap.

ABC 11 consistently spent about 27–30 seconds importing before curved-face
refusal. A separate diagnostic profile attributed approximately 94.8% of its
measured import phases to OCCT `TransferRoots`, rather than mesh generation.
That one diagnostic profile overlapped later validation activity and is not an
additional clean timing benchmark. Volume/baking phases were never reached.

ABC 5454 imports as a valid planar solid (67 faces, 193 edges, 126 vertices),
but its current Tet4 surface candidate stops on face tag 34 because Steiner
insertion is disabled. The separate planar constraint-recovery surface option
can triangulate all faces; this does not qualify its subsequent volume route.
The failure preserves the destination and creates no archive. Native witnesses
and detailed diagnostics retain this case for future restoration work. A
private experiment feeding the planar constraint-recovery candidate into the
current volume route still stops at `NonStarEdgeCavity` / `NonCycleLink`, so
switching the surface algorithm alone cannot qualify this body. Production
routing was retained. Failed generation time and nested face/recovery details
are now retained instead of reporting zero time and a coarse volume failure.

The first full native run completed 54 tests: 52 passed, with the two unit
work-count assertions and the second allocation termination described above.
After the final repairs, **all 54 CTest entries pass**, with **475 unit cases**
completing in 719.35 seconds and total CTest wall time 719.41 seconds. The
general, complex and spatial allocation tests and all four opt-in topology
campaigns pass. The profiled-sizing test passes 75 meshes and 58 expected
refusals. Reference-toolchain hash/oracle skips described above still apply.
Logs are `final-repair-build.log` and `final-repaired-tests.log/.xml`.
The interrupted default sandbox run exposed Python/shared-library launch
restrictions; those invocation failures are kept separate from the native runs.

The actual final CAD executable also reran both ABC5454 failure-publication
checks. It records positive failed-generation time and the named nested face
recovery diagnostic, preserves an existing archive byte-for-byte, and creates
no absent/partial archive. `abc-planar-cli-check.json` records that executable;
the earlier patched-object witness receipt is retained separately. The 153-run
screening retains its original input/executable hashes rather than being
relabeled as a new timing run.

Raw corpus input hashes, commands, outcomes and repeat records are in
`abc-screening/`. Independent reload checks, construction-floor refusals and
the import profile are in `cad-native-evidence/`. The concise CAD findings are
also in `evidence/local-cad-inspection/native-results.md`.

## Reviewed source snapshot

`out/paacs-mesh-1.1.0-local-reviewed.zip` is a locally modified source snapshot,
not an official re-release. It excludes build outputs, binaries, downloaded
SDKs and local evidence. The local report references those evidence records
in this workspace; they are deliberately absent from the compact archive.
Its embedded `SHA256SUMS.txt` describes the reviewed contents. The original
root manifest is preserved unchanged for comparison with the downloaded release.

```powershell
python scripts/package_source.py --preserve-source-manifest `
  --output out/paacs-mesh-1.1.0-local-reviewed.zip
python evidence/local-2026-10-02/final-validation-receipt.py
```

Packaging checks every archived file hash and ZIP CRC. The receipt additionally
compares the archive with the current source and original manifest, records
changed/added/excluded paths, verifies the 54/14 test receipts, and writes
`out/local-reviewed-SHA256SUMS.txt` plus
`evidence/local-2026-10-02/final-validation-receipt.json`.
