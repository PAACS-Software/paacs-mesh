# 1.2.0-rc.1 validation

## Core

The Windows and Linux core CLIs generated and audited the same 43-Tet4 box
with mesh hash `16723556976053534829`. The Windows run used only system
runtime paths. The Linux executable linked only the standard compiler and system
runtimes.

The standalone macOS core job also generated and audited the 43-Tet4 box
with that hash on Apple Silicon.

The Linux reference CI build passed all 53 CTest suites in 740 seconds,
including the combined unit suite, and all 12 installed-consumer tests.
It used GCC 16.2.0, C++26, CMake 4.4.3, warnings as errors and the core-only
configuration with the UI bridge disabled. The installed license, notice
and citation files matched the source files byte for byte.

The earlier local Linux C++26 build used GCC 16.2.0 and CMake 4.4.3 with
warnings as errors.
Fifty-five of 56 CTest suites passed. The combined unit suite reached its
600-second time limit. All 25 focused reruns covering its remaining cases
passed, including the provider fault-injection case that was running at the
deadline. That case alone took about 309 seconds. The Linux release CI job
gives the combined suite 1,200 seconds.
All 13 downstream tests passed against the installed Linux SDK.

The full ThreadSanitizer run passed 54 of 55 suites and found a lifetime
defect in the SDK job test's blocking helper. The helper now owns the promise
and future used by its task and signals readiness after publishing progress.
A regression case destroys the helper while its task is active. The repaired
SDK suite passed ten consecutive repetitions under ThreadSanitizer and ten
under AddressSanitizer/UndefinedBehaviorSanitizer, including leak detection.

The full AddressSanitizer/UndefinedBehaviorSanitizer run passed 50 of 55
suites, including the combined unit suite. Four recovery tests mixed their
replacement allocation hooks with the sanitizer's default nothrow allocator.
They now use matching test-only hooks. A small regression reproduced the
original mismatch and checked both scalar-only and scalar/array replacement
patterns with GCC 13 and GCC 16 under both sanitizers and leak detection.
The exhaustive refill audit reached its 180-second test deadline. Sanitized
builds now allow 900 seconds for that audit; ordinary builds retain 180 seconds.
The mesher's work and resource limits are unchanged.
Both edge-recovery suites passed two consecutive ASan/UBSan executions after
the hook fix. The cavity and region recovery suites also passed. The exhaustive
refill audit completed in 466 seconds with its full orientation, boundary
incidence and tetrahedron-intersection checks retained.
All five affected suites passed the Windows release build, including the full
Steiner allocation-failure campaign and the exhaustive refill audit.

## Native STEP

The Windows build used MSVC 19.44, CMake 4.0.2, C++23 compatibility mode and
OCCT 8.0.1. The [51-file run](step-corpus-1.2.0-rc.1.csv) reused the same
input hashes and settings as the previous matched run. Each process used
advanced Tet4 recovery, a size factor of 0.15, one meshing worker, one CAD
validation worker, a 90-second deadline and a 2 GiB sampled working-set
watchdog. Up to three processes ran concurrently.

| Outcome | Files |
|---|---:|
| Published and independently audited | 32 |
| Deadline | 14 |
| Unresolved body interface | 2 |
| Working-set stop | 1 |
| Invalid CAD shape | 1 |
| Meshing failure | 1 |
| Total | 51 |

All 51 outcomes matched the previous run. The 32 published archives were
byte-identical to their previous archives. A separate screen verified positive
Tet4 orientation, opposite shared-face winding, boundary incidence, body tags
and agreement with CAD volume.

After the GCC warning cleanup, the rebuilt Windows STEP executable differed
from the corpus executable only in its PE header and debug timestamps.

The Windows CLI check covers invalid-input and archive preservation cases.
The committed cylinder fixture exercised a saved face and box sizing recipe.
The export check recovered all three face groups and body group 1 from the
C3D4 include and matched their oriented connectivity against the native archive.

Five focused Windows CTest suites passed. They cover the native CLI,
boundary, volume validation, meshing adapter and version command.
All 15 downstream tests passed against the installed SDK with its DLL and
the OCCT runtime on `PATH`.

The release-preparation GitHub Actions jobs ran the
[automatic STEP build](BUILDING.md#automatic-step-build) on Windows Server
2022, Ubuntu 24.04 and macOS 15 on Apple Silicon. Each compiled pinned OCCT
8.0.1 from source and passed the cylinder STEP smoke test and mesh audit using
its locally built shared libraries. Windows and macOS also passed all 49 CLI
checks. Windows passed the saved sizing recipe and C3D4 face and body tag check
against its own OCCT DLLs. The macOS build used Xcode 16.4.

An earlier local Ubuntu 24.04 build with GCC 13.3 also passed the 49 CLI checks.
The release-preparation Linux STEP job ran the cylinder check rather than the
full CLI script.

In ABC corpus case CAD06, the box reached a median Tet4 edge size of 0.294 mm
against a 0.25 mm clamped field request. Its 95th percentile was 0.608 mm.
The [plate-with-hole example](../examples/elasticity_plate_hole/README.md)
illustrates combined sizing on the project's generated STEP geometry. The
[asset provenance](ASSET_PROVENANCE.md) identifies the source for CAD06.

## Source and SDK preparation

A clean Windows build from the source ZIP passed all 15 installed C/C++
consumer tests and all 49 native STEP CLI checks. The five installed-package
component checks covered a core-only consumer, optional OCCT availability
and rejection of a required unavailable OCCT dependency.

The installed SDK includes `LICENSE`, `NOTICE`, `LICENSING_STATUS.md` and
`CITATION.cff` under `share/doc/PAACSMesh`. Their bytes matched the files in
the tested source archive. The citation passed CFF 1.2.0 schema validation
and identifies Paul LeDewhurst with his ORCID.
A parent-project configuration check confirmed that the PAACS metadata uses
its own documentation directory when the library is embedded.

A saved cylinder sizing recipe combined face sizing, a refinement box,
minimum size and growth. Reopening the recipe and exporting C3D4 preserved
face groups 1, 2 and 3 and body group 1, including their oriented connectivity.

The source ZIP passed CRC and per-file checksum checks. A credential scan
found no findings in its contents or the private preparation history.
Relative documentation file links were checked against the package contents.
Packaging equivalent LF and CRLF source trees produced byte-identical ZIPs,
with all 32 binary figures unchanged. The file manifest was also checked
against Git's stored source files.

## Earlier qualification

[Core and SDK validation](../VALIDATION.md) records the checks from earlier
releases. The native STEP recovery experiments are documented in the
[CAD development notes](NATIVE_CAD_HISTORY.md).
