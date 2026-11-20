# V0.10 package verification

## Provenance

- Authoritative base: attached PAACS Mesh V0.9.3 deterministic-parallel
  Delaunay candidate.
- The package contains a complete final source tree, a modified-files overlay
  and a standalone Git patch.
- `BASE_FILE_SHA256SUMS.txt` pins every pre-existing V0.9.3 file touched by
  the patch.

## Clean patch verification

Completed successfully:

1. Copied the authoritative V0.9.3 source to a clean directory.
2. Verified all pinned base-file SHA-256 values.
3. Ran `git apply --check` on the V0.10 patch.
4. Applied the patch.
5. Compared the patched tree byte-for-byte with `source/`.
6. Independently applied `overlay/` to another clean copy.
7. Compared the overlay result byte-for-byte with `source/`.

Both comparisons were empty.

## Packaged-source build

The exact `source/` directory in this package was configured and compiled with:

```text
Clang 17
C++20
Release
-Wall -Wextra -Wpedantic -Werror
```

The full static production library, V0.10 CLI, retained benchmarks and focused
probes compiled successfully.

## Packaged-source tests

```text
paacs.v010.pre_occt_core          PASS
paacs.v010.reusable_mapped_hex    PASS
paacs.v010.retained_v8_release    PASS
paacs.v010.delaunay_goldens       PASS
--------------------------------------
4/4 PASS
```

A package-built CLI then generated and audited the unit-box classified Tet4
mesh with the expected hash `16723556976053534829`.

## Sanitizers

The complete production static library was built with Clang ASan+UBSan at
`-O0 -g`, and these package-equivalent focused routes completed under:

```text
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1
UBSAN_OPTIONS=halt_on_error=1
```

```text
paacs.v010.pre_occt_core     PASS
paacs.v010.delaunay_goldens  PASS
```

No ASan, UBSan or LeakSan report was emitted.

## Installed consumer

The package-built library was installed to a clean prefix. A separate CMake
project found `PAACSMeshV010`, linked `PAACS::mesh` and built canonical box
curves successfully:

```text
0.10.0-final-pre-occt-candidate 12 2052521375206730775
```

A separate C-language consumer also called the narrow C ABI and produced:

```text
status=0 nodes=26 cells=115 error=
```

## Additional focused evidence

- Prepared-surface allocation probe: passed, zero escaped failures.
- Volume-engine probe: passed.
- T3-A refinement probe: one committed insertion, zero per-insertion global
  scans, zero context rebuilds, structured incomplete result.
- T4-A and T4-B route probes: passed their focused non-fault runs.
- One/four-thread box output: identical archive and strict hash.
- Tet10 box audit: zero invalid Tet10.
- Cylinder boundary audit: passed with zero inverted cells.

## Explicitly not claimed

- Complete historical project test umbrella.
- Full allocation-failure sweep for every newly added public wrapper.
- TSan for the new orchestration layer.
- Windows/macOS compiler qualification.
- Real-CAD/OCCT qualification.

## Archive integrity

The final ZIP and tar.gz were both extracted independently. Each extracted
package verified every entry in `SHA256SUMS.txt`; 309 internal files were
checked. The extracted ZIP and tar trees were byte-identical, and no build
directory was present under `source/`.
