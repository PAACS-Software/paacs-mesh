# PAACS Mesh 1.0.0 — release verification

Official source release finalized 2026-10-02 UTC. Reference configuration:
Linux x86_64, GCC 16.2.0, C++26, CMake 4.4.3, Ninja 1.13.2, Release,
strict floating point, warnings as errors, portable CPU flags, no IPO or OCCT.

This finalization changes release identity, packaging and documentation. Every
file under `src/` and `tests/` is byte-identical to the qualified input package.
Across the 327 retained C/C++/CMake/version inputs, only CMake/VERSION metadata,
the SDK version string and one C API comment changed. The
source comparison (`evidence/1.0.0/source-promotion.json`) records the input hash.
Original historical evidence is unchanged.

## Fresh checks on 1.0.0

| Check | Result |
|---|---|
| Clean Release library and both CLIs | Built successfully with GCC 16 / C++26 |
| Selected release CTest checks | **4/4 passed in 4.50 s**: CLI version, hardening, control fields and spatial controls |
| Installed C/C++ examples | **11/11 passed in 0.45 s**; installed library hash matches the fresh Release library |
| Version entry points | VERSION, CMake, SDK, CLI, C API, installed package and CPack identify **1.0.0** |
| Field checks within fresh run | **62,342 independent checks, 27 compilation refusals** |
| Mesh checks within fresh run | **19 configurations, 23 refusal/diagnostic cases** |
| Output continuity | **14/14 recorded mesh hashes** and **10/10 gallery archive SHA256 hashes** match the input release |
| Current visuals | Eight independently measured meshes re-rendered with final release labeling |

See fresh release evidence (`evidence/1.0.0/README.md`). Durations above are
correctness checks, not isolated comparative speed measurements.

## Retained full qualification

The [original implementation qualification](../docs/VALIDATION_RC11.md) remains
included with its actual version labels and logs:

- 38/38 full Release suites, including 475 inherited unit tests.
- Six complete ASan/UBSan suites; leak detection disabled.
- Exhaustive observed compiler/adaptive allocation injection and explicitly
  stratified general-geometry injection.
- 207 unchanged earlier mesh-output hashes.
- Isolated local-versus-uniform sizing measurements and documented quality/memory
  limits; these are not external-mesher benchmarks.

Those complete runs were performed on the input implementation. The full inherited
suite and sanitizers were **not rerun in full for this metadata-only promotion**.
The fresh clean build, installed consumers, selected integration checks and source
comparison above verify final packaging without relabeling old tests as new runs.

## Scope and limits

The [release policy](../docs/RELEASE_POLICY.md) and API guides define the 1.0.0 scope.
Protected-boundary refinement, general curved CAD, proximity/curvature-driven
sizing, periodic/embedded constraints and curved CFD layers remain development
areas. General Tet4 boundary sizing does not by itself guarantee interior size
compliance; use the achieved-size report and the application's acceptance policy.

OCCT integration/qualification and matched comparisons are deferred as requested.
No additional Windows/macOS, GPU, solver-convergence or industrial-CAD coverage
is claimed by assigning the official version.

The supplied source has no distribution license; see [licensing status](../LICENSING_STATUS.md).

## Reproduce

Use the root README build and installed-consumer instructions. To reproduce the
selected promotion checks after building the relevant targets:

```sh
ctest --test-dir build/gcc16-release --output-on-failure \
  -R '^(paacs\.cli\.version|paacs\.hardening|paacs\.sizing\.(control_fields|spatial_controls))$'
python3 scripts/verify_version.py --final
python3 scripts/package_source.py --output out/paacs-mesh-1.0.0.zip
```

The packager validates version agreement, all ZIP CRCs and each archived file's
SHA256 against its manifest, then writes an external archive checksum. Historical
figures/logs retain their labels; current figures use the final version and the
freshly reproduced, byte-identical meshes.
