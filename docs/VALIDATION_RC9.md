# PAACS Mesh 1.0.0-rc.9 qualification

Reference: Linux x86_64, GCC 16.2.0, C++26, CMake 4.4.3, Ninja 1.13.2,
strict floating point, library warnings treated as errors, portable CPU flags.
Container allowance: 8 CPU equivalents and 8 GiB RAM. This remains a release
candidate; general CAD and solver readiness are not established.

The latest saved rc.8 ZIP was restored after the workspace reverted. This
revision was implemented and validated from that archive; unpublished results
from the interrupted earlier attempt are not used as evidence here. Exact
input identity is recorded in baseline.json (`evidence/rc9/baseline.json`).
The prior [rc.8 report](VALIDATION_RC8.md) is retained.

## Implemented

- Shared immutable global/local axial size profiles with overlapping controls,
  finite bands, cap selection, bounded compilation and distribution planning.
- Whole-interval axial size caps, independently checked adjacent growth,
  finite-domain closure and constant-plateau compaction.
- Profile consumers for Cylinder/Tube Tet4/Tet10/Hex8 sweeps and prepared
  Tri3/Quad4/mixed prism/hex layers.
- Opt-in varying-section solid-cylinder Tet4/Tet10 generation and a reusable
  conforming Tet4-to-Hex8 converter with explicit source exterior-winding checks.
- Reuse of identical transverse slab patterns, with cached/uncached diagnostics,
  complete world-space geometry validation and unchanged resource admission.
- Public installed SDK example, independent geometry/topology/field tests,
  exhaustive allocation-failure paths, measured benchmarks and actual mesh cuts.

The [API guide](REUSABLE_SIZE_CONTROLS.md) defines contracts and limitations.
Existing uniform defaults are preserved. C++ request/report structs gained
fields, so consumers must rebuild. The mesh archive format remains v1.

## Validation gates

| Gate | Observed result |
|---|---|
| Release CTest coverage | **31/31 suites covered successfully**: the full run passed 30/31, then the converter fix and all three affected suites passed targeted reruns |
| Inherited unit suite | **475 tests passed**, within the full run |
| Capped-profile oracle | **221 profile configurations, 20 tracked refusals**, plus moved-from safety checks |
| Profile consumers | **75 mesh configurations, 58 refusal/rollback cases**; independent topology/volume, full-reference high-order certificates, cached/uncached validation, 1/4-thread identity |
| Generic Tet4-to-Hex8 converter | **3 scale variants, 51 refusal/rollback cases**, including reversed exterior faces and independent outward normals |
| Existing primitive matrix | **87/87 output hashes unchanged** from rc.8 |
| Existing primitive CLI suite | **74 checks passed** |
| Installed C/C++ consumers | **9/9 passed** after refreshing the final installed library |
| Debug ASan + UBSan | **6 distinct suites covered successfully** across initial passes and post-fix reruns; leak detection disabled |
| Isolated benchmark sweep | **156 fresh-process samples, 52 summary rows**; every repeat and tested 1/4-thread hash agrees |
| Illustration provenance | **7/7 illustrated mesh archives** match final Release output byte-for-byte |

The profile corpus covers start/end/both caps, central and overlapping bands,
uniform fine controls, cylinders/tubes, all three element families, scales
1e-6/1/1e6 with rotation/translation, structured quads/triangles and mixed patches.
Independent checks include the analytic distance-field minimum over complete
intervals, boundary incidence, boundary Euler characteristic, divergence-theorem
volume, geometry tolerances, shared quadratic edge identity and reference-cell
Jacobian certification. Resource limits are checked at the reported budget and
one byte below, alongside quality gates, cancellation and invalid controls.

Allocation injection visits every observed allocation point until success.
Profile Sweep uses 14/24/17 points for Tet4/Tet10/Hex8; Volume uses 54/64/70;
profile-driven patches use 25. These fixtures exercise repeated slab reuse.
Every injected failure preserves the prior destination. Existing model,
geometry, adaptive-box and layer fault coverage is retained as well.

## Failures found and resolved

The first converter regression run exposed a real issue: a reversed exterior
Tri3 could pass the generic face-incidence audit and be converted. That audit
checks internal volume-face orientation, but does not require outward exterior
winding. The converter now explicitly compares source boundary signs with the
oriented tetrahedral faces, before publication. The malformed-source regression
and independent output-normal checks pass at all three scales. The original
failure and fix record (`evidence/rc9/converter-winding-fix.json`) are retained.

The full Release run took 356.90 s and included that failing converter test.
After the fix, converter, profile-consumer and allocation suites all passed in
39.91 s. Unchanged inherited suites were not repeated solely to obtain a clean
single-run percentage. These are correctness-run times, not isolated benchmarks.

The initial sanitizer build contained two zero-byte generated object files.
Only those generated objects were removed and rebuilt from unchanged source;
source hashes and the failed linker log are retained. The 75-mesh sanitizer
corpus initially exceeded its 300-second limit. Its sanitizer-only limit was
raised to 900 seconds without dropping checks. It then passed in 562.13 s;
the post-fix allocation suite passed in 64.10 s. The initial passed suites plus
the successful reruns cover six distinct ASan/UBSan suites. No ASan/UBSan defect
was reported in the successful runs. This does **not** establish leak freedom.

## Measured cost and trade-offs

Cylinder radius 0.6, height 4.8, global target 0.36, selected start-cap target
0.12, field growth 1.25, chord tolerance 0.015 and relative-volume tolerance
0.03. Medians of three isolated fresh processes at **one thread**, measuring
complete mesh calls including planning, validation and boundary audit.
Field compilation, hashing and optional archive I/O are outside this interval.

| Route | Volume cells | Complete call | Output array payload | Minimum quality |
|---|---:|---:|---:|---:|
| Tet4 sweep | 10,368 | 1.925 ms | 0.692 MiB | 0.288 |
| Tet4 volume | 2,973 | 10.398 ms | 0.240 MiB | 0.110 |
| Tet10 volume | 2,973 | 23.838 ms | 0.539 MiB | 0.110 |
| Hex8 sweep | 2,560 | 1.980 ms | 0.351 MiB | 0.150 |
| Hex8 volume | 11,892 | 27.004 ms | 1.644 MiB | 0.055 |

Quality is tetrahedron mean ratio, or the minimum fixed sub-tetrahedron mean
ratio for Hex8; ideal 1. The converted Hex8 example flags **10** cells below
its default warning threshold 0.1. Tet4 has three cells below 0.2, with minimum
about 0.110; its default 0.1 threshold flags none. Use stricter thresholds when
required. Converted Hex8 has four times the tetrahedron count and costs more
memory/time than the direct Hex8 sweep. It is a flexible conforming conversion,
not a claim of superior all-hex quality.

For the same Volume route, uniform fine sizing produces **18,240**
tetrahedra versus **2,973** for local refinement: **83.7% fewer**.
They are different meshes, not an equal-error solver comparison. Most small
process RSS readings hit the launcher floor near **11.75 MiB**;
the uniform fine converted-Hex case reached **24.81 MiB**.
Array payload, conservative working admission and process RSS are different
measures. The refined Volume Tet4 working admission is about
**33.82 MiB**, while its observed Delaunay arena/scratch payload
is about **96.6 KiB**.

### Slab reuse

| Tet4 case | Uncached full call | Reused full call | Measured speedup | Built / reused slabs |
|---|---:|---:|---:|---:|
| Refined start cap | 21.121 ms | 10.398 ms | 2.03× | 6 / 10 |
| Uniform fine sections | 129.872 ms | 13.233 ms | 9.81× | 1 / 39 |

All 12 cached/uncached family/thread comparisons in this benchmark have
identical mesh hashes. The reuse contract itself guarantees validated affine
slabs; it does not claim that arbitrary geometries always have the same Delaunay
connectivity under scaling. Four-thread samples are included in the raw data;
planning/Delaunay insertion remain serial, so no broad many-core claim is made.

## Scope and reproduction

The Volume adapter is for a solid straight cylinder. Tube controls use the
fixed-section Sweep strategy; prepared prism/hex profiles change normal
spacing only. The outer polygon stays common across adaptive cylinder sections.
The global/profile cap and growth guarantee constrain axial slab widths,
not every 3D edge or face-neighbor equivalent-volume ratio. A short global
plateau can remain below the cap to close the domain and honor growth.

Generic conversion requires a nonoverlapping single-region Tet4 source with a
complete oriented Tri3 exterior. It has no global intersection detector and
performs no curved CAD projection. Cross-platform builds, curved CAD import,
healing, multi-region conversion, solver convergence and general hex refinement
are unqualified here. The existing optional OCCT adapter was not enabled in
this reference build. Uniform defaults and feature-analysis opt-in behavior
remain unchanged.

See raw evidence (`evidence/rc9/README.md`), benchmark data (`evidence/rc9/profiled-scaling.json`),
[reproduction commands](REUSABLE_SIZE_CONTROLS.md#reproduce-and-interpret-evidence),
and the [public example](../examples/profiled_controls.cpp). The source ZIP includes
all source/tests/docs/rendering scripts and SHA256SUMS.txt; toolchain binaries,
build trees and regenerable mesh archives are excluded.
