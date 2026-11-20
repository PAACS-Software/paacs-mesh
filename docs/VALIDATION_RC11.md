# PAACS Mesh 1.0.0-rc.11 qualification

Validated 2026-10-02 UTC on Linux x86_64 with GCC 16.2.0, C++26, CMake 4.4.3,
Ninja 1.13.2, strict floating point and library warnings treated as errors.
Portable CPU flags; no IPO. The container allows 8 CPU equivalents and 8 GiB RAM.
This remains a release candidate, not a qualified general CAD or solver product.

The input is the rc.10 ZIP identified in baseline.json (`evidence/rc11/baseline.json`).
The [rc.10 qualification](VALIDATION_RC10.md) is retained. Complete API
semantics and route coverage are in [the controls guide](SPATIAL_CONTROL_FIELDS.md).

## Implemented

- Immutable scalar fields: point, segment, finite triangle, box and ball distance;
  linear/smooth Threshold; Minimum/Maximum; Scale/Clamp; global bounds/factor.
  Validate all input nodes, prune unused nodes, bound retained memory and evaluate
  without heap allocations.
- Shared `BoundedSizeField` contract and conservative whole-box minima for adaptive
  Tet4 generation, with deterministic sizing-work accounting. Existing region
  fields retain their behavior; tiny controls are not lost to sparse sampling.
- Global/per-body fields for constrained geometry Tet4, explicit incompatible-control
  refusals and optional bounded interior refinement.
- Independent linear-volume size audit: conservative and sampled ratios, warning
  or strict policies, aggregate audit budgets and callback-failure diagnostics.
- Tests, installed example, reproducible benchmark and actual mesh sections.
  Uniform defaults and archive v1 remain intact. Rebuild C++ clients after public
  request/report layout and field-base-class changes.

## Completed gates

| Gate | Observed result |
|---|---|
| Full Release CTest | **38/38 passed**, complete run in **270.80 s** |
| Inherited unit suite | **475 tests passed**, inside that full run |
| New field oracles | **62,342 independent checks, 27 compilation refusals** |
| New spatial meshes | **19 configurations, 23 refusal/diagnostic cases** |
| Corrected allocation harness | Passed separately in Release; coverage detailed below |
| Debug ASan + UBSan | **6/6 passed**, complete run in **157.42 s**; leak detection disabled |
| Installed C/C++ examples | **11/11 passed**; installed and qualified Release library hashes match |
| Unchanged rc.10 outputs | **207 hashes match**: 87 primitives, 64 complex profiles, 42 geometry, 10 adaptive, 4 layered |
| Isolated benchmark | **60 fresh-process samples, 20 summaries**; repeats and one/four-thread hashes match |
| Visual provenance | **10/10 gallery archives** match final Release hashes; eight unique meshes illustrated |
| Source identity | **327 C/C++/CMake/version inputs** match the final frozen manifest |

Counts are checks/configurations, not thousands of distinct CAD models. Eight
new base fields act on a unit cube; four scaled point/segment variants and seven
geometry-route configurations bring the new mesh count to 19. The inherited
engineering profiles, holes, mapped/swept choices and general-route tests also run.

Field oracles use independent long-double geometry calculations, five scales
from 1e-9 to 1e9, rotated/translated frames, increasing/decreasing thresholds and
composed fields. Random boxes/interior samples check bounds; a tiny region
between corners checks against sampling loss. Other checks cover malformed
dependencies, invalid geometry, exact/one-below budgets, moved-from fields and
concurrent deterministic queries.

Mesh checks independently cover divergence-theorem volume, reference boundary
incidence, topology, face tags and twenty barycentric samples inside each
adaptive tetrahedron. One/four-thread archives match. Refusal checks cover strict
size policy, incompatible methods, per-body overrides, aggregate multi-body work,
cancellation and rollback. Standalone audits cover Tet4/Hex8/Prism6/Pyramid5,
reject Tet10 volume cells and handle nonfinite/throwing fields without terminating.

Sanitizers cover new fields, spatial integration and allocation faults, plus
existing adaptive boxes, geometry controls and complex controls. The entire
inherited unit suite was not also run under sanitizers. Leak detection was
disabled; no leak-freedom claim is made. Correctness runs overlapped some other
build/test work, so their durations are not isolated performance measurements.

## Allocation coverage and failures found

The compiler fixture has **one** observed ordinary allocation; the adaptive
fixture has **17**. Injection exhausts those observed paths. One thousand
point/whole-box evaluations and a complete size audit also succeed with ordinary
allocations disabled.

The inherited constrained-Tet4 pipeline has **191,938** observed ordinary
allocations on the small constant-field cube. Its test uses **460 deterministic
prefix/tail/stratified settings**: 330 rejected calls retain the original mesh;
130 successful calls produce the exact baseline archive. This is sampled, not
exhaustive, general-route coverage. Some temporary-buffer failures recover.
The harness covers throwing/nothrow `operator new` and `new[]`; aligned/external
allocators are outside it. The high general-route allocation count needs profiling.

The first sanitizer run found an allocator mismatch in the **test harness**:
standard-library nothrow allocations reached ASan's default allocator, while
deallocation reached the harness replacement. Matching nothrow new/delete
overloads fixed it with mismatch checking still enabled. The corrected harness
passed separately in Release and in the final complete six-suite sanitizer run.
This test is the only source changed after the initial full Release build;
production sources remained unchanged. Both source snapshots are retained.

An early library build also had a zero-byte generated object. Removing that
object and rebuilding fixed the link failure. An initial probe initializer was
corrected to pass an owned vector to a span API. Initial allocation experiments
hit a 20,000-setting cap, prompting measurement and the explicitly stratified
test. One gallery call used incorrect arguments; the corrected `--gallery` call
and full CTest passed. Failed logs remain, clearly distinguished from final gates.
No library warning or sanitizer check was disabled to obtain a pass.

## Isolated benchmark

Three fresh-process repeats per field and requested thread count. These are
**one-thread** adaptive Tet4 medians. The timer covers `meshAdaptiveBox`, including
mandatory boundary auditing. It excludes field compilation, the extra standalone
size audit, serialization and process startup. No project builds or correctness
tests ran concurrently with these measurements.

| Field | Tetrahedra | Median call ms | Output MiB | Planner working MiB |
|---|---:|---:|---:|---:|
| Point | 5,704 | 3.17 | 0.368 | 0.592 |
| Segment | 17,176 | 9.69 | 1.077 | 1.762 |
| Triangle patch | 22,484 | 13.10 | 1.403 | 2.302 |
| Box | 24,672 | 14.33 | 1.556 | 2.536 |
| Ball | 22,320 | 12.97 | 1.402 | 2.291 |
| Minimum | 24,112 | 14.01 | 1.510 | 2.472 |
| Maximum | 4,352 | 2.80 | 0.280 | 0.451 |
| Smooth + scale + clamp | 50,736 | 33.95 | 3.111 | 5.162 |
| Uniformly fine, 0.065 | 98,304 | 58.40 | 6.205 | 10.106 |

The **point** case and uniform baseline share finest target 0.065; the point
case has far target 0.4. It uses **94.2% fewer cells**, **94.6% lower median
mesh-call time** and **94.1% smaller logical output payload**. This is the benefit
of localized refinement on this fixture, not an advantage over another mesher.
Minimum mean-ratio quality is 0.560 locally and 0.800 in the uniform baseline;
this is not an equal-solver-error or equal-quality-distribution comparison.
The baseline's minimum is not matched to every other field in the table.

Field payloads are 240 bytes for two-node controls, 600 for two-influence controls
and 480 for the smooth example. Mesh working/output payloads exclude allocator
metadata and are **not peak process RSS**. All eight adaptive fields pass the
conservative size checks and have minimum mean-ratio quality at least 0.560.
Four requested threads show mixed timing changes; no general parallel-scaling
claim is made. Raw repeats/ranges/work/hashes are in
spatial-scaling.json (`evidence/rc11/spatial-scaling.json`).

## General-geometry frontier

A field applied to boundary generation does not guarantee interior compliance.
The `.3/.6` point-field cube has a valid boundary-conforming mesh but conservative
size ratio **3.3347**, with 44 cells above warning ratio 1. The audit exposes this.

For a `.4/.6` field, optional interior refinement improves the conservative ratio
from **2.5092 to 2.0753**, with one insertion (26 to 27 nodes), passing publication
limit 2.1. Tightening the refinement criterion to 1.5 returns
`BoundaryRefinementRequired` and preserves the destination. Engine sampling and
the conservative publication audit are separate criteria. The
sixteen-case probe (`evidence/rc11/refinement-probe.log`) includes these refusals.

Coupled boundary subdivision/recovery, curved trimmed-CAD distance/projection,
curvature/proximity sizing, spatial indexing, periodic/embedded constraints and
interfaces remain priorities. Prepared planar prisms do not establish curved
CFD layers joined to a general tet core. This revision does not claim to fix
the retained rc.10 pocket/junction probes; the tapered-duct gate still passes.

OCCT, Windows/macOS, GPUs, million-element general-CAD scaling and solver
convergence were not qualified here. Use the explicit API coverage
matrix and [roadmap](GENERAL_MESHER_ROADMAP.md) to judge the next stage.

## Reproduce and inspect

Use the root README, the [API reproduction commands](SPATIAL_CONTROL_FIELDS.md#reproduce)
and the evidence index (`evidence/rc11/README.md`). Images are exact archived-cell
sections with independently measured volume/quality. Heatmaps are display samples;
the mesher uses conservative whole-box bounds rather than those samples.
