# PAACS Mesh 1.0.0-rc.10 qualification

Validated 2026-10-02 UTC on Linux x86_64 with GCC 16.2.0, C++26, CMake 4.4.3,
Ninja 1.13.2, strict floating point and library warnings treated as errors.
Portable CPU flags; no IPO. The container allows 8 CPU equivalents and 8 GiB RAM.
This is a release candidate, not a qualified general CAD or solver product.

The input is the saved rc.9 source archive identified in
baseline.json (`evidence/rc10/baseline.json`). The
[previous qualification](VALIDATION_RC9.md) is retained.

## Implemented

- Reused capped global/local axial profiles on geometry-model polygonal sweeps
  and mapped blocks. Select either cap, internal bands or both ends; limit
  local size and growth independently from the transverse target.
- Added conforming Tet4 output from convex Prism6/Hex8 cells, without new nodes.
  The public converter also accepts mixed Tet4/Prism6/Hex8 source meshes.
- Added optional quality warnings, retained strict quality rejection, and
  bounded profile planning/conversion payloads, with transactional publication.
- Exposed deeper constrained-boundary diagnostics for failed general Tet4 routes.
- Fixed roundoff in exactly filled wall-layer schedules; meaningful overrun
  still rejects. Added source finiteness/classification preflight to the converter.
- Added independent geometry/conversion tests, exhaustive observed allocation
  fault paths, an installed SDK example, reproducible benchmarks and real mesh cuts.

Contracts and limitations are in the
[API guide](COMPLEX_GEOMETRY_CONTROLS.md). Uniform defaults are retained.
C++ request/report layouts changed: rebuild clients. Archive format remains v1.

## Completed validation gates

| Gate | Observed result |
|---|---|
| Full Release CTest | **35/35 passed in 362.73 s**, in one complete run |
| Inherited unit suite | **475 tests passed**, within that full run |
| New complex geometry corpus | **64 mesh configurations, 24 refusal/rollback cases** |
| New convex-cell converter | **48 scaled/permuted mixed meshes, 21 refusal/rollback cases** |
| New allocation injection | **6,044 observed allocation points** across six fixtures, each failure retaining the prior destination |
| Existing primitive matrix | **87/87 output hashes unchanged** from rc.9 |
| Existing primitive CLI | **74 checks passed** |
| Installed C/C++ consumers | **10/10 passed**, installed library hash matches the qualified Release library |
| Debug ASan + UBSan | **7/7 passed in 70.54 s**, in one complete run; leak detection disabled |
| Isolated benchmark | **168 fresh-process samples, 56 summary rows**; all repeats and 1/4-thread output hashes match |
| Mesh image provenance | **56/56 gallery archives** match final Release hashes; all **11 unique illustrated meshes** match their recorded SHA256 |
| Source identity | **316 C/C++/CMake/version inputs** unchanged from the frozen source manifest |

These are suite/configuration counts, not a count of unrelated CAD models.
The complex corpus covers seven base profiles: I-beam, six-bore block, star duct,
twelve cooling channels, notched frame, thin bracket and mapped block. Four
profile modes and two element choices give 56 cases; eight rotated/scaled
variants of two profiles bring the total to 64. Tests independently check
volume by the divergence theorem, full reference boundary incidence, Euler
characteristic, CAD classifications, whole axial interval size caps, adjacent
growth and identical output with one/four requested threads.

The converter corpus permutes node IDs on connected mixed-cell meshes across
three scales. Checks include shared quad diagonals, source exterior orientation,
planarity, convexity, positive volume, boundary audit, in-place operation,
exact budget admission, one-byte-below rejection, cancellation and rollback.
Five-scale regressions exercise wall-layer roundoff from 1e-9 to 1e9.

Allocation injection covers profiled Hex8 (245 points), profiled Hex8 to Tet4
(255), direct hex conversion (10), profiled Prism6 (2,757), profiled Prism6 to
Tet4 (2,767) and direct prism conversion (10). It exhausts allocations observed
for these fixtures; it does not establish every possible allocation path.

Sanitizer suites: capped profiles, layered patches, existing geometry controls,
new conversion, complex geometry, complex allocation failures and tapered duct.
The complete Release suite was run; the entire inherited unit suite was not
also run under sanitizers. Leak detection was disabled, so no leak-freedom claim
is made. Correctness test times above are not isolated speed measurements.

## Benchmark observations

Three repeats per fixture/mode/family/thread setting; median mesh-call times.
No project builds or tests ran concurrently. The following is the **one-thread
Tet4** subset. Local refinement is at the start cap; the baseline uses fine
spacing along the entire axis. Both retain the same transverse resolution.

| Geometry | Cells: local / fine | Fewer cells | Call ms: local / fine |
|---|---:|---:|---:|
| i beam | 2,730 / 7,140 | 61.8% | 11.08 / 14.53 |
| multi bore | 5,220 / 12,528 | 58.3% | 12.78 / 18.93 |
| star duct | 4,224 / 10,368 | 59.3% | 11.31 / 17.69 |
| cooling channels | 4,914 / 10,920 | 55.0% | 13.78 / 18.93 |
| notched frame | 4,890 / 10,269 | 52.4% | 15.16 / 18.82 |
| thin bracket | 1,215 / 2,565 | 52.6% | 13.08 / 13.64 |
| mapped block | 2,310 / 5,460 | 57.7% | 1.51 / 3.38 |

Across both element choices the local profiles use 52.4–61.8% fewer volume cells.
One-thread Tet4 median call savings range from 4.1% to 55.3%; the benefit varies
because source-surface preparation and other fixed work remain. Some differences
are small relative to timing variability. Four requested threads provide no
consistent benefit on this small corpus; no strong parallel-scaling claim is made.

The timed call includes geometry analysis, profile planning, generation,
conversion, requested quality assessment and boundary auditing. Model/profile
construction, archive serialization and process startup are outside the timer.
Reported working/output payloads are in the raw JSON and exclude allocator
metadata. Conversion holds the native and converted meshes at the same time.

Fresh-process peak RSS from `wait4` ranges from 12,928 to 13,312 KiB. It includes
startup/inherited launcher high-water and is floor-dominated for these small
meshes; it is **not** a measurement of mesh heap allocation or evidence of a
memory saving. Process-RSS limits and million-element scalability are not
established by this corpus.

## Failures found and handled

The first mixed-cell fixture exposed wall-layer summation roundoff: three
0.2-thickness layers could be rejected for a depth of 0.6. The planner now accepts
only a small relative roundoff excess, clamps the terminal boundary exactly,
and still rejects real overrun. Converter preflight now rejects nonfinite or
invalidly classified source nodes before calling the boundary auditor, making
the reported reason predictable. The original failing logs are retained.

Generated build artifacts also had recoverable problems: two zero-byte Release
object files, two zero-byte sanitizer object files, and missing execute bits on
some sanitizer/benchmark ELF binaries. The affected generated objects were
rebuilt and executable permissions restored; source, checks, compiler warnings
and test gates were not weakened. Initial failures, recovery records, and the
subsequent complete passing runs are retained. Some earlier intermediate build
logs are incomplete; the final gate logs/status files are authoritative.

## General-meshing frontier and remaining limits

A separate tapered hollow duct passes the general constrained-Tet4 route,
including independent volume, boundary and classification checks. The dedicated
CTest gate requires a published mesh. The pocketed-block and branching-junction
probes instead refuse with `NonStarEdgeCavity` during edge discovery, preserving
the destination. Their exit code zero indicates a completed probe, not meshing
success. The exact fixtures, target size 0.7, work limit 500,000 and stage codes
are recorded. The pocket/junction boundaries include coplanar subdivisions;
these observations apply to those configurations, not every such shape.

The new control is axial on recognized straight sweeps. It is not arbitrary
3D isotropic sizing, side-face refinement on any solid, automatic CAD decomposition
or a general all-hex route. The new converter requires convex cells with planar
faces, one region and a complete oriented exterior. It does not detect global
cell overlap or add general geometry-model Tet10 support. Primitive Tet10 remains.

The optional OCCT adapter was not built in this run. Curved trimmed STEP models,
healing/defeaturing, assemblies, curved CFD layers joined to a tet core, solver
convergence, Windows/macOS and GPU behavior are not qualified by these results.
The [prioritized roadmap](GENERAL_MESHER_ROADMAP.md) starts with the observed
boundary-recovery failures, then curved CAD, shared spatial sizing, interfaces
and wall layers, followed by solver and cross-platform release qualification.

## Reproduce and inspect

Use the root README reference build and `ctest --test-dir build --output-on-failure`.
The [API guide](COMPLEX_GEOMETRY_CONTROLS.md#reproduction) provides commands
for the corpus, failure probes, gallery and benchmark. The complete
evidence index (`evidence/rc10/README.md`) distinguishes final gates, experiments
and initial failures. Build flags, source hashes, raw JSON and logs are included.
