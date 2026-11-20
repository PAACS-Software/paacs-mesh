# Geometry qualification and strategy selection

`algo::meshGeometryModel` inspects an immutable `geo::Model` and creates a
private result. Recognition uses topology, declared analytic types and endpoint
correspondence. Sampling a surface does not certify it as a plane. The result
replaces the caller's mesh only after the requested audits pass.

| Recognized geometry | Automatic route | Boundary cells |
|---|---|---|
| Affine hexahedral block | Direct structured Hex8 | Quad4 |
| Straight polygonal extrusion, including concavity and through-holes | Mesh one cap, translate through layers as Prism6 | Tri3 caps, Quad4 sides |
| Other admissible closed geometry with supported charts | Existing constrained Tet4 engine | Tri3 |
| Separated bodies with disjoint bounding boxes | Choose independently for each body, append with separate ownership | Per-body boundary |
| Touching/overlapping bounding boxes | Refuse with `UnresolvedBodyInterface` | Requires an explicit conformal hybrid plan |

This is recognition of whole bodies, not automatic decomposition of an
arbitrary solid into blocks and sweeps. A general Tet4 proposal is not a promise
that recovery will succeed. Existing explicit hybrid, layer and transition
interfaces remain available for prepared regions. Periodic/seam charts are
refused by this new route.

## Dependency-free geometry

Include `paacs/geometry/PolyhedralModel.hpp`. `PolyhedralInput` contains vertices
and faces; each face contains an outer vertex loop followed by oppositely wound
hole loops. Loops omit the repeated closing vertex. Indices are zero-based;
the built model assigns positive entity tags.

`makePolyhedralModel` validates finiteness, index validity, planar faces, loop
intersections, hole containment, two-sided edge incidence, consistent winding,
connected shell topology and nonzero signed volume. Globally inward shells
are normalized. It does not certify arbitrary inter-face self-intersections,
perform tolerance-based sewing, or repair defective CAD. The convenience
`extrudedPolygon` constructor produces closed shells from a profile and a
translation vector. Its allocation exceptions follow ordinary C++ conventions;
the validating builder returns structured failure outcomes.

The model owns its geometry. Keep it alive and immutable for the entire
meshing call. The analytic backend contract must be truthful: `Line` and
`Plane` declarations enable specialized construction. Unknown remains the
safe default for third-party backends. Source/target face, vertical edge and
side face classifications are preserved by recognized sweeps.

## C++26 usage

```cpp
#include <paacs/algo/GeometryMesher.hpp>
#include <paacs/geometry/PolyhedralModel.hpp>

std::vector<std::vector<paacs::Vec2>> rings{
    {{0,0},{3,0},{3,3},{0,3}},
    {{1,1},{1,2},{2,2},{2,1}}};
auto geometry = paacs::geo::makePolyhedralModel(
    paacs::geo::extrudedPolygon(rings, {0.2, 0, 1}));
if (!geometry.built()) { /* inspect geometry.reason */ }
else {
    paacs::algo::GeometryMeshingOptions options;
    options.targetSize = 0.25;
    options.minimumElementQuality = 0.05;
    paacs::mesh::Mesh mesh;
    auto result = paacs::algo::meshGeometryModel(*geometry.model, options, mesh);
    if (!result.published) { /* inspect result.reason; mesh is unchanged */ }
}
```

`analyzeGeometryStrategies` reports proposed routes without meshing. Mapped
block counts are exact. Sweep counts are conservative source estimates; exact
output counts are checked after cap construction. Reports include actual
output sizes, timings and underlying engine diagnostics. Geometry hashes are
diagnostic fingerprints, not cryptographic identities or universal CAD cache
keys. Meshing re-analyzes the input rather than trusting a caller-edited report.

## Body sizing and graded layers

`GeometryMeshingOptions::bodyControls` contains optional controls keyed by the
model's solid tag. Tags must exist and be unique. `targetSize = 0` inherits the
global size; an unset `strategy` inherits the global strategy. An explicit
`GeometryStrategy::Auto` requests automatic selection for that body. Bodies
still need disjoint bounding boxes; this does not build touching interfaces.

```cpp
paacs::algo::GeometryBodyControl control;
control.solidTag = geometry.model->solids().front()->tag();
control.targetSize = 0.25;  // cap resolution, and other block directions
control.layers = 12;      // exact count along the recognized translation
control.endToStartLayerWidthRatio = 0.25; // last layer is one quarter of first
control.maximumAdjacentLayerRatio = 1.3;
options.bodyControls = {control};
```

The layer ratio is the **last width / first width**, not the ratio of adjacent
widths. A ratio below one refines towards the target; above one refines towards
the source. The adjacent-width bound applies in both directions. Uniform is
one. A single layer cannot realize a nonunit ratio and is refused. Numerical
collapse or excessive growth also causes refusal. An exact layer count can
produce high aspect ratios: use a quality floor suitable for your solver.

`sourceFaceTag` chooses a starting cap and may reverse the automatic direction.
Use `analyzeGeometryStrategies` to inspect source/target tags and sweep vectors.
For mapped blocks this controls the W axis; U/V orientation remains determined
by CAD topology. `layers = 0` derives the count from size. Counts and ratios
are reported per body. These controls require a verified affine block or
straight sweep; the Tet route refuses layer controls instead of ignoring them.
Entity tags are local to a CAD import and should not be reused across changed
files without resolving them again. The ordinary no-control path is retained.

## CAD import and command line

Configure `PAACS_WITH_OCCT=ON` with Open CASCADE 8+. Link `PAACS::occt` and call
`geo::importPlanarCad(path, options)`. OCCT parses and checks the shape; the
adapter accepts one closed valid solid with planar faces and straight edges,
extracts ordered wires, validates the polyhedral representation, and compares
its volume with OCCT's volume integral. The result owns no OCCT objects.

Assemblies, curved surfaces/edges and extra non-solid geometry are refused.
The adapter does not preserve STEP product names, colors or assembly metadata;
entity tags are local to an import. STEP transfers explicitly to millimetres.
`targetSize` uses the resulting model units. BREP carries
coordinates without this API imposing a physical unit. No custom healing is
requested; OCCT's normal STEP transfer behavior still applies.

```sh
paacs-cad-mesh part.step 0.3 part.paacsmesh --threads 4
paacs-cad-mesh part.brep 0.3 part.paacsmesh --strategy sweep --min-quality .05
paacs-cad-mesh part.step 0.3 layered.paacsmesh --layers 12 --layer-ratio .25
paacs-cad-mesh part.step 0.3 bounded.paacsmesh --max-audit-bytes 67108864
```

CLI output reports strategy, cells, CAD/mesh volumes, minimum quality, logical
output bytes, source/target face tags, layer controls and audit time/storage. `--max-nodes`, `--max-cells` and `--max-bytes` lower output limits.
CAD input byte limits apply before parsing; they are not an OCCT parser memory
or time sandbox. Archive writes stage an exclusive file beside the destination
and replace it only after the complete write and close succeed.

`--report-json` emits one structured record on success or refusal, including
input topology, unit receipts, publication flags and separate import,
analysis, generation, audit, quality and archive timings. Generation includes
the attempted work of a failed/refused body and currently combines surface and
volume work. A positive `generation_seconds` does not imply publication.
`--max-file-bytes` and `--max-work` expose
the corresponding input and meshing budgets. `--size-factor F` chooses a
target of `F` times the imported solid's bounding-box diagonal.

General Tet failures retain nested surface status/reason, the first stopped
face's artifact persistent ID, processed/unattempted face counts and Steiner
recovery codes in `GeometryBodyReport`. The JSON adds `surface_reason_name` and
`recovery_reason`; for example, `candidate_provider_error` with
`steiner_disabled` identifies a boundary constraint that the current surface
route cannot recover without adding prohibited boundary points. The artifact
persistent ID is distinct from a source/target CAD face tag. A successful planar
import does not certify that every subsequent surface/volume recovery is supported.

STEP transfer uses explicit millimetres and reports the source unit names.
BREP contains no declared physical unit and retains its coordinate units.

`scripts/run_cad_corpus.py` runs deterministic selections in isolated sequential
processes, records input/executable hashes, repeated timings and Windows peak
working set, and enforces a timeout. Its reports distinguish unsupported
geometry, malformed input, exhausted budgets, process failures and published
archives. Process memory is measured; input/output budgets are not an OCCT
process-memory limit. See the local validation record for the measured ABC sample.

## Bounds, audits and performance

- Nodes, cells and logical mesh bytes are checked with overflow-safe arithmetic.
  These are not RSS limits: the input, old destination, private output, source
  triangulation, adapter topology and audits can coexist temporarily.
- `maximumAnalysisWork` bounds counted topology/recognition work.
  `maximumSourcePoints` conservatively bounds surface candidate preparation.
  `maximumMeshingWork` bounds planar constraint insertion and optional Tet
  refinement/Tet10 phase counters, and checks an exact node-plus-cell
  construction floor for each mapped/swept body. Baseline constrained Tet4
  construction has no complete work counter; this is not a total-work cap.
- `maximumAuditFaceUses` and `maximumAuditTemporaryBytes` bound the independent
  boundary audit before allocating either temporary array. The audit groups
  by minimum node and sorts exact cyclic keys inside each group; it checks
  exterior bijection, interior orientation, duplicates and nonmanifold incidence.
  It preserves quad cyclic adjacency, rather than comparing vertex sets alone.
  It is not a general geometric self-intersection or Jacobian proof.
- Specialized engines certify their affine/straight construction and switch
  to emitted-cell geometric checks near coordinate-resolution limits. Direct
  mapped/sweep entry points check count budgets before layer-table allocation. The Tet
  route retains its existing geometric validation. `minimumElementQuality > 0`
  additionally requires a minimum SubTetMeanRatio before publication; its
  distribution audit allocates O(volume cells) and sorts quality values.
  Thin valid cells can still have poor isotropic quality. This pass adds
  rejection, not automatic sliver removal to a guaranteed quality floor.
- Cancellation is cooperative at analysis/stage boundaries and inside the
  boundary audit's scans/bucket loop. A single bucket sort is not interruptible.
  OCCT parsing and some legacy kernels still lack inner-loop interruption.
- The worker pool retains its one-external-dispatcher contract. Simultaneous
  independent meshing sessions are not newly qualified. GPU acceleration is
  not implemented by this change.

The planar recovery routine starts from the existing exact-predicate Delaunay
candidate and uses deterministic edge flips to insert protected segments.
No boundary Steiner points are added. It refuses crossing constraints,
unsplit collinear vertices, unrecoverable strips and exceeded work bounds.
The output is a constrained triangulation; it is not re-legalized to a
constrained Delaunay optimum. Existing surface topology/ownership checks still
run before publication. The old public surface call keeps its prior default;
the overload with `ProductionSurfaceOptions` enables the new path explicitly.

## Reproduce the qualification

```sh
ctest --test-dir build/gcc16-release -L geometry --output-on-failure
python3 benchmarks/run_geometry_corpus.py \
  build/gcc16-release/tests/paacs-geometry-corpus --output geometry-results
build/gcc16-release/tests/paacs-occt-corpus --fixtures generated-cad
python3 benchmarks/run_audit_comparison.py \
  build/gcc16-release/tests/paacs-audit-benchmark --output audit-results
```

`tests/geometry_fixtures.hpp` defines independent analytic volumes. Tests check
positive oriented cell volume, boundary incidence, boundary entity positions,
trim containment, complete face coverage, strategy choice and stable output
across thread counts. CAD fixtures are constructed by OCCT and round-tripped
through both formats. They are controlled regression models, not an external
industrial CAD corpus. The support matrix records the qualified cases.
