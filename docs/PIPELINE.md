# Staged boundary and volume API

PAACS Mesh **1.2.1** exposes the classical construction stages through
`paacs_mesh.pipeline`. Install its matching wheel; no local SDK, compiler or
machine-learning package is needed:

```sh
python -m pip install paacs-mesh==1.2.1
# Optional coordinate/connectivity arrays:
python -m pip install "paacs-mesh[numpy]==1.2.1"
```

The first scope is uniform sizing on one original solid, a conforming Tet4
volume, Standard recovery and the Normal advancing interior constructor.
Assemblies, face/box controls, growth controls and alternative constructors are
not accepted by this staged contract. Ordinary `CadModel.mesh()` keeps its
existing controls and defaults. No learned decisions are made by this API.

## Construct and inspect stages

```python
import paacs_mesh as pm
from paacs_mesh import pipeline
from pathlib import Path

Path("run").mkdir(exist_ok=False)
model = pm.read_step("part.step", retain_source=True)
boundary = pipeline.prepare_boundary(model, size=0.4, max_seconds=300)
boundary.save("run/boundary")

# Real CAD-edge chains and the qualified oriented Tri3 surface:
boundary.edges().save("run/edges")
boundary.surface().save("run/surface")

# Classical sites are available before any constraint recovery or volume fill.
cloud = boundary.generate_initial_points()
cloud.save("run/points")
print(cloud.point_count, cloud.sha256, cloud.report)

mesh = boundary.mesh(cloud)
mesh.write_archive("run/mesh.paacsmesh")
mesh.write_vtk("run/mesh.vtk")
print(mesh.report["audit"], mesh.report["quality"])
```

Use new output directories. Boundary and point bundles contain private CAD or
CAD-derived data. Surface and edge `.save()` outputs include real VTK/archive
files with their metadata and hashes. `.write("surface.vtk")` and
`.write("edges.paacsmesh")` export a single file; use `.save()` for the complete
provenance bundle. Neither format is an AI-generated interpretation.
Boundary, point and geometry bundles are staged in a private sibling directory
before directory rename publishes the final destination. Existing destinations
are checked before publication and ordinary failures remove only owned staging
files. Atomic visibility relies on the local filesystem's rename semantics;
file flushing is not a cross-filesystem or crash-durability guarantee.

`boundary.edges()` exports retained ordered CAD-edge chains, including shared or
seam uses. It does not label every triangle edge as a CAD edge. Explicit CAD
body/face ownership and senses are available through
`boundary.cad_topology(max_bytes=...)` and export metadata. A face may have
multiple owners; incidence is taken from CAD inspection rather than inferred
from proximity. Surface coordinates, triangle winding, original tags and units
remain attached to the frozen boundary.

Preparation qualifies a source-CAD boundary and stops before interior placement.
Standalone point generation performs ordinary placement and auxiliary numerical
topology finalization, then stops **before constraint recovery and material
volume construction**. Auxiliary work can be substantial and must be timed.
The point-generation receipt is not an audited volume mesh, even when its sites
are admissible. A final `mesh()` result must pass the Native publication checks.

## Save, reload and reuse points

```python
from paacs_mesh import pipeline

boundary = pipeline.load_boundary("run/boundary")
cloud = pipeline.load_initial_points("run/points")
mesh = boundary.mesh(cloud)
```

Boundary loading snapshots bounded files, reimports the exact saved STEP and
revalidates the tagged shell against CAD. Point loading verifies a fixed manifest,
point bytes, count, report and hashes; no CAD or volume operation is performed by
point loading. Passing a loaded cloud to a boundary then checks its exact source,
boundary, build and normalized sizing/resource options before Native admission.
Changing a report or matching its checksum cannot bypass geometric checks.

Exports from reloaded boundaries and point clouds protect their saved
input members and their file aliases, including completed volume mesh writers.
Write results to separate output files or new artifact directories.
Hashes identify artifact bytes and do not attest who executed a recorded
generation; benchmark provenance must also retain its controlled execution.

The first bundle versions require the **same native-extension binary hash**.
Recreate boundaries and clouds for another package build, Python ABI or platform;
portable archive encoding alone does not establish cross-build CAD equivalence.
Published 1.2.0 and older development bundles are not silently upgraded.

`InitialPointCloud` owns immutable little-endian IEEE64 coordinates in
millimetres. Its public properties are `point_count`, `sha256`, `source_sha256`,
`boundary_sha256`, `units`, `options`, `provenance` and `report`. Dictionary
properties return independent copies. The following extraction is explicit:

```python
points = cloud.points(max_bytes=32 * 1024**2)  # tuple of three-number tuples
array = cloud.to_numpy(max_bytes=32 * 1024**2)  # optional float64 N×3 copy
surface_arrays = boundary.to_numpy(max_bytes=128 * 1024**2)
```

Array copies are initially read-only and own their storage. Editing a copy does
not edit the saved cloud, boundary or Native model. Point order is its stable
local identity: index 0 remains the first supplied site. The final mesh may use
different export node IDs; use explicit coordinate/classification evidence when
matching them. Native may add recovery points, so supplied count is distinct from
final mesh node count.

For an external point-placement experiment, keep the same prepared boundary and
pass a bounded list/array of proposed finite triples to `boundary.mesh(points)`.
No callback or model executes inside the engine. Count inference, containment
queries, rejected proposals and reconstruction costs in that experiment.

## Admission and failures

| Request | Meaning |
| --- | --- |
| `boundary.mesh()` or `boundary.mesh(None)` | Classical initial placement followed by the common constructor. |
| `boundary.mesh([])` | Explicitly supplied empty cloud; no automatic placement fallback. |
| `boundary.mesh(cloud)` | Exact-bound staged cloud, passed as owned encoded bytes. |
| `boundary.mesh(points)` | Explicit bounded sequence/array of finite millimetre triples. |

Duplicate, outside, nonfinite or unretained supplied points refuse the attempt.
The engine does not clip, replace or silently drop them. Every supplied point
must remain fixed and participate in a published Tet4. The tagged oriented
boundary, original CAD checks and complete Native audit/quality gates remain
mandatory. Shape-quality preferences are distinct from geometric validity and
solver adequacy.

A failed generation call returns no partial cloud. A failed volume call returns
no partial valid volume mesh. `MeshingError.report` retains the relevant Native
diagnostics. Existing `pipeline.initial_points(mesh_or_refusal, max_bytes=...)`
still supports the earlier replay snapshot interface; a captured empty snapshot
is distinct from an unavailable snapshot. Prefer standalone generation when a
study needs points without running the volume constructor first.

Resource options are frozen with a boundary: `max_nodes`, `max_cells`, `max_work`,
`max_working_bytes`, `max_output_bytes`, `max_audit_bytes` and `max_seconds`.
Copies and reads have their own `max_bytes` allowances. Point tuple extraction
reserves 64 KiB plus a conservative 256 bytes per decoded point; point array
extraction reserves 64 KiB plus its float64 payload. Encoded clouds are capped
at one million points and their resource/byte limits. Oversized, malformed or
identity-mismatched artifacts refuse rather than being partially loaded.

`max_seconds` is a cooperative Native deadline after the operation gate is
acquired. It does not cover STEP import, all Python file/fingerprint work or a
sequence of retries. Memory counters model named allocations rather than OS
process RSS. A benchmark requiring a hard overall allowance must supervise the
complete process and count preparation, generation, every attempt, checks and
exports. Named stage clocks are documented spans, not an additive partition of
total wall time.

`paacs_mesh.experimental` remains a compatibility alias for the same functions
and classes. The [boundary experiment reference](EXPERIMENTAL_VOLUME.md) describes
the archive contract and the original replay API in more detail.
