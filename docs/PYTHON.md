# Python STEP interface

The Python distribution calls the C++ STEP importer and Tet4 mesher in process.
Development wheels target **CPython 3.12** on these platforms:

| System | Architecture | Compatibility target |
|---|---|---|
| Windows | x64 | Qualified on Windows Server 2022 |
| macOS | ARM64 | macOS 15 or newer |
| Linux | x64 | glibc 2.28 or newer, manylinux 2.28 |

Each repaired wheel contains the shared CAD libraries and runtime dependencies.
Installation needs Python and pip, with no compiler, OCCT SDK or `PATH` setup.

Wheels are development artifacts from the
[Python STEP wheel workflow](../.github/workflows/python-wheel.yml). Choose the
artifact matching the system and architecture above.
This interface has not been published to PyPI. Download and extract the workflow
artifact, then install its `.whl` file in a Python 3.12 environment. For example,
on Windows:

```powershell
python -m pip install --force-reinstall --only-binary=:all: .\paacs_mesh-1.2.0rc1-cp312-cp312-win_amd64.whl
```

Development artifacts can share the candidate version; `--force-reinstall`
ensures an updated wheel replaces an earlier artifact with that version.

## Import, mesh and export

```python
import paacs_mesh as pm

model = pm.read_step("part.step")
print(model.info)  # CAD counts, source units, transferred units and volume
mesh = model.mesh(size=0.7)
mesh.write("part.vtk")
mesh.write("part.paacsmesh")
print(mesh.report)
```

STEP coordinates, target sizes and reported distances are in **millimetres**;
volumes are in cubic millimetres. Source length-unit names remain in
`model.info["source_length_units"]`. Meshes own their data and remain usable
after the imported model is discarded. A returned report is a copy.

`Mesh` exposes `node_count`, `cell_count`, `tetrahedron_count`, `face_tags`,
`body_tags` and `report`. Legacy VTK exports include full-precision coordinates
and both `entity_dimension` and `entity_tag` for nodes and cells. Native archive
v1 retains classifications and connectivity; it does not store the STEP source,
units or the meshing report. Export cannot overwrite the input STEP file.

The same operation is available from the command line:

```powershell
python -m paacs_mesh part.step part.vtk --size 0.7
```

## Face sizing, volume regions and growth

Use `model.face_tags` and `model.body_tags` to inspect the imported IDs before
assigning controls. Face tags belong to that imported model; a different or
edited STEP file can assign different IDs.

`FaceSize` refines a selected CAD face. `BoxSize` defines a sizing influence
region inside the existing material; it does not create another body or an
internal material boundary. These immutable controls are copied into the native
request. Overlapping influences use the smaller requested size.

For the project-generated [cylinder fixture](../tests/fixtures/cad/cylinder.step),
face 1 is the curved wall:

```python
model = pm.read_step("cylinder.step")
print(model.face_tags)  # [1, 2, 3]

mesh = model.mesh(
    size=0.7,
    face_sizes=[pm.FaceSize(tag=1, size=0.16)],
    volume_controls=[pm.BoxSize(
        lower=(-0.3, -0.3, 0.7), upper=(0.3, 0.3, 1.3), size=0.15,
    )],
    min_size=0.12,
    growth=1.55,
    assembly="conforming",
    max_work=100_000_000,
)
mesh.write("sized-cylinder.vtk")
print(mesh.report["sizing"])
```

The spatial route currently supports **one Tet4 solid with
`assembly="conforming"`**. The assembly-preserving default remains unchanged;
unsupported combinations return a `MeshingError` with a hint. No parts are
fused or dropped to accommodate a sizing request. The example explicitly raises
the logical work budget for the spatial route.

`native_sizing=None` automatically selects spatial sizing when controls are
present. Set `native_sizing=True` to use that route even without controls, for
example when comparing a uniform baseline against a refined mesh. Set it to
`False` to retain the existing analytic end-cap sizing on supported cylinders
and tubes; box controls require spatial sizing. Calls without controls retain
the existing route selection and defaults.

`min_size` clamps the requested field, rather than imposing a minimum on every
produced edge. Curvature, boundary recovery and connectivity can require smaller
edges. A control's `growth=None` inherits the global `growth`; an explicit
per-control value overrides it:

```python
wall = pm.FaceSize(tag=1, size=0.16, growth=1.25)
```

Growth must be at least 1. It sets the distance-based transition from a fine
control to the global size, rather than a strict adjacent-cell size ratio.
**Growth 1 extends that control's fine size throughout the domain.** Sizes remain
best-effort targets. Inspect achieved medians, percentiles and geometric error
when qualifying a mesh; a good median alone can conceal larger cells.

The command-line equivalent is:

```powershell
python -m paacs_mesh cylinder.step sized-cylinder.vtk --size .7 `
  --face-size 1 .16 --size-box -.3 -.3 .7 .3 .3 1.3 .15 `
  --min-size .12 --growth 1.55 --assembly conforming --max-work 100000000
```

Repeat `--face-size` and `--size-box` for additional controls. The CLI also
accepts `--native-sizing` and `--no-native-sizing`.
The [Python sizing example](../examples/python_step_sizing/README.md) compares
four requests on the same cylinder, with measured sizes and mesh sections.

## NumPy mesh arrays

NumPy is optional. Install it alongside the development wheel, or use the
`numpy` extra when installing a wheel from a local directory:

```powershell
python -m pip install "numpy>=1.26.4,<3"
python -m pip install --force-reinstall --find-links . "paacs-mesh[numpy]==1.2.0rc1"
```

Import, meshing and export continue to work without NumPy. Array extraction
loads it on demand and gives an installation hint when it is absent.

```python
import numpy as np

arrays = mesh.to_numpy()
points = arrays["points"]
tets = arrays["cells"]["tet4"]
tet_points = points[tets["connectivity"]]
body_1 = tets["connectivity"][tets["entity_tags"] == 1]

triangles = arrays["cells"].get("tri3")
if triangles is not None:
    face_1 = triangles["connectivity"][triangles["entity_tags"] == 1]

# Signed Tet4 volumes in mm^3, retaining the native node ordering.
edges = tet_points[:, 1:] - tet_points[:, :1]
volumes = np.einsum("ij,ij->i", edges[:, 0],
                    np.cross(edges[:, 1], edges[:, 2])) / 6
print(volumes.min(), volumes.sum())
```

The returned dictionary contains:

| Key | dtype | Shape / meaning |
|---|---|---|
| `points` | `float64` | `(N, 3)`, coordinates in mm |
| `node_entity_dimensions` | `int32` | `(N,)`, native owning entity dimension |
| `node_entity_tags` | `int32` | `(N,)`, native owning entity tag |
| `cells` | dictionary | One block per cell type present |
| `units` | string | `"mm"` |

Each cell block contains:

| Key | dtype | Shape / meaning |
|---|---|---|
| `connectivity` | `int64` | `(M, nodes_per_cell)`, zero-based rows into `points` |
| `entity_dimensions` | `int32` | `(M,)`, native cell entity dimensions |
| `entity_tags` | `int32` | `(M,)`, native cell entity tags |
| `cell_indices` | `int64` | `(M,)`, original zero-based native cell positions |

Block names are `line2`, `line3`, `tri3`, `tri6`, `quad4`, `tet4`, `tet10`,
`prism6`, `hex8` and `pyramid5`; absent types have no block. This describes
the array format, rather than adding mesh families to the STEP interface.
Rows within a block retain native order, connectivity and orientation.
`cell_indices` maps grouped rows back to the original mixed cell order.
Coordinates and nodes are never welded, including coincident assembly nodes.
Node classification may refer to a vertex, edge, face or body; it is not a
derived node-to-body assignment. Face and body tags remain model-local IDs.

Each call makes independent, C-contiguous NumPy-owned copies. Arrays are
initially read-only to catch accidental edits, and survive deletion of the mesh,
model or input file. NumPy can re-enable writing to an owned array; doing so
changes only that snapshot. It cannot edit the native mesh, its audit or later
exports. Use `.copy()` when preparing writable solver data.

`mesh.to_numpy(max_bytes=536_870_912)` defaults to a 512 MiB aggregate numeric
array budget. The complete payload is checked before allocation; an exceeded
limit raises `MemoryError`. The limit excludes Python object overhead and the
already-owned native mesh, and is separate from meshing budgets. It does not
bound process memory. No partially built snapshot is returned on failure.

The [array example](../examples/python_mesh_arrays/README.md) computes volumes,
selects a face and saves the arrays with their classifications.

## Recovery, assemblies and limits

`model.mesh` accepts `recovery="standard"` or the opt-in `"advanced"` recovery
mode. Both use the native publication checks. This wrapper does not expand the
engine's documented [CAD coverage](NATIVE_CAD.md).

The default `assembly="preserve_parts"` retains the original bodies and distinct
body nodes. Each part must pass its native checks. `mesh.report` records
`mesh_scope="independent_parts"`, interference pairs and their native overlap
volumes. This result does not certify a globally non-overlapping material domain.
Choose `assembly="conforming"` to request the engine's conforming-domain policy;
unsupported or unresolved interfaces are refused.

Meshing and validation each default to one worker. Native import and meshing
calls are serialized within a Python process because the engine uses a shared
worker dispatcher. The GIL is released during those calls, so other Python work
can continue. Worker counts control parallelism inside the active meshing call.
Advanced controls are keyword arguments:

```python
mesh = model.mesh(
    size=0.7,
    recovery="advanced",
    assembly="preserve_parts",
    workers=1,
    validation_workers=1,
    max_nodes=1_000_000,
    max_cells=5_000_000,
    max_seconds=90,
)
```

Other limits are `max_work`, `max_output_bytes`, `max_working_bytes` and
`max_audit_bytes`; `read_step` accepts `max_file_bytes`. Unspecified limits retain
the native defaults. `max_working_bytes` bounds the primitive-route working
budget. Byte limits cover their named native budgets, not process RSS.
`max_seconds` requests cooperative cancellation during meshing; an OCCT
operation can overrun the deadline. The deadline starts after the native
operation gate is acquired; time queued behind another job is excluded. Import
has no deadline in this interface.

`CadImportError` and `MeshingError` expose `stage`, `reason` and `report`.
`ExportError` is an `OSError` with the same fields. A refused mesh is not returned
as a successful result. Reports include native boundary/connectivity audit,
quality, per-body routes, CAD error and stage timings. Quality warnings remain
distinct from validity failures.

```python
try:
    mesh = model.mesh(0.7, max_cells=100)
except pm.MeshingError as error:
    print(error.reason, error.report)
```

## Inspect CAD and select faces

An imported model exposes copied JSON metadata for bodies, faces, edges and
vertices. Inspection is paginated (`limit` up to 1,000); each response identifies
the exact imported STEP bytes, the CAD kernel version and the units.

```python
import paacs_mesh as pm

model = pm.read_step("housing.step")
page = model.inspect("faces", limit=100, topology=True)
for face in page["entities"]:
    print(face["tag"], face["surface"], face["area"], face["owners"])

bore = model.select_faces(surface="cylinder", radius=2.5, expect=1)
caps = model.select_faces(surface="plane", normal=[0, 0, 1], expect=1)
model.define_group("bore_wall", bore, group_id="housing-bore-wall")
model.save_intent("housing.paacs.json", revision="design-1")

mesh = model.mesh(
    1.0, assembly="conforming", min_size=.15, growth=1.3,
    face_sizes=[pm.FaceSize(tag, .2) for tag in bore.tags],
    volume_controls=[pm.BoxSize([-3, -3, 0], [3, 3, 5], .3)],
)
```

Predicates combine with AND. Available filters include explicit face tags, a
companion group, a directly mapped STEP representation name, surface category,
owning body, area, analytic radius, area-centroid box, planar outward normal and
shared-edge adjacency. Radius applies to cylinders and spheres. A normal filter
applies to planes; a shared face needs an explicit owning body to define its
outward direction. Face area centroids can lie outside a trimmed face or off a
curved surface; they are not projection samples.

Geometric/name queries expect one face by default. Use `expect=2` or an explicit
`expect=(minimum, maximum)` for multiple matches. Queries using only explicit
tags, and group queries, retain their full count. `SelectionError` reports
unknown tags, incomplete required metadata, budget limits or a cardinality
mismatch. Queries scan the complete
bounded candidate set rather than choosing the first match. Returned records,
query descriptions and group dictionaries are independent copies. A selection
belongs to its imported model; applying it to another model is refused.

`read_step(..., names=True)` retains nonempty names from directly transferred,
shape-backed STEP representation items. This is not full product/XDE metadata:
product labels, colors and assembly paths are not provided. Names must match the
same placed topology; a prototype name is not guessed onto another instance.
Conflicting candidates and unmatched face-name mappings remain diagnostic.
`names=False` skips this optional metadata scan without changing meshing.

`SelectionError` and `IntentError` are `ValueError` subclasses with `stage`,
`reason` and a copied JSON-compatible `report`.

## Save groups and reopen a design

```python
reopened = pm.read_step("housing.step", intent="housing.paacs.json")
bore = reopened.select_faces(group="bore_wall")
print(reopened.groups["bore_wall"]["id"], bore.tags)
```

The `paacs-cad-intent-v1` companion file records the source SHA-256 and byte count,
kernel/tag-scheme versions, caller-assigned group IDs, queries, face tags and
geometry/orientation guards. Reopening verifies all groups before applying any
of them. Edited sources, incompatible kernels and mismatched guards are refused
by default. Tags identify entities within one imported snapshot; they do not
claim universal correspondence across CAD edits or STEP re-export.

`load_intent(path, rebind=True)` explicitly reruns independent geometric/name
queries for a changed source and reports changes while preserving group IDs.
It requires the previous exact match count and surface categories. Queries tied
to old tags, bodies, adjacency or other groups cannot be rebound automatically.
Review and preview the resulting selections before using a changed design.

The source hash is captured from the same immutable bytes passed to the parser.
Later edits to the input path cannot change a live model's identity. Ordinary
imports release those bytes; `retain_source=True` keeps them for later project
export even when the original file has been deleted. This increases retained
memory by the input size. Input/topology and metadata-output budgets do not
impose a hard process-memory limit on the CAD parser.

## Preview selections without a GUI

```python
bore = model.select_faces(surface="cylinder", radius=2.5, expect=1)
preview = model.preview(bore, linear_deflection=.02)
preview.write_svg("housing-selection.svg", title="Housing bore wall")
display = preview.to_dict()
```

The preview contains owned points, triangles and imported face tags. Selected
faces are gold in the SVG. It tessellates a deep CAD copy, checks the source
snapshot for changes, and leaves the geometry used by the mesher intact. It is
display geometry, not a conforming surface mesh or a validity/error certificate.
`linear_deflection` uses millimetres. `max_points` and `max_triangles` bound the
returned arrays; serialized display data is limited to 8 MiB. Those limits do
not impose a hard runtime or memory limit on the CAD tessellator.

SVG is a headless selection aid. Its orthographic painter ordering can be
approximate for interlocking parts or intersecting display triangles; label
anchors have a separate bounded depth test. A 3D application can consume the
same tagged display arrays for interactive inspection.

## Reproducible, bounded meshing workflows

Recipes capture source identity, companion groups, sizing controls, worker and
resource limits, and an explicit retry policy. They use the same native mesher
and validity checks as `CadModel.mesh()`.

```python
model = pm.read_step("housing.step", retain_source=True)
bore = model.select_faces(surface="cylinder", radius=2.5, expect=1)
model.define_group("bore_wall", bore, group_id="housing-bore-wall")

recipe = pm.create_recipe(
    model, 1.0,
    selected_face_sizes=[(bore, .2)],
    volume_controls=[pm.BoxSize([-3, -3, 0], [3, 3, 5], .3)],
    assembly="conforming", min_size=.15, growth=1.3,
    workers=1, validation_workers=1, max_work=100_000_000,
    allow_advanced=True, finer_factors=(.7,),
    max_attempts=3, max_seconds=90,
)
pm.save_recipe(model, recipe, "housing.recipe.json")

try:
    result = pm.run_recipe(model, recipe)
    print(result.report["attempts"])
    pm.publish_result(result, "housing-run-001")
except pm.WorkflowError as error:
    print(error.stage, error.reason, error.report)
```

The default is one standard attempt. Additional attempts need an explicit
`max_attempts` and an enabled alternative. Advanced recovery is eligible only
for diagnosed seed/cavity/material recovery failures. A finer attempt needs an
eligible surface, CAD error or size/quality failure. Unknown failures, exhausted
budgets, cancellation, allocation failures, unresolved interfaces and unsupported
controls stop the workflow. Assembly policy is preserved; interference is
reported rather than silently fusing parts. Quality warnings can accompany a
valid native result and do not trigger an extra attempt by themselves.

`finer_factors` are strictly decreasing fractions of the original global size.
They change only that global size. Existing face/box controls, growth, minimum
size, assembly policy and resource budgets remain unchanged. A candidate below
`min_size` is skipped and recorded. This is not a claim that every achieved mesh
edge becomes uniformly smaller.

Every native attempt records its exact options, trigger, elapsed time and full
native report, including failures, while the diagnostics fit the 8 MiB report
budget. If that budget is exhausted, the workflow stops with `report_budget`,
records the attempted work and explicitly marks omitted or summarized data;
it does not relabel a successful native mesh as invalid. Reports distinguish
independent validated
parts from a globally conforming domain. A workflow cannot turn an incomplete
or refused native mesh into a successful result.

`run_recipe()` times validation, preparation, queueing and all retries after the
caller has imported the model. For import-inclusive timing, use
`run_saved_recipe("housing.recipe.json", "housing.step")`. The source path must
be supplied explicitly; a recipe's display name is never opened as a path.
The total deadline is cooperative: queue waits and CAD calls can overrun it.
The remaining time is checked before each attempt and passed to the native
operation; any overrun remains visible in the report. Node, cell, work and
memory/output budgets apply per attempt. Publication is separate and timed
separately.

`publish_result()` stages six files in a sibling directory and publishes them
together with an exclusive atomic rename:

- `source.step`: the exact import-time bytes, verified against the source hash;
- `intent.json`: captured group IDs, queries and property guards;
- `recipe.json`: controls and retry policy;
- `mesh.paacsmesh` and `mesh.vtk`: native tagged exports;
- `report.json`: all attempts, native audits and publication metadata.

An existing destination is refused, including an empty directory or symlink.
Export or publication failures leave the owned result usable and previous
outputs intact. Filesystems without exclusive atomic directory rename are
refused. The default embedded-source limit is 64 MiB. Ordinary imports need the
original source to remain unchanged; `retain_source=True` keeps its exact bytes
for later export. Companion groups are captured when the recipe is created, so
later edits to a live model's group dictionary do not change the recipe.

```python
reopened = pm.read_step("housing-run-001/source.step",
                        intent="housing-run-001/intent.json")
saved = pm.load_recipe(reopened, "housing-run-001/recipe.json")
again = pm.run_recipe(reopened, saved)
```

Recipes and intent files reject edited sources or incompatible kernel/tag
schemes. Rebind and visually confirm selections explicitly before creating a
new recipe for a changed design.

Pin the source commit or qualified wheel checksum when comparing runs. Recipes
record inputs and engine/kernel versions; different builds or platforms may
produce different meshes. Every replay reruns the native validity checks.

## Build and qualify a wheel

Windows builds require Visual Studio 2022 C++ build tools and CPython 3.12 x64.
Run from the repository root:

```powershell
python -m venv build/python-dev
build/python-dev/Scripts/python.exe -m pip install build==1.6.1 delvewheel==1.13.2 cmake==4.4.3 ninja==1.13.2
build/python-dev/Scripts/python.exe scripts/build_python_wheel.py --jobs 4
```

On macOS 15 ARM64, install Xcode command-line tools and ARM64 CPython 3.12, then:

```sh
python3.12 -m venv build/python-dev
build/python-dev/bin/python -m pip install build==1.6.1 delocate==0.13.0 cmake==4.4.3 ninja==1.13.2
build/python-dev/bin/python scripts/build_python_wheel.py --jobs 3
```

Linux wheels build inside the digest-pinned manylinux 2.28 x64 container in the
[workflow](../.github/workflows/python-wheel.yml). Its compiler and glibc target
establish the compatibility floor before auditwheel repairs the dependency
closure. The builder refuses a different glibc build host. See the workflow for
the complete container command and image digest.

The script builds the SHA-512-pinned OCCT 8.0.1 source, builds the extension and
a matching native CLI, then repairs the wheel with delvewheel, delocate or
auditwheel. Windows runtime files are checked against their unmodified
redistributable originals. macOS architecture and deployment target propagate
through both OCCT bootstrap layers.

The wheel is installed in a fresh virtual environment with SDK, loader and CAD
variables removed. Qualification exercises STEP import, Tet4 publication,
failures, sizing, tagged export and both supported NumPy series. A separate
check compares complete Python/native CLI archives for a cylinder, an assembly
and combined face/volume sizing at matching settings.

Actions also downloads each artifact onto a separate runner or fresh container
that has never built OCCT, and repeats the installed-package checks. The Linux
consumer uses glibc 2.28 and a different mount path from the builder. This checks
that repaired libraries resolve without absolute paths into the build SDK.

`dist/` receives the repaired wheel, matching OCCT source archive, dependency
audit, native parity receipt, sizing measurements and validation receipt. Keep
the source archive and notices available when redistributing the wheel. Use
this script for distributable wheels; backend output requires platform repair.

For local iteration, reuse a source-built shared OCCT 8.0.1 installation:

```powershell
build/python-dev/Scripts/python.exe scripts/build_python_wheel.py `
  --occt-dir C:/path/to/occt/lib/cmake/opencascade --jobs 4
```

The CMake option `PAACS_BUILD_PYTHON` defaults to `OFF`. Enabling it requires
`PAACS_WITH_OCCT=ON`; `PAACS_PYTHON_BUILD_OCCT=ON` builds the pinned SDK when an
existing SDK is not supplied. Ordinary C++ SDK builds do not require Python or
binding dependencies. The wheel installs only the `Python` component, excluding
the C++ SDK headers and static libraries.

## Next milestones

Extend the qualified platform and Python-version matrix.
Extend the qualified workflow and selection tools for remote automation.
A thin tool adapter can build on the bounded Python functions without replacing
the native engine or its validity checks.
