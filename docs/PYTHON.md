# Python STEP interface

The first Python distribution targets **CPython 3.12 on Windows x64**. It calls
the C++ STEP importer and Tet4 mesher in process. The repaired wheel contains
the shared CAD libraries and runtime dependencies; installation does not need
a compiler, an OCCT SDK or changes to `PATH`.

Wheels are development artifacts from the
[Python STEP wheel workflow](../.github/workflows/python-wheel.yml).
This interface has not been published to PyPI. Download and extract the workflow
artifact, then install its `.whl` file in a Python 3.12 environment:

```powershell
python -m pip install --only-binary=:all: .\paacs_mesh-1.2.0rc1-cp312-cp312-win_amd64.whl
```

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

## Build and qualify a wheel

Building requires Visual Studio 2022 C++ build tools and CPython 3.12 x64. Run
from the repository root:

```powershell
python -m venv build/python-dev
build/python-dev/Scripts/python.exe -m pip install build==1.6.1 delvewheel==1.13.2
build/python-dev/Scripts/python.exe scripts/build_python_wheel.py --jobs 4
```

The script builds the SHA-512-pinned OCCT 8.0.1 source, builds the extension,
repairs the wheel's DLL dependencies and verifies the bundled Microsoft runtime
files. It installs the wheel into a fresh virtual environment and tests STEP
import, Tet4 publication, failures and tagged exports with SDK paths removed
and CAD environment variables cleared. `dist/` receives the repaired wheel,
matching OCCT source archive and a validation receipt. Keep the source archive
and notices available when redistributing the wheel.
Use this script for a distributable Windows wheel: the build backend produces
an intermediate wheel that still needs DLL repair.

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

Extend this interface with combined face/volume sizing and growth, NumPy mesh
arrays, then macOS and Linux wheels. CAD inspection, semantic selections and
an inspect–select–mesh–audit loop build on those foundations. These features
will use the native engine's controls and classification rather than a separate
Python meshing implementation.
