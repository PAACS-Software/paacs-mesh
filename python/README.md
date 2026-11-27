# PAACS Mesh

PAACS Mesh imports STEP geometry, builds audited Tet4 volume meshes and exports
tagged VTK or native mesh archives. The Python package calls the C++ engine in
process and bundles its CAD runtime.

Version **1.2.1** targets regular, **64-bit CPython 3.12–3.15**
on Windows x64 and Windows 11 ARM64, Linux x64/ARM64 with glibc 2.28 or newer,
and macOS ARM64 with a minimum of macOS 15. The macOS wheels must pass
installation on macOS 15, 26 and 27.
Wheel installations need no C++ compiler, separate OCCT installation or
loader-path configuration. NumPy is optional. Intel macOS, Python 2, PyPy, debug
interpreters and free-threaded Python ABIs are excluded.

The expanded implementation passed 20 wheel builds, 28 independent installations
and the native regression jobs in
[PR #11](https://github.com/PAACS-Software/paacs-mesh/pull/11).
The [compatibility guide](https://github.com/PAACS-Software/paacs-mesh/blob/main/docs/PYTHON_COMPATIBILITY.md)
records that evidence and the required release checks. Official releases are
built and qualified from their own immutable tagged source before publication.
Version 1.2.1 requires fresh qualification from its own source; the matching
GitHub release carries the publication receipts. Historical PR or earlier
release results do not substitute for the versioned release checks.

Install the pinned release:

```sh
python -m pip install paacs-mesh==1.2.1
```

Published 1.2.0 and release-candidate files and tags remain immutable. The earlier
`1.2.0rc2` release supports CPython 3.12 on Windows x64,
macOS ARM64 and Linux x64 only; its files remain unchanged. Use
`python -m pip install --pre paacs-mesh` to allow prereleases when selecting a version.
Stable releases use `python -m pip install paacs-mesh`.

```python
import paacs_mesh as pm

model = pm.read_step("part.step")
mesh = model.mesh(size=0.7)
mesh.write("part.vtk")
print(mesh.report["quality"])
```

Sizes and transferred STEP coordinates use millimetres. The report records
native validity checks, quality, timings, CAD tags and assembly interference.
Parts are preserved by default; interference is reported rather than silently
fusing overlapping material.

The same operation is available from a terminal:

```sh
paacs-mesh part.step part.vtk --size 0.7
```

The package includes bounded mesh quality inspection, SVG views, comparisons, Normal
resource admission and default-off experimental thin-region ordering. The
[release notes](https://github.com/PAACS-Software/paacs-mesh/blob/main/RELEASE_NOTES.md)
record their scope and the required qualification. Capability discovery and
CAD inspection are available before selecting controls:

```sh
paacs-mesh capabilities --json
paacs-mesh inspect part.step --topology --json
```

Version 1.2.1 adds optional explicit boundary, point-generation and volume stages:

```python
from paacs_mesh import pipeline

boundary = pipeline.prepare_boundary(model, size=0.7)
boundary.save("saved-boundary")
boundary.surface().write("surface.vtk")
boundary.edges().save("cad-edges")
cloud = boundary.generate_initial_points()
cloud.save("initial-points")
mesh = boundary.mesh(cloud)
mesh.write("replayed.vtk")
```

Use new output directories. The initial scope is one original solid with uniform
sizing, native Tet4 construction and the documented recovery route. Surface
artifacts retain the qualified Tri3 boundary; edge artifacts contain retained
ordered CAD-edge chains, with ambiguous extraction reported explicitly. Display
previews are separate from qualified meshing artifacts.

Saved boundaries and clouds carry exact source, boundary, build and sizing
identities. Generating points stops before constraint recovery and volume fill,
but includes the classical placement's auxiliary tetrahedral work. A cloud is
not a valid volume mesh; final construction retains native publication checks.
No learned model or training dependency is included. See the
[pipeline guide](https://github.com/PAACS-Software/paacs-mesh/blob/main/docs/PIPELINE.md)
for resource limits, saved-cloud replay and NumPy copies. Ordinary `model.mesh()`
algorithms and defaults are unchanged.

Saved recipes, readable diagnostics and the complete inspect/select/mesh/compare
example are covered in the
[inspection guide](https://github.com/PAACS-Software/paacs-mesh/blob/main/docs/AGENT_WORKFLOW.md).

Face refinement, volume regions, growth controls, optional NumPy arrays,
geometric selections and bounded saved recipes are described in the
[Python guide](https://github.com/PAACS-Software/paacs-mesh/blob/main/docs/PYTHON.md).
Coverage and geometry limitations are recorded in the
[release qualification](https://github.com/PAACS-Software/paacs-mesh/blob/main/docs/RELEASE_READINESS.md).

PAACS Mesh is licensed under MPL-2.0 by PAACS Software Ltd. The separately bundled
OCCT libraries use LGPL-2.1 with the Open CASCADE exception. Notices are included
in the wheel; the matching verified CAD source and build recipe accompany each
[release](https://github.com/PAACS-Software/paacs-mesh/releases).
