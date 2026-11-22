# PAACS Mesh

PAACS Mesh imports STEP geometry, builds audited Tet4 volume meshes and exports
tagged VTK or native mesh archives. The Python package calls the C++ engine in
process and bundles its CAD runtime.

The qualified wheels support **64-bit CPython 3.12** on Windows x64, macOS 15 or
newer on ARM64, and Linux x64 with glibc 2.28 or newer. These wheels need no C++
compiler, separate OCCT installation or loader-path configuration. NumPy is
optional.

```sh
python -m pip install paacs-mesh
```

For release candidates, explicitly include prereleases:

```sh
python -m pip install --pre paacs-mesh
```

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

Discover the available operations and inspect CAD faces before selecting controls:

```sh
paacs-mesh capabilities --json
paacs-mesh inspect part.step --topology --json
```

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
