# PAACS Mesh

PAACS Mesh imports STEP geometry, builds audited Tet4 volume meshes and exports
tagged VTK or native mesh archives. The Python package calls the C++ engine in
process and bundles its CAD runtime.

The **1.2.0rc3 release candidate** supports regular, **64-bit CPython 3.12–3.15**
on Windows x64 and Windows 11 ARM64, Linux x64/ARM64 with glibc 2.28 or newer, and macOS ARM64
with a minimum of macOS 15. The macOS wheels are tested on macOS 15, 26 and 27.
Wheel installations need no C++ compiler, separate OCCT installation or
loader-path configuration. NumPy is optional. Intel macOS, Python 2, PyPy, debug
interpreters and free-threaded Python ABIs are excluded.

The expanded implementation passed 20 wheel builds, 28 independent installations
and the native regression jobs in
[PR #11](https://github.com/PAACS-Software/paacs-mesh/pull/11).
The [compatibility guide](https://github.com/PAACS-Software/paacs-mesh/blob/main/docs/PYTHON_COMPATIBILITY.md)
records that evidence and the required release checks. Official releases are
built and qualified from their own immutable tagged source before publication.

Install the pinned release candidate:

```sh
python -m pip install paacs-mesh==1.2.0rc3
```

The previous published `1.2.0rc2` release supports CPython 3.12 on Windows x64,
macOS ARM64 and Linux x64 only; its files remain unchanged. Use
`python -m pip install --pre paacs-mesh` to select the newest available candidate.
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

RC3 and the qualified development wheels provide discovery and CAD inspection
before selecting controls:

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
