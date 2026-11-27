# PAACS Mesh 1.2.1

Version 1.2.1 adds optional Python access to the CAD boundary, initial-point and
volume stages. Applications can save a checked boundary, export its surface and
CAD edges, obtain classical interior points before volume construction, and
replay prescribed points through the native mesher.

Ordinary `CadModel.mesh()` algorithms and defaults are unchanged. No learned
model, training code or machine-learning dependency is included. Published
1.2.0 remains an immutable reference; its tags, wheels and qualification records
are not replaced or relabelled.

## Installation and compatibility

```sh
python -m pip install paacs-mesh==1.2.1
```

The source, SDK, CLI and Python distribution report `1.2.1`. The wheel matrix
remains regular 64-bit CPython 3.12–3.15 on Windows x64 and Windows 11 ARM64,
Linux x64/ARM64 with glibc 2.28 or newer, and macOS ARM64 with minimum macOS 15.
The same macOS wheels must pass installation on macOS 15, 26 and 27. CAD runtime
dependencies are bundled; NumPy remains optional. Intel macOS, PyPy, 32-bit,
debug and free-threaded Python ABIs remain outside this matrix.

## Explicit meshing stages

```python
import paacs_mesh as pm
from paacs_mesh import pipeline

model = pm.read_step("part.step", retain_source=True)
boundary = pipeline.prepare_boundary(model, size=0.7)
boundary.save("saved-boundary")
boundary.surface().write("surface.vtk")
boundary.edges().save("cad-edges")

cloud = boundary.generate_initial_points()
cloud.save("initial-points")
mesh = boundary.mesh(cloud)
mesh.write("part.vtk")
```

Use fresh output paths. The first staged scope is one original solid, uniform
sizing and native conforming Tet4 construction. Unsupported controls refuse.
Assemblies and local face or box sizing retain their ordinary meshing
interfaces; they are not admitted by this initial pipeline contract.

Boundary preparation completes native CAD, winding, geometric-error and exact
embedded-shell checks before interior placement. Saving first preserves that
surface even if a later volume attempt fails. Loading reimports the exact STEP
source and repeats native admission; hashes identify bytes rather than replacing
geometry checks.

`boundary.surface()` exposes the qualified, tagged Tri3 complex.
`boundary.edges()` exposes ordered chains on actual CAD edges, with edge tags
and CAD topology metadata. It does not substitute all triangle edges or infer
ownership from proximity. Ambiguous edge extraction can refuse. Display previews
remain separate from qualified meshing artifacts.

Point generation includes the classical placement's auxiliary tetrahedral
construction and topology finalization, then stops before constraint recovery,
material fill and volume validation. Its immutable cloud preserves point order
and exact coordinates. A generated empty cloud differs from an unavailable
capture. Neither a point cloud nor a prepared surface certifies a volume mesh.

Saved clouds bind to the source, frozen boundary, build and sizing contract.
`pipeline.load_initial_points()` and `pipeline.load_boundary()` support replay
with that same pinned build. Prescribed points are checked all-or-nothing;
invalid, duplicate or unretained points cannot silently fall back to a new cloud.
Final construction preserves the tagged, oriented boundary and native CAD,
volume, audit and publication checks. Recovery may insert additional points,
which remain visible in the final counts.

See the [pipeline guide](docs/PIPELINE.md) for save/reopen, exports, optional NumPy
copies, provenance and resource controls. The earlier development spelling
`paacs_mesh.experimental` remains available for compatibility.

## Scope and release checks

Work, modeled memory, count, output, cancellation and cooperative time limits
remain authoritative. Modeled memory does not bound all process RSS, and an
operation deadline is not a global retry budget. Complete-workflow comparisons
must count import, preparation, point generation, replay, validation, exports
and all attempted retries; named stage spans are not an additive decomposition.

Native STEP meshing remains a preview with unsupported geometries and measured
quality limitations. Package metadata remains Beta. This release adds no claim
of arbitrary-CAD completion, solver accuracy or improved learned meshing. The
[historical 51-file results](docs/RELEASE_READINESS.md),
[owned hollow study](docs/THIN_HOLLOW_QUALIFICATION.md) and
[1.2.0 qualification protocol](docs/STABLE_1.2.0_QUALIFICATION.md) retain their
original versions, settings and evidence.

The final 1.2.1 tree must pass native regression jobs, the staged native
checkpoint regression and the complete 20-wheel/28-consumer qualification before
publication. Installed tests cover staged save/reopen, surface and CAD-edge
exports, point identity, supplied replay and transactional failures alongside
existing interfaces. Historical checks and private development wheels do not
qualify new release artifacts.

The matching GitHub release carries source, dependency and qualification
receipts. Publication uses immutable tagged source and a fresh hash-pinned PyPI
installation check. See [Python compatibility](docs/PYTHON_COMPATIBILITY.md)
and [publication gates](docs/PYPI_RELEASE.md).

Mesh archive format v1, C bridge ABI v1 and the optional nature of OCCT remain
unchanged. Rebuild C++ consumers with matching headers, libraries and toolchain.
The project remains MPL 2.0, copyright PAACS Software Ltd.

The [1.2.0 notes](https://github.com/PAACS-Software/paacs-mesh/blob/v1.2.0/RELEASE_NOTES.md)
and [changelog](CHANGELOG.md) retain earlier release history.
