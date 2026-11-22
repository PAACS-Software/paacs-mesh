# PAACS Mesh 1.2.0-rc.2

This candidate adds a Python STEP package and improves initial Tet4 point
placement. The C++ core continues to build independently of OCCT.

- Compatible single-solid Tet4 requests construct wall sites before filling
  interior gaps. C++, Python and the STEP CLI provide explicit on/off controls.
  Multi-part assemblies and specialized routes retain their documented policies.
- Local cavity certification and cached source-contact checks reduce repeated
  geometry work. Native validity, prescribed surfaces and CAD tags remain checked.
  Automatic source-stage fallback records its failed work and elapsed time.
- CPython 3.12 wheels bundle the CAD runtime on Windows x64, macOS 15 ARM64
  and manylinux 2.28 x64. Python supports face and volume sizing, NumPy snapshots,
  geometric selections, bounded recipes, tagged exports and a terminal command.
  The PyPI workflow requires qualified builders and separate clean consumers
  from one source revision before publication.

The native STEP features introduced in the previous candidate remain available:

- Native Tet4 combines a global size, a selected CAD face and one or more
  refinement boxes. Shared CAD edges, trimmed surfaces and volume refinement
  use the same size field. The native checks run before publication.
- A saved recipe binds face tags to the input STEP hash. The C3D4 example
  checks the exported volume and boundary groups against the native archive.
- The optional STEP build downloads pinned OCCT 8.0.1 and builds it as shared
  libraries. CI uses the same build on Windows, Linux and macOS, with STEP CLI
  checks on Windows and macOS and a Windows solver-export check.
- The plate-with-hole and hybrid straight-duct examples document local
  refinement, stress convergence and conforming Prism6/Tet4 wall layers.

Native STEP remains a preview. Single-solid Tet4 sizing is covered by the
fixture tests. Trimmed surfaces, assemblies and volume recovery have bounded
failure routes. The [candidate validation](docs/RELEASE_READINESS.md) records
which cases completed and the settings used.

The CMake package version is `1.2.0`. The SDK and CLI report `1.2.0-rc.2`.
Rebuild downstream C++ applications with the candidate's headers and libraries.
See [release policy](docs/RELEASE_POLICY.md) and [licensing status](LICENSING_STATUS.md).
The mesh archive remains format v1. The [Python guide](docs/PYTHON.md) records
platform limits and installation; [package publication](docs/PYPI_RELEASE.md)
describes the qualified release workflow.

Earlier changes are in the [1.1.0 notes](docs/RELEASE_NOTES_1.1.0.md) and [changelog](CHANGELOG.md).
