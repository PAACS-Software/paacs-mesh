# PAACS Mesh 1.2.0-rc.1

This candidate adds a documented native STEP workflow alongside the existing
C++ mesh generators and SDK. The core builds independently of OCCT.

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

The CMake package version is `1.2.0`. The SDK and CLI report `1.2.0-rc.1`.
Rebuild downstream C++ applications with the candidate's headers and libraries.
See [release policy](docs/RELEASE_POLICY.md) and [licensing status](LICENSING_STATUS.md).

Earlier changes are in the [1.1.0 notes](docs/RELEASE_NOTES_1.1.0.md) and [changelog](CHANGELOG.md).
