# Release scope

Version `1.2.1` uses the [publication gates](PYPI_RELEASE.md). Its
matching GitHub release carries the qualification and publication receipts.
Existing release tags and published files remain immutable.

The package uses stable version numbering for its documented interfaces and
supported routes. Native STEP meshing retains its preview scope, and the Python
package is classified as Beta. Unsupported CAD, quality tails and unqualified
solver workflows remain explicit limitations. This release does not promise
completion for arbitrary geometry.

## Core library and SDK

The core includes the primitive, mapped, sweep and constrained meshing routes
in the [README](../README.md). Each route documents its input geometry, element
families and controls. C++26/GCC 16 is the reference configuration. The
C++23 compatibility build is also used by the Windows STEP job.

The `PAACS::mesh` target depends on the C++ standard library and platform
threads. Applications can link it directly. OCCT and the debug UI bridge
are optional targets.

Use `find_package(PAACSMesh 1.2.1 EXACT CONFIG REQUIRED)` for this package.
The SDK and CLI report `1.2.1`. Rebuild C++ consumers with the matching
headers, libraries and toolchain. The mesh archive remains format v1.
The `paacs_ui_` C bridge keeps ABI version 1. The Python and JSON examples
cover the routes listed in [UI integration](UI_INTEGRATION.md).

## Native STEP preview

The optional OCCT adapter handles supported STEP and BRep models. It retains
native CAD face and body tags and validates geometry before publishing an
archive. Joint face and box sizing is available for a single Tet4 solid.
The [native CAD guide](NATIVE_CAD.md) describes recovery options and assembly
policies. The [historical validation](RELEASE_READINESS.md) records measured
coverage and open cases.

Straight-sided Tet10, mapped or swept Hex8 and prepared Prism6 layers have
their own documented scope. Broader curved Tet10, native CAD wall layers and
Hex27 remain on the [roadmap](ROADMAP.md).

## Optional pipeline access

Version 1.2.1 exposes [staged Python access](PIPELINE.md) without changing
ordinary meshing algorithms or defaults. The initial scope is one imported
original solid, uniform sizing and native conforming Tet4 construction. Local
face or box sizing, assemblies and unsupported construction controls refuse in
this staged API; their existing ordinary meshing interfaces retain their scope.

A saved boundary has native CAD, winding, geometric-error and embedded-shell
checks. Its surface export retains that qualified Tri3 complex; CAD-edge exports
retain actual ordered CAD-edge chains and may refuse ambiguous extraction.
Display previews do not establish those meshing checks. Initial-point generation
uses the classical placement, including auxiliary tetrahedral topology work,
then stops before constraint recovery and material fill. Neither artifact
certifies a volume mesh or simulation accuracy.

Saved clouds bind their ordered coordinates to exact source, boundary, build
and sizing identities. Prescribed points are admitted explicitly; invalid,
duplicate or unretained points cannot silently fall back to another cloud.
Final volume construction keeps native validity gates and explicit resources.
The engine contains no learned model or training dependency. External studies
can pin its published wheel and keep learning code in a separate project.
Published 1.2.0 and earlier tags, wheels and reports remain unchanged.

## Source package

```sh
python3 scripts/verify_version.py --final
python3 scripts/package_source.py --preserve-source-manifest --output out/paacs-mesh-1.2.1.zip
```

The source package includes build definitions, tests, shareable fixtures and
examples. It carries a file manifest, a ZIP checksum, the [MPL 2.0 license](../LICENSE)
and the [PAACS Software Ltd copyright notice](../NOTICE). The
[licensing notes](../LICENSING_STATUS.md) record the terms for dependencies
and dataset material.

Text files use LF line endings in Git and in the source ZIP. The packager
normalizes Windows text line endings before calculating the file manifest;
PNG figures retain their original bytes.

Historical validation pages record paths to private development logs and
figures under `evidence/` and `out/`. Those records are retained separately
from the source package.
