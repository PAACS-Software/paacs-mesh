# Release scope

## Core library and SDK

The core includes the primitive, mapped, sweep and constrained meshing routes
in the [README](../README.md). Each route documents its input geometry, element
families and controls. C++26/GCC 16 is the reference configuration. The
C++23 compatibility build is also used by the Windows STEP job.

The `PAACS::mesh` target depends on the C++ standard library and platform
threads. Applications can link it directly. OCCT and the debug UI bridge
are optional targets.

Use `find_package(PAACSMesh 1.2.0 EXACT CONFIG REQUIRED)` for this package.
The SDK and CLI report `1.2.0-rc.1`. Rebuild C++ consumers with the matching
headers, libraries and toolchain. The mesh archive remains format v1.
The `paacs_ui_` C bridge keeps ABI version 1. The Python and JSON examples
cover the routes listed in [UI integration](UI_INTEGRATION.md).

## Native STEP preview

The optional OCCT adapter handles supported STEP and BRep models. It retains
native CAD face and body tags and validates geometry before publishing an
archive. Joint face and box sizing is available for a single Tet4 solid.
The [native CAD guide](NATIVE_CAD.md) describes recovery options and assembly
policies. The [candidate validation](RELEASE_READINESS.md) records measured
coverage and open cases.

Straight-sided Tet10, mapped or swept Hex8 and prepared Prism6 layers have
their own documented scope. Broader curved Tet10, native CAD wall layers and
Hex27 remain on the [roadmap](ROADMAP.md).

## Source package

```sh
python3 scripts/verify_version.py
python3 scripts/package_source.py --output out/paacs-mesh-1.2.0-rc.1.zip
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
