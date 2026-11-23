# Python compatibility

The published `1.2.0rc2` release supports regular, 64-bit CPython 3.12 on Windows
x64, macOS ARM64 with a minimum of macOS 15, and Linux x64 with glibc 2.28 or newer.
It bundles the CAD runtime and installs without a compiler or a separate OCCT SDK.

The following expanded matrix is a **development qualification plan**. Its new
interpreter and architecture combinations are not yet qualified or available on
PyPI. Development metadata accepts `>=3.12,<3.16` so the candidate wheels can be
built and tested; metadata alone does not establish binary compatibility.

## Required artifacts

Each row requires a wheel for every regular CPython version listed. Every wheel
must pass its builder tests and a separate clean installation on the specified
consumer environments.

| Platform | Python versions | Wheel platform tag | Independent consumer environments |
|---|---|---|---|
| Windows x64 | 3.12, 3.13, 3.14, 3.15 | `win_amd64` | Windows Server 2022 |
| macOS ARM64 | 3.12, 3.13, 3.14, 3.15 | `macosx_15_0_arm64` | macOS 15, 26 and 27 |
| Linux x64 | 3.12, 3.13, 3.14, 3.15 | `manylinux_2_28_x86_64` | Native x64, glibc 2.28 container |
| Linux ARM64 | 3.12, 3.13, 3.14, 3.15 | `manylinux_2_28_aarch64` | Native ARM64, glibc 2.28 container |
| Windows ARM64 | 3.12, 3.13, 3.14, 3.15 | `win_arm64` | Windows 11 ARM64 |

The complete set contains **20 wheels and 28 independent consumer results**.
The extra eight consumer results install the existing macOS 15 wheels on macOS 26
and 27; they do not rebuild those wheels. macOS 15 remains the deployment target.
Linux repair may also produce an audited additional manylinux tag; the required
2.28 tag and matching architecture must remain present.

macOS 27 currently uses GitHub's `xcode-27` **public-preview** runner. Qualification
checks the actual operating system major version, not just the workflow label.
Its checks remain required. See the [runner announcement](https://github.blog/changelog/2026-09-10-xcode-27-runner-image-now-runs-on-macos-27/)
and [installed image details](https://github.com/actions/runner-images/blob/main/images/macos/xcode-27-arm64-Readme.md).

This plan excludes Intel macOS, Python 2, 32-bit interpreters, PyPy, debug builds
and free-threaded CPython ABIs. Native C++ builds remain independent of the Python
wheel matrix and continue through their existing regression checks.

## One source of configuration

[python/compatibility.json](../python/compatibility.json) defines the interpreter
versions, native platforms, pinned Linux images, NumPy checks and consumer
environments. [scripts/python_compatibility.py](../scripts/python_compatibility.py)
validates that plan and supplies the build, installation and release gates. An
omitted, duplicated or extra target is an error.

The generated matrices can be inspected locally:

```sh
python scripts/python_compatibility.py --targets
python scripts/python_compatibility.py --consumers
python scripts/python_compatibility.py --platforms
```

Python extension modules use their matching CPython ABI, such as
`cp314-cp314`. A wheel for one minor version is not reused for another minor
version. Builds run on native architectures, with Linux builds inside the
appropriate digest-pinned manylinux 2.28 image. The official
[GitHub runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners),
[interpreter manifest](https://github.com/actions/python-versions/blob/main/versions-manifest.json)
and [manylinux project](https://github.com/pypa/manylinux) document the upstream
runner and interpreter availability. Availability there still requires a passing
PAACS qualification result here.

## Installation checks

Every clean consumer installs the exact repaired wheel identified by the builder
hash, with no SDK paths, compiler configuration or separate CAD installation. The
checks exercise STEP import, Tet4 meshing and native validity checks, sizing
controls, deterministic repeatability, face and body tags, exports, saved
workflows and the installed command. Builder receipts also include native C++/Python
parity and dependency inspection.

The base installation runs without NumPy. Optional array tests then require
NumPy 2.5.5 on every target. CPython 3.12 also tests NumPy 1.26.4 where that release
provides a native wheel: Windows x64, macOS ARM64 and both Linux architectures.
Windows ARM64 and Python 3.13–3.15 use the current NumPy test only. NumPy's
[1.26.4 release notes](https://numpy.org/doc/stable/release/1.26.4-notes.html)
limit that release to Python 3.9–3.12; its absence on another target is not a test
to silently skip. The expected test set is explicit in the compatibility plan.

## Reusing the CAD SDK

The shared OCCT SDK contains no Python extension, so one verified SDK per native
platform and toolchain can serve all four interpreter versions. Its cache
identity includes the pinned upstream source, build recipes, compiler/toolset,
operating system or container, architecture, Release/shared configuration and
macOS deployment target. Restored files and source archives are checked against
their receipts before use. ARM64 and x64 SDKs are distinct artifacts.

This reuse reduces repeated CAD compilation. Wheel repair, interpreter-specific
builds, dependency checks and independent installations still run for every
required combination. Matching OCCT source and redistribution notices remain
part of the release process.

### Windows ARM64 compiler

The Windows ARM64 qualification plan pins MSVC **14.44** for both the shared OCCT
SDK and the Python bindings. The runner's default MSVC 14.52 compiler crashed
while compiling OCCT, with a `UtcBldCompileSerializedFile` stack followed by D8040.
An [upstream runner issue](https://github.com/actions/runner-images/issues/15328)
reports the same compiler version and stack; selecting 14.44 compiled the same
source on the same image, and the runner maintainer recommends that workaround
while investigating.

MSVC 14.52 is a stable release; this pin addresses an observed ARM64 compiler
regression. Both build stages select the same installed toolset and record the
resolved compiler in their SDK context and qualification receipts. The changed
ARM64 compiler produces a distinct SDK cache fingerprint. Other platforms retain
their existing toolchain selection.

## Publication gate

No new package-index release is published as part of this expansion. Before a
future release, all 20 wheels and 28 independent consumers must pass from one
clean immutable source revision, along with the native regression checks and
release assembly checks. Assembly verifies the exact target set, interpreter and
platform tags, metadata, source identity, dependency receipts and artifact hashes.
A missing or failed combination blocks publication, including a failed macOS 27
preview check. Existing release tags and PyPI files remain immutable.

The [release guide](PYPI_RELEASE.md) describes staging, Trusted Publishing and
verification of the published files. Support documentation should be updated to
qualified status only after the corresponding evidence exists.
