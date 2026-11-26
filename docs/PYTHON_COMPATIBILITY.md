# Python compatibility

Version **1.2.0** targets the expanded matrix below: regular,
64-bit CPython 3.12–3.15 on five native platforms. Wheels bundle the CAD runtime and
install without a compiler or a separate OCCT SDK. Package metadata restricts
Python to `>=3.12,<3.16`; each interpreter still needs its matching binary wheel.

Stable 1.2.0 requires its own complete qualification. Its matching GitHub release carries
the publication receipts. The expanded distribution implementation passed the
historical qualification recorded below; that evidence does not qualify the
current versioned release.
Official releases pass the release gates from their own immutable tagged source
before publication. Published release-candidate files and tags remain immutable. The
earlier `1.2.0rc2` release remains limited to
CPython 3.12 on Windows x64, macOS ARM64 and Linux x64; its files are unchanged.

## Qualification evidence

[PR #11](https://github.com/PAACS-Software/paacs-mesh/pull/11) qualified the expanded
implementation at commit `44a22fe00333191e2900d8053958259198d8ca4f`:

- [Python qualification](https://github.com/PAACS-Software/paacs-mesh/actions/runs/45095800374)
  passed all **55 jobs**: 20 wheel builders, 28 independent consumers, five CAD SDK
  jobs, matrix generation and release-assembly regression checks.
- [Native qualification](https://github.com/PAACS-Software/paacs-mesh/actions/runs/45095800447)
  passed all **seven jobs**, including GCC thread and address/undefined
  sanitizers and native STEP checks on Windows, Linux and macOS.
- A local release-assembly dry run verified all 48 downloaded builder/consumer
  artifacts against their GitHub ZIP digests, assembled the 20 wheels and passed
  strict Twine metadata checks for every wheel.

The qualified workflow checkout was clean virtual merge
`3d79164059bfe274d20c810e8c8255296c5fbc7c`, with source tree
`eca306ca6e2d9b623f9f58af48832749c3a7187f`. That tree matched the reviewed candidate.
These implementation artifacts used the earlier development version; they are
evidence for the expanded matrix, not published RC3, RC4 or stable files. RC3
and RC4 subsequently repeated qualification from their own versioned sources.
These historical runs do not establish a stable 1.2.0 result.

### RC4 publication evidence

RC4 used source commit `7aa215890268a64d63760b2497d82ffd7544ba29`, whose tree
`670c0df9976d3fbe20965aab120308973511d17c` matched the pre-tag qualified tree.

- [Native qualification](https://github.com/PAACS-Software/paacs-mesh/actions/runs/45752886648)
  passed all seven jobs.
- [Python qualification](https://github.com/PAACS-Software/paacs-mesh/actions/runs/45752886693)
  passed all 55 jobs.
- [Tagged-source publication](https://github.com/PAACS-Software/paacs-mesh/actions/runs/45762427424)
  passed all 60 jobs, including a fresh complete wheel/consumer matrix, assembly,
  source assets, Trusted Publishing and the package-index installation check.

The [RC4 release](https://github.com/PAACS-Software/paacs-mesh/releases/tag/v1.2.0-rc.4)
retains source and qualification receipts, including the independent Windows
PyPI installation record. These receipts qualify RC4 only; stable 1.2.0 must
repeat the required checks with its own versioned artifacts.

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

This matrix excludes Intel macOS, Python 2, 32-bit interpreters, PyPy, debug builds
and free-threaded CPython ABIs. Native C++ builds remain independent of the Python
wheel matrix and continue through their existing regression checks.

## One source of configuration

[python/compatibility.json](../python/compatibility.json) defines the interpreter
versions, native platforms, pinned Linux images, NumPy checks and consumer
environments. [scripts/python_compatibility.py](../scripts/python_compatibility.py)
validates that matrix and supplies the build, installation and release gates. An
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

Windows ARM64 builds pin the MSVC **14.44 toolset** (compiler 19.44) for both the
shared OCCT SDK and the Python bindings. The runner's default MSVC 14.52 toolset failed
while compiling OCCT, with a `UtcBldCompileSerializedFile` stack followed by D8040.
An [upstream runner issue](https://github.com/actions/runner-images/issues/15328)
reports the same compiler version and stack; selecting 14.44 compiled the same
source on the same image, and the runner maintainer recommends that workaround
while investigating.

MSVC 14.52 is a stable toolset release; this pin addresses an observed ARM64 compiler
regression. Both build stages select the same installed toolset and record the
resolved compiler in their SDK context and qualification receipts. The changed
ARM64 compiler produces a distinct SDK cache fingerprint. Other platforms retain
their existing toolchain selection.

## Publication gate

Stable publication requires all 20 wheels and 28 independent consumers to pass from one
clean immutable source revision, along with the native regression checks and
release assembly checks. Before tagging, record all seven native and all 55 Python
jobs and verify their tested tree matches merged main and the tag target. Assembly
also requires named v2 native/Python resource-policy parity evidence and verifies
the exact target set, interpreter and platform tags, metadata, source identity,
dependency receipts and artifact hashes.
A missing or failed combination blocks publication, including a failed macOS 27
preview check. Existing release tags and PyPI files remain immutable.

The [release guide](PYPI_RELEASE.md) describes staging, Trusted Publishing and
verification of the published files. A passing implementation PR does not replace
qualification of the versioned release source or establish package-index
availability. Published versions and their checksums identify the files offered
through the package index; development artifacts retain their own source and
qualification receipts.
