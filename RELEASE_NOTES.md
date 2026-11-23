# PAACS Mesh 1.2.0-rc.3

This candidate expands the Python distribution and adds inspection and diagnostic
commands for engineers, scripts and agents. It preserves the native meshing
algorithms and validity checks from the previous candidate.

## Python installation

```sh
python -m pip install --pre --upgrade paacs-mesh
```

For a reproducible installation, use `paacs-mesh==1.2.0rc3`. The source version is
`1.2.0-rc.3`; Python package tools normalize it to `1.2.0rc3`.

The release requires a wheel for every regular 64-bit CPython version listed:

| Platform | Python | System requirement |
|---|---|---|
| Windows x64 | 3.12, 3.13, 3.14, 3.15 | Windows x64; independently tested on Windows Server 2022 |
| Windows ARM64 | 3.12, 3.13, 3.14, 3.15 | Windows 11 ARM64 |
| Linux x64 | 3.12, 3.13, 3.14, 3.15 | glibc 2.28 or newer |
| Linux ARM64 | 3.12, 3.13, 3.14, 3.15 | glibc 2.28 or newer |
| macOS ARM64 | 3.12, 3.13, 3.14, 3.15 | macOS 15 or newer; same wheels tested on 15, 26 and 27 |

The CAD runtime is bundled: wheel installations need no compiler or separate
OCCT installation. NumPy remains optional. Intel macOS, Python 2, PyPy, 32-bit,
debug and free-threaded Python builds are outside this wheel matrix. The macOS 27
CI runner is a public preview; its actual OS version and result are checked.

## Inspection and diagnostics

The terminal now exposes `capabilities`, paginated `inspect`, saved-recipe `run`
and `diagnose`. Each supports readable output and `--json`; machine output stays
on stdout while native messages go to stderr. Unicode CAD names round-trip even
when the terminal cannot represent their characters directly.

Python adds `capabilities()`, `diagnose()` and `format_diagnostics()`. Diagnostics
report the native checks already completed, bounded body evidence, stable issue
codes and suggested actions. Missing measurements are explicit. They do not
invent per-cell or per-face quality localization or independently recertify a
supplied report. Assembly body identities and interference reporting are retained.

See the [inspection workflow](docs/AGENT_WORKFLOW.md) and its owned cylinder example.
Face sizing, refinement boxes, growth controls, NumPy snapshots, guarded selections,
source-bound recipes and tagged exports remain available.

## Build and release validation

The complete release gate requires **20 repaired wheels and 28 independent clean
installation results** from one immutable source revision. Each target exercises
STEP import, meshing, auditing, controls, tags, exports, deterministic repeats and
the CLI, followed by the NumPy versions assigned to that target. Native C++/Python
parity, binary architecture, compiler configuration and dependency hashes are checked.

OCCT SDK reuse checks the source pin, build recipe, toolchain and installed files.
Windows ARM64 selects the validated MSVC 14.44 toolset. A Linux host-directory fix
allows verified SDK caches to be restored before container builds.

Before publication, the native regression jobs and sanitizers must also pass.
The GitHub release carries the matching OCCT source, build and qualification
receipts, source archive and checksums. A final index check verifies the published
wheel hashes and performs a fresh pip installation, import, mesh and export.
See [compatibility](docs/PYTHON_COMPATIBILITY.md) and [publication](docs/PYPI_RELEASE.md).

Native STEP remains a preview with bounded failure routes for unsupported geometry
and recovery limits. This release changes the interface and packaging rather than
the meshing algorithms; the [recorded CAD corpus results](docs/RELEASE_READINESS.md)
are the existing engine baseline, not a new corpus run.

The CMake package version remains `1.2.0`; SDK and CLI report `1.2.0-rc.3`.
Rebuild C++ consumers with matching headers and libraries. Mesh archive format v1
and C bridge ABI version 1 are unchanged. The C++ core builds independently of OCCT.
The project remains MPL 2.0, copyright PAACS Software Ltd.

Earlier changes are in the [rc.2 notes](docs/RELEASE_NOTES_1.2.0-rc.2.md) and
[changelog](CHANGELOG.md).
