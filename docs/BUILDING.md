# Build and qualification

Use the root CMakeLists.txt. `standalone/` is a forwarding entry point; old
versioned options/package names have been superseded by `PAACS_*`/`PAACSMesh`.
The source list is explicit in `cmake/Sources.cmake`.

## Reference compiler installation on Linux x86-64

The reference compiler was obtained through conda-forge with micromamba:

```sh
micromamba create -y -n paacs-tools -c conda-forge \
  python=3.12 pip gxx_linux-64=16.2.0 libboost-headers=1.92.0
micromamba activate paacs-tools
python -m pip install cmake==4.4.3 ninja==1.13.2
cmake --preset gcc16-release -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
cmake --build --preset gcc16-release --parallel 4
ctest --preset gcc16-release
```

An equivalent system GCC installation is supported. Keep the GCC 16 runtime
with programs compiled against its standard library. `docs/TOOLCHAIN_LOCK.json`
records exact installed reference dependency versions; it is an inventory,
not a complete content-addressed environment lockfile.

## Options

| Option | Default | Meaning |
|---|---|---|
| PAACS_CXX_STANDARD | 26 | 26 or explicit compatibility level 23 |
| PAACS_REQUIRE_GCC16 | OFF | Reference preset sets ON; fails for other compilers |
| PAACS_BUILD_TESTS | ON at root | Full suite, probes, fault injection; needs Boost headers |
| PAACS_BUILD_CLI | ON | Build/install `paacs-mesh` |
| PAACS_WITH_OCCT | OFF | Optional OCCT 8+ planar STEP/BREP importer, CAD CLI and round-trip tests |
| PAACS_BUILD_BENCHMARKS | OFF | Delaunay, Hex, prepared workflows, archive microbenchmark |
| PAACS_BUILD_EXTENDED_PROBES | OFF | Retained compaction, cavity, BND and surface allocation/rollback campaigns |
| PAACS_WARNINGS_AS_ERRORS | OFF | Project library warnings are build errors |
| PAACS_ENABLE_IPO | OFF | Probe support, then enable LTO/IPO |
| PAACS_NATIVE_ARCH | OFF | Host-specific optimization; do not distribute as portable binaries |
| PAACS_SANITIZER | none | none, address-undefined, thread |
| PAACS_UNIT_TEST_TIMEOUT | 600, or 7200 with a sanitizer | Whole inherited suite timeout, in seconds |

Preset `gcc16-ipo` enables IPO. Presets `asan` and `tsan` use Debug and the
corresponding sanitizer; select one per build directory. Pass compiler/path
overrides to each configure preset if the tools are not named gcc-16/g++-16.

```sh
cmake --preset asan
cmake --build --preset asan --parallel 2
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1 ctest --preset asan
```

Do not enable fast-math, unsafe reassociation, or floating-point contraction
for this kernel. The build applies `-fno-fast-math -ffp-contract=off` (or
`/fp:strict`) to the library and public inline code. This is required for the
filtered/exact predicate design; no numerical kernel was converted to SIMD
without a correctness/performance study. C++26 modules, reflection, and
contracts are not prerequisites.

`compile_commands.json` records translation-unit commands. A per-configuration
`build-info-*.txt` records high-level settings. Consumer metadata is separate
from the library's reference build metadata.

## Tests and install

`paacs-tests` links a separate test-instrumented library, keeping fault-injection
hooks out of the installed production library. It accepts an optional
substring filter, e.g. `paacs-tests predicates`. The hardening binary exercises
archive mutations, integer limits, nested parallel calls, spatial queries, and
invalid C requests. The allocation-failure binary overrides global new only in
its isolated process; production allocation behavior is unchanged.

```sh
cmake --install build/gcc16-release --prefix out/install
cmake -S examples/consumer -B build/consumer -G Ninja \
  -DCMAKE_C_COMPILER=gcc-16 -DCMAKE_CXX_COMPILER=g++-16 \
  -DCMAKE_PREFIX_PATH="$PWD/out/install"
cmake --build build/consumer
ctest --test-dir build/consumer --output-on-failure
```

The SDK installs `LICENSE`, `NOTICE`, `LICENSING_STATUS.md` and `CITATION.cff`
under `share/doc/PAACSMesh` by default. The directory is namespaced to PAACS
when another CMake project embeds the library.

Installed core-only consumers use `find_package(PAACSMesh 1.2.0 EXACT CONFIG
REQUIRED)` and link `PAACS::mesh` (the earlier `PAACS::Mesh` spelling remains an
alias). This lookup never searches for OpenCASCADE, even when the SDK includes
the optional adapter. CAD consumers request the adapter explicitly:

```cmake
find_package(PAACSMesh 1.2.0 EXACT CONFIG REQUIRED COMPONENTS OCCT)
target_link_libraries(cad_application PRIVATE PAACS::occt)
```

`PAACS::MeshOCCT` is an alias for the loaded `PAACS::occt` target. These aliases
preserve CMake target names for local builds only; they do not restore the
previous V1/CAD source APIs or establish cross-version binary compatibility.

Use `OPTIONAL_COMPONENTS OCCT` to enable CAD only when both the installed
adapter and its OCCT dependency are available. Check `PAACSMesh_OCCT_FOUND`
before using `PAACS::occt`. The installed consumer example uses this optional
lookup. OCCT is exported separately in `PAACSOCCTTargets.cmake`.

The isolated package check disables OpenCASCADE discovery and compiles/runs a
core consumer. It also checks that an optional OCCT request succeeds without
the adapter target and that a required request fails:

```sh
cmake -DPAACS_INSTALL_PREFIX="$PWD/out/install" \
  -DPAACS_COMPONENT_TEST_BINARY_DIR="$PWD/build/package-components" \
  -DPAACS_COMPONENT_TEST_GENERATOR=Ninja \
  -DPAACS_COMPONENT_TEST_CXX_COMPILER=g++-16 \
  -P scripts/check_installed_components.cmake
```

For an OCCT-enabled installation, additionally pass
`-DPAACS_COMPONENT_TEST_OCCT_DIR=/path/to/occt/lib/cmake/opencascade` to verify
optional and required adapter discovery and compile/link its real toolkit
closure with the actual dependency present. The installed consumer example
supplies the separate CAD runtime and round-trip check.

The [candidate validation](RELEASE_READINESS.md) records completed CI and local
checks, including platform and sanitizer results. `VALIDATION.md` retains
earlier qualification records. Use those results when assessing a configuration;
a workflow or preset alone does not establish that its checks passed.

The geometry suite includes exact-audit differential checks and graded
body-control cases. `tests/cad_cli_checks.py` exercises the optional CAD CLI;
`benchmarks/run_audit_comparison.py` records paired audit measurements. Both
scripts use only Python's standard library. The gallery renderer optionally
needs matplotlib. On memory-limited builders, use `--parallel 1` for sanitizer
compilation; compilation memory is separate from meshing runtime memory.

## Automatic STEP build

Install CMake 3.30 or newer and your platform's C++23 build tools. Use Visual
Studio with the C++ workload on Windows, GCC or Clang on Linux, or the Xcode
command-line tools on macOS. The first build needs network access to download
OCCT and time to compile it. Subsequent builds reuse the local dependency.
The macOS build was checked on Apple Silicon with Xcode 16.4. CMake enables
the SDK's thread cancellation support when it requires `-fexperimental-library`.

```sh
cmake -S cmake/step -B build/step
cmake --build build/step --config Release
ctest --test-dir build/step -C Release --output-on-failure
```

The build verifies the OCCT 8.0.1 source archive against its pinned SHA-512.
It builds the STEP toolkits and their dependencies as shared libraries in
`build/step/occt`. The PAACS executables go in `build/step/bin`. Set
`-DPAACS_STEP_JOBS=2` when configuring to reduce parallel build memory. Standard
CMake generator selection applies. Use `-G Ninja` if Ninja and your compiler
are already on `PATH`.

For an installed PAACS C++ SDK with STEP support, follow the regular SDK
build and install steps with `PAACS_WITH_OCCT=ON`. Point `OpenCASCADE_DIR`
to the absolute path of `build/step/occt/lib/cmake/opencascade` to reuse the
locally built OCCT SDK. CAD consumers can use that same path.

For direct use of the STEP tools on Windows, put the local OCCT DLL directory
on `PATH` in the current PowerShell session. The CTest fixture sets it itself.

```powershell
$env:PATH = "${PWD}/build/step/occt/bin;$env:PATH"
build/step/bin/paacs-cad-mesh.exe tests/fixtures/cad/cylinder.step .7 cylinder.paacsmesh --native-cad --report-json
```

The [native STEP job](../.github/workflows/ci.yml) uses these same commands for
Windows, Linux and macOS. Its Windows run also checks sizing options and the
tagged C3D4 export. The root CMake project still offers the standalone core
build with `PAACS_WITH_OCCT=OFF`.

### Build with an existing OCCT SDK

Install Visual Studio 2022 with the C++ workload, CMake and an x64 OCCT 8.0.1
SDK. Set `OpenCASCADE_DIR` to the directory containing the SDK's CMake config.
Prepend the SDK and dependency DLL directories to `PATH`.

```powershell
cmake -S . -B build/step-sdk -G 'Visual Studio 17 2022' -A x64 `
  -DPAACS_CXX_STANDARD=23 -DPAACS_WITH_OCCT=ON -DPAACS_BUILD_TESTS=OFF `
  -DOpenCASCADE_DIR="$env:OpenCASCADE_DIR"
cmake --build build/step-sdk --config Release --target paacs-cad-mesh paacs-native-cad-preview
python tests/native_cad_cli_checks.py build/step-sdk/Release/paacs-cad-mesh.exe
```

When running a Windows application linked to the installed UI bridge, add
the installation's `bin` directory to `PATH` so it can load `paacs_ui.dll`.

### Local Windows validation

The [local validation record](LOCAL_VALIDATION_2026-10-03.md) contains the exact
MSVC/OCCT commands and results. CMake 4.0.2 with MSVC 19.44 does not recognize
`cxx_std_26`; use an explicit `-DPAACS_CXX_STANDARD=23` compatibility build on
that toolchain. The default C++26 requirement is retained.

`scripts/setup_occt_sdk.ps1` extracts the required subset of the supplied OCCT
8.0.1 Windows x64 archive into `build/deps/occt-8.0.1`. With `-ActivateRuntime`,
it also sets the matching DLL paths for the current PowerShell process. SDK
binaries and third-party dependencies are build inputs and are not included
in the mesher source package.

When using this SDK with Visual Studio, configure with
`"-DCMAKE_VS_GLOBALS=VcpkgEnabled=false;VcpkgApplocalDeps=false"`.
Global MSBuild vcpkg integration otherwise copied an unrelated OCCT 7.9 runtime
beside the executables on this machine, overriding the selected 8.0.1 DLLs in
`PATH`. Reconfigure and remove those verified stale build-output copies before
retesting; changing `PATH` alone cannot override adjacent DLLs.

Enable `PAACS_BUILD_EXTENDED_PROBES` to compile the retained targeted allocation
campaigns; run them with `ctest -L extended`. Full default tests and the new
archive/projection/CAD-runner regressions remain separate from these campaigns.
