#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Build, repair and qualify a CPython 3.12 STEP wheel with its CAD runtime.

Supports Windows x64, macOS 15 ARM64 and manylinux 2.28 x64. Build isolation
installs the pinned backend/binding headers. The platform repair tool bundles
the shared CAD dependencies; no package index publication is performed.
"""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version as distribution_version
import json
import os
from pathlib import Path
import platform
import re
import shutil
import struct
import subprocess
import sys
import urllib.request
import venv
import zipfile


CRT_NAMES = {
    "msvcp140.dll", "msvcp140_1.dll", "msvcp140_2.dll",
    "msvcp140_atomic_wait.dll", "msvcp140_codecvt_ids.dll", "concrt140.dll",
    "vcruntime140.dll", "vcruntime140_1.dll", "vcruntime140_threads.dll",
}
NUMPY_TEST_VERSIONS = ("1.26.4", "2.5.5")
PLATFORM_TAGS = {
    "windows-x86_64": "win_amd64",
    "macos-arm64": "macosx_15_0_arm64",
    "linux-x86_64": "manylinux_2_28_x86_64",
}


def source_revision(source: Path) -> dict:
    """Record source identity; unversioned developer builds cannot be published."""
    try:
        def git(*args):
            return subprocess.check_output(["git", "-c", "safe.directory=" + source.resolve().as_posix(),
                                           "-C", str(source), *args],
                                           text=True, stderr=subprocess.DEVNULL).strip()
        return {"commit": git("rev-parse", "HEAD"), "tree": git("rev-parse", "HEAD^{tree}"),
                "clean": not git("status", "--porcelain", "--untracked-files=normal")}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "tree": None, "clean": False}


def host_platform(*, building: bool = False) -> str:
    if (sys.implementation.name != "cpython" or sys.version_info[:2] != (3, 12)
            or struct.calcsize("P") != 8):
        raise RuntimeError("Wheel qualification requires 64-bit CPython 3.12")
    architecture = platform.machine().lower()
    if sys.platform == "win32" and architecture in {"amd64", "x86_64"}:
        return "windows-x86_64"
    if sys.platform == "darwin" and architecture in {"arm64", "aarch64"}:
        if int(platform.mac_ver()[0].split(".")[0]) < 15:
            raise RuntimeError("The supported macOS wheel requires macOS 15 or newer")
        return "macos-arm64"
    if sys.platform == "linux" and architecture in {"amd64", "x86_64"}:
        libc, release = platform.libc_ver()
        if libc != "glibc" or tuple(map(int, release.split(".")[:2])) < (2, 28):
            raise RuntimeError("The supported Linux wheel requires glibc 2.28 or newer")
        if building and release != "2.28":
            raise RuntimeError("Build Linux wheels inside the pinned manylinux 2.28 x64 image")
        return "linux-x86_64"
    raise RuntimeError("Supported wheels are Windows x64, macOS 15 ARM64 and Linux x64")


def wheel_platform(wheel: Path) -> str:
    fields = wheel.name.removesuffix(".whl").rsplit("-", 3)
    if len(fields) != 4 or fields[1:3] != ["cp312", "cp312"]:
        raise RuntimeError("Expected a CPython 3.12 wheel with its native ABI")
    return fields[3]


def run(command: list[str], *, cwd: Path, env: dict[str, str] | None = None) -> str:
    print("+ " + subprocess.list2cmdline(command), flush=True)
    lines = []
    with subprocess.Popen(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True,
                          encoding="utf-8", errors="replace") as process:
        for line in process.stdout:
            print(line, end="", flush=True)
            lines.append(line)
        status = process.wait()
    output = "".join(lines)
    if status:
        raise subprocess.CalledProcessError(status, command, output=output)
    return output


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def sdk_prefix(config_dir: Path) -> Path:
    for candidate in (config_dir.parent.parent.parent, config_dir.parent):
        if ((candidate / "bin" / "TKernel.dll").is_file()
                or any((candidate / "lib").glob("libTKernel.so*"))
                or any((candidate / "lib").glob("libTKernel*.dylib"))):
            return candidate
    raise RuntimeError("Use a shared OCCT SDK with its TKernel library and CMake configuration")


def runtime_directory(explicit: Path | None) -> Path:
    candidates = [explicit] if explicit else []
    if os.environ.get("VCToolsRedistDir"):
        candidates.append(Path(os.environ["VCToolsRedistDir"]) / "x64/Microsoft.VC143.CRT")
    program_files = Path(os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)"))
    vswhere = program_files / "Microsoft Visual Studio/Installer/vswhere.exe"
    if vswhere.is_file():
        installation = subprocess.check_output(
            [str(vswhere), "-latest", "-products", "*", "-property", "installationPath"],
            text=True).strip()
        if installation:
            redist = Path(installation) / "VC/Redist/MSVC"
            candidates.extend(sorted(redist.glob("*/x64/Microsoft.VC143.CRT"), reverse=True))
    for candidate in candidates:
        if candidate is not None and (candidate / "msvcp140.dll").is_file():
            return candidate.resolve()
    raise RuntimeError("Visual C++ redistributable files not found; set --vc-runtime-dir")


def occt_pin(source: Path) -> dict[str, str]:
    pin = (source / "cmake/OCCTVersion.cmake").read_text(encoding="utf-8")
    values = {}
    for key in ("VERSION", "URL", "SHA512"):
        match = re.search(rf'set\(PAACS_OCCT_{key}\s+"?([^"\s)]+)', pin)
        if not match:
            raise RuntimeError(f"Missing OCCT pin {key}")
        values[key] = match.group(1)
    # The URL may reference the shared version variable rather than repeat it.
    values["URL"] = values["URL"].replace("${PAACS_OCCT_VERSION}", values["VERSION"])
    values["URL"] = values["URL"].replace("${PAACS_OCCT_TAG}", "V" + values["VERSION"].replace(".", "_"))
    return values


def source_archive(source: Path, sdk: Path, output: Path) -> Path:
    values = occt_pin(source)
    name = f"occt-{values['VERSION']}.tar.gz"
    target = output / name
    original = sdk.parent / "downloads" / name
    if original.is_file():
        shutil.copyfile(original, target)
    elif not target.is_file():
        print(f"Downloading matching OCCT source: {values['URL']}", flush=True)
        with urllib.request.urlopen(values["URL"], timeout=120) as response, target.open("wb") as stream:
            shutil.copyfileobj(response, stream)
    with target.open("rb") as stream:
        actual = hashlib.file_digest(stream, "sha512").hexdigest()
    if actual.lower() != values["SHA512"].lower():
        raise RuntimeError("OCCT source archive does not match the pinned SHA-512")
    return target


def crt_dependencies(wheel: Path, sdk: Path, runtime: Path) -> list[str]:
    import pefile  # Installed by delvewheel; inspect actual binary imports.
    required = set()
    binaries = list((sdk / "bin").glob("*.dll"))
    with zipfile.ZipFile(wheel) as archive:
        payloads = [archive.read(name) for name in archive.namelist() if name.endswith(".pyd")]
    payloads.extend(binary.read_bytes() for binary in binaries)
    for payload in payloads:
        image = pefile.PE(data=payload, fast_load=True)
        image.parse_data_directories(directories=[pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_IMPORT"]])
        for entry in getattr(image, "DIRECTORY_ENTRY_IMPORT", ()):
            name = entry.dll.decode("ascii").lower()
            if name in CRT_NAMES:
                required.add(name)
        image.close()
    for name in required:
        if not (runtime / name).is_file():
            raise RuntimeError(f"Required redistributable is missing: {runtime / name}")
    return sorted(required)


def validate_contents(wheel: Path, runtime: Path | None = None,
                      crt: list[str] | tuple = ()) -> dict:
    current_platform = host_platform()
    if PLATFORM_TAGS[current_platform] not in wheel_platform(wheel).split("."):
        raise RuntimeError("The wheel platform tag does not match the supported host")
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        extension = ".pyd" if sys.platform == "win32" else ".so"
        if not any(name.startswith("paacs_mesh/_native") and name.endswith(extension) for name in names):
            raise RuntimeError("The wheel has no native STEP binding")
        if any(name.endswith((".lib", ".a", ".hpp", ".hxx", ".cmake", ".exe"))
               or Path(name).name in {"paacs-mesh", "paacs-cad-mesh", "paacs-native-cad-preview"}
               for name in names):
            raise RuntimeError("The wheel unexpectedly contains C++ SDK development or tool files")
        dlls = sorted(name for name in names if name.lower().endswith(".dll"))
        libraries = sorted(name for name in names if name.lower().endswith((".dll", ".dylib"))
                           or re.search(r"\.so(?:\.|$)", name) and "/_native" not in name)
        if not any(Path(name).name.lower().startswith(("tkdestep", "libtkdestep")) for name in libraries):
            raise RuntimeError("The repaired wheel does not contain the OCCT STEP library")
        if not any(Path(name).name.lower().startswith(("tkernel", "libtkernel")) for name in libraries):
            raise RuntimeError("The repaired wheel does not contain the OCCT kernel library")
        # Microsoft redist files must remain unmodified: exclude them from name mangling.
        for name in crt:
            matches = [entry for entry in dlls if Path(entry).name.lower() == name]
            if (runtime is None or len(matches) != 1
                    or hashlib.sha256(archive.read(matches[0])).hexdigest() != sha256(runtime / name)):
                raise RuntimeError(f"Redistributable was omitted or changed: {name}")
        for required in ("LICENSE_LGPL_21.txt", "OCCT_LGPL_EXCEPTION.txt", "OCCT_NOTICE.txt",
                         "pybind11-LICENSE.txt", "LICENSE", "NOTICE"):
            if not any(name.endswith("/" + required) for name in names):
                raise RuntimeError(f"The wheel omits {required}")
        for resource in ("BOPAlgo", "SHMessage", "XSMessage", "XSTEPResource", "UnitsAPI", "StdResource"):
            prefix = f"paacs_mesh/_occt/resources/{resource}/"
            if not any(name.startswith(prefix) and not name.endswith("/") for name in names):
                raise RuntimeError(f"The wheel omits the OCCT {resource} resources")
        result = {"wheel": wheel.name, "sha256": sha256(wheel), "platform": current_platform,
                  "platform_tag": PLATFORM_TAGS[current_platform], "libraries": libraries}
        if sys.platform == "win32":
            result.update(dlls=dlls, unmodified_microsoft_runtime=list(crt))
        return result


def repair_wheel(wheel: Path, sdk: Path, output: Path, source: Path,
                 runtime: Path | None) -> tuple[Path, dict]:
    environment = dict(os.environ)
    crt = []
    if sys.platform == "win32":
        if runtime is None:
            raise RuntimeError("The Windows repair requires the Microsoft runtime directory")
        crt = crt_dependencies(wheel, sdk, runtime)
        paths = os.pathsep.join((str(sdk / "bin"), str(runtime)))
        run([sys.executable, "-m", "delvewheel", "show", str(wheel), "--add-path", paths], cwd=source)
        command = [sys.executable, "-m", "delvewheel", "repair", str(wheel), "--add-path", paths,
                   "--wheel-dir", str(output)]
        if crt:
            command.extend(["--include", os.pathsep.join(crt), "--no-mangle", os.pathsep.join(crt)])
        audit = [sys.executable, "-m", "delvewheel", "show"]
        tool = "delvewheel"
    elif sys.platform == "darwin":
        environment["DYLD_LIBRARY_PATH"] = str(sdk / "lib")
        environment["MACOSX_DEPLOYMENT_TARGET"] = "15.0"
        command = ["delocate-wheel", "--require-archs", "arm64", "-v", "-w", str(output), str(wheel)]
        audit = ["delocate-listdeps", "--all"]
        tool = "delocate"
    else:
        environment["LD_LIBRARY_PATH"] = str(sdk / "lib")
        command = [sys.executable, "-m", "auditwheel", "repair", "--plat", "manylinux_2_28_x86_64",
                   "--wheel-dir", str(output), str(wheel)]
        audit = [sys.executable, "-m", "auditwheel", "show"]
        tool = "auditwheel"
    run(command, cwd=source, env=environment)
    prefix = wheel.name.rsplit("-", 1)[0] + "-"
    candidates = [path for path in output.glob(prefix + "*.whl")
                  if PLATFORM_TAGS[host_platform()] in wheel_platform(path).split(".")]
    if len(candidates) != 1:
        raise RuntimeError("Expected one repaired wheel for this build; use a dedicated output directory")
    repaired = candidates[0]
    # Inspect the repaired dependency closure without adding SDK search paths.
    audit_environment = {key: value for key, value in os.environ.items()
                         if key.upper() not in {"LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT"}
                         and not key.upper().startswith("DYLD_")}
    audit_output = run([*audit, str(repaired)], cwd=source, env=audit_environment)
    if str(sdk) in audit_output:
        raise RuntimeError("The repaired wheel still reports a dependency on the build SDK")
    (output / "dependency-audit.txt").write_text(audit_output, encoding="utf-8")
    receipt = validate_contents(repaired, runtime, crt)
    receipt["repair"] = {"tool": tool, "version": distribution_version(tool),
                         "audit": "dependency-audit.txt",
                         "audit_sha256": sha256(output / "dependency-audit.txt")}
    return repaired, receipt


def isolated_environment(environment: Path) -> dict[str, str]:
    clean = {key: value for key, value in os.environ.items()
             if key.upper() not in {"PYTHONPATH", "PYTHONHOME", "CASROOT",
                                    "OPENCASCADE_DIR", "LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT",
                                    "LIBRARY_PATH"}
             and not key.upper().startswith(("CSF_", "DYLD_"))}
    if sys.platform == "win32":
        windows = Path(os.environ.get("SystemRoot", "C:/Windows"))
        clean["PATH"] = os.pathsep.join([str(environment / "Scripts"), str(windows / "System32"), str(windows)])
    else:
        clean["PATH"] = os.pathsep.join([str(environment / "bin"), "/usr/bin", "/bin"])
    clean["PYTHONNOUSERSITE"] = "1"
    return clean


def prepare_native_runtime(native_cli: Path, sdk: Path, runtime: Path | None) -> None:
    if sys.platform != "win32":
        return
    # An unrelated system-wide MSBuild package integration can copy another
    # OCCT installation beside the executable, ahead of the requested SDK.
    for binary in native_cli.parent.glob("TK*.dll"):
        original = sdk / "bin" / binary.name
        if not original.is_file() or sha256(binary) != sha256(original):
            raise RuntimeError("The native CLI directory contains a different OCCT runtime: "
                               f"{binary.name}. Use a clean build directory with MSBuild "
                               "package integration disabled.")
    if runtime is None:
        raise RuntimeError("The native CLI requires the matching Microsoft runtime")
    # The CLI is a build-time parity tool. Give it the same unmodified runtime
    # used by the repaired wheel, independent of the machine's system version.
    for name in CRT_NAMES:
        original = runtime / name
        if original.is_file():
            shutil.copyfile(original, native_cli.parent / name)


def installed_smoke(source: Path, wheel: Path, build_dir: Path) -> list[dict]:
    environment = build_dir / "installed-test"
    venv.EnvBuilder(with_pip=True, clear=True).create(environment)
    python = environment / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    run([str(python), "-m", "pip", "install", "--no-index", "--no-deps", str(wheel)], cwd=environment)
    console = environment / ("Scripts/paacs-mesh.exe" if sys.platform == "win32" else "bin/paacs-mesh")
    if not console.is_file():
        raise RuntimeError("The installed package omits the paacs-mesh command")
    run([str(console), "--help"], cwd=environment, env=isolated_environment(environment))
    fixture = environment / "cylinder.step"
    shutil.copyfile(source / "tests/fixtures/cad/cylinder.step", fixture)
    shutil.copyfile(source / "tests/fixtures/cad/assembly.step", environment / "assembly.step")
    clean = isolated_environment(environment)
    clean["PAACS_TEST_STEP_FIXTURE"] = str(fixture)
    clean["PAACS_EXPECT_INSTALLED"] = "1"
    clean["PAACS_EXPECT_NO_NUMPY"] = "1"
    clean.pop("PAACS_REQUIRE_NUMPY", None)
    sizing_results = build_dir / "sizing-validation.json"
    if sizing_results.is_file():
        sizing_results.unlink()
    clean["PAACS_SIZING_RESULTS"] = str(sizing_results)
    run([str(python), "-I", "-c", "import paacs_mesh; print(paacs_mesh.__file__)"], cwd=environment, env=clean)
    run([str(python), "-I", "-c",
         "import importlib.util; assert importlib.util.find_spec('numpy') is None"],
        cwd=environment, env=clean)
    result = run([str(python), "-I", "-m", "unittest", "discover", "-s", str(source / "python/tests"), "-v"],
                 cwd=environment, env=clean)
    if "Ran 0 tests" in result:
        raise RuntimeError("Installed-wheel qualification did not discover any tests")
    if (source / "python/tests/test_sizing.py").is_file() and not sizing_results.is_file():
        raise RuntimeError("Installed-wheel sizing checks did not produce their measurements")
    count = re.search(r"^Ran (\d+) tests? in ", result, re.MULTILINE)
    if count is None or int(count.group(1)) == 0:
        raise RuntimeError("Installed-wheel qualification did not report its test count")
    skipped = re.search(r"^OK \(skipped=(\d+)\)", result, re.MULTILINE)
    base_tests = {"discovered": int(count.group(1)),
                  "skipped": int(skipped.group(1)) if skipped else 0}
    base_tests["passed"] = base_tests["discovered"] - base_tests["skipped"]
    (build_dir / "installed-validation.json").write_text(
        json.dumps(base_tests, indent=2) + "\n", encoding="utf-8")
    numpy_results = []
    if (source / "python/tests/test_numpy.py").is_file():
        clean.pop("PAACS_EXPECT_NO_NUMPY", None)
        clean["PAACS_REQUIRE_NUMPY"] = "1"
        for version in NUMPY_TEST_VERSIONS:
            run([str(python), "-m", "pip", "install", "--disable-pip-version-check",
                 "--only-binary=:all:", "--no-deps", "--force-reinstall", f"numpy=={version}"],
                cwd=environment, env=clean)
            run([str(python), "-I", "-c",
                 f"import numpy; assert numpy.__version__ == {version!r}; print(numpy.__version__)"],
                cwd=environment, env=clean)
            result = run([str(python), "-I", "-m", "unittest", "discover",
                          "-s", str(source / "python/tests"), "-p", "test_numpy.py", "-v"],
                         cwd=environment, env=clean)
            count = re.search(r"^Ran (\d+) tests? in ", result, re.MULTILINE)
            if count is None or int(count.group(1)) == 0 or re.search(r"^OK \(skipped=", result, re.MULTILINE):
                raise RuntimeError("NumPy qualification must discover and pass every array test")
            numpy_results.append({"version": version, "tests": int(count.group(1)), "status": "passed"})
    return numpy_results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--occt-dir", type=Path, help="Reuse a source-built shared SDK's CMake config directory")
    parser.add_argument("--vc-runtime-dir", type=Path, help="Licensed x64 Microsoft.VC143.CRT directory")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--no-isolation", action="store_true", help="Use already installed pinned build requirements")
    args = parser.parse_args()
    try:
        current_platform = host_platform(building=True)
    except RuntimeError as error:
        parser.error(str(error))
    if args.vc_runtime_dir and sys.platform != "win32":
        parser.error("--vc-runtime-dir applies only to Windows wheels")
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    source = args.source.resolve()
    build_dir = (args.build_dir or source / "build/python-wheel").resolve()
    output = (args.out_dir or source / "dist").resolve()
    raw = build_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, "-m", "build", "--wheel", str(source), "--outdir", str(raw),
               f"-Cbuild-dir={build_dir / 'cmake'}",
               "-Ccmake.define.PAACS_BUILD_CLI=ON",
               "-Cbuild.targets=_native;paacs-cad-mesh",
               f"-Ccmake.define.PAACS_PYTHON_OCCT_JOBS={args.jobs}"]
    if args.no_isolation:
        command.append("--no-isolation")
    if args.occt_dir:
        config_dir = args.occt_dir.resolve()
        if not (config_dir / "OpenCASCADEConfig.cmake").is_file():
            parser.error("--occt-dir must contain OpenCASCADEConfig.cmake")
        command.extend([f"-Ccmake.define.OpenCASCADE_DIR={config_dir}", "-Ccmake.define.PAACS_PYTHON_BUILD_OCCT=OFF"])
    else:
        config_dir = build_dir / "cmake/_occt-build/occt/lib/cmake/opencascade"
    env = dict(os.environ, CMAKE_BUILD_PARALLEL_LEVEL=str(args.jobs))
    if sys.platform == "win32":
        env.update(CMAKE_GENERATOR="Visual Studio 17 2022", CMAKE_GENERATOR_PLATFORM="x64")
        command.append("-Ccmake.define.CMAKE_VS_GLOBALS=VcpkgEnabled=false;VcpkgApplocalDeps=false")
    else:
        env["CMAKE_GENERATOR"] = "Ninja"
        env.pop("CMAKE_GENERATOR_PLATFORM", None)
        env.pop("CMAKE_GENERATOR_TOOLSET", None)
    if sys.platform == "darwin":
        env["MACOSX_DEPLOYMENT_TARGET"] = "15.0"
        env["ARCHFLAGS"] = "-arch arm64"
        command.extend(["-Ccmake.define.CMAKE_OSX_ARCHITECTURES=arm64",
                        "-Ccmake.define.CMAKE_OSX_DEPLOYMENT_TARGET=15.0"])
    run(command, cwd=source, env=env)
    wheels = list(raw.glob("*.whl"))
    if len(wheels) != 1:
        raise RuntimeError(f"Expected one raw wheel in {raw}; use a dedicated build directory")
    sdk = sdk_prefix(config_dir)
    # Both targets are built before the backend's isolated environment is
    # removed. Only the Python component is installed into the wheel.
    native_cli = (build_dir / "cmake" / "Release/paacs-cad-mesh.exe" if sys.platform == "win32"
                  else build_dir / "cmake/paacs-cad-mesh")
    if not native_cli.is_file():
        raise RuntimeError("The native parity executable was not built")
    runtime = runtime_directory(args.vc_runtime_dir) if sys.platform == "win32" else None
    prepare_native_runtime(native_cli, sdk, runtime)
    wheel, receipt = repair_wheel(wheels[0], sdk, output, source, runtime)
    receipt["source_revision"] = source_revision(source)
    if current_platform == "linux-x86_64":
        receipt["minimum_glibc"] = "2.28"
        receipt["build_image"] = os.environ.get("PAACS_MANYLINUX_IMAGE")
    elif current_platform == "macos-arm64":
        receipt["minimum_macos"] = "15.0"
    archive = source_archive(source, sdk, output)
    receipt["occt_source"] = {"file": archive.name, "sha256": sha256(archive)}
    numpy_results = installed_smoke(source, wheel, build_dir)
    environment = build_dir / "installed-test"
    python = environment / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    parity = output / "native-parity.json"
    run([str(python), "-I", str(source / "scripts/check_python_parity.py"),
         "--native-cli", str(native_cli), "--sdk-lib-dir", str(sdk / ("bin" if sys.platform == "win32" else "lib")),
         "--output", str(parity)], cwd=environment, env=isolated_environment(environment))
    if not parity.is_file():
        raise RuntimeError("Native CLI parity checks did not produce their receipt")
    receipt["native_parity"] = {"file": parity.name, "sha256": sha256(parity)}
    receipt["installed_tests"] = json.loads((build_dir / "installed-validation.json").read_text(encoding="utf-8"))
    if numpy_results:
        receipt["numpy_validation"] = numpy_results
    sizing_results = build_dir / "sizing-validation.json"
    if sizing_results.is_file():
        sizing_output = output / sizing_results.name
        shutil.copyfile(sizing_results, sizing_output)
        receipt["sizing_validation"] = {"file": sizing_output.name, "sha256": sha256(sizing_output)}
    receipt["installed_smoke"] = "passed with isolated Python and cleared SDK, loader and CAD paths"
    (output / "wheel-validation.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
