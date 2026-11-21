#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Build, repair and test the first Windows x64 CPython 3.12 STEP wheel.

Requires the ``build`` and ``delvewheel`` build tools. Build isolation installs
the pinned backend/binding headers; consumers install only the repaired wheel.
No package index publication is performed.
"""
from __future__ import annotations

import argparse
import hashlib
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
        if (candidate / "bin" / "TKernel.dll").is_file():
            return candidate
    raise RuntimeError("Use a shared OCCT SDK with bin/TKernel.dll and lib/cmake/opencascade")


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


def source_archive(source: Path, sdk: Path, output: Path) -> Path:
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


def validate_contents(wheel: Path, runtime: Path, crt: list[str]) -> dict:
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        if not any(name.startswith("paacs_mesh/_native") and name.endswith(".pyd") for name in names):
            raise RuntimeError("The wheel has no native STEP binding")
        if any(name.endswith((".lib", ".hpp", ".hxx", ".cmake")) for name in names):
            raise RuntimeError("The wheel unexpectedly contains C++ SDK development files")
        dlls = sorted(name for name in names if name.lower().endswith(".dll"))
        if not any(Path(name).name.lower().startswith("tkdestep") for name in dlls):
            raise RuntimeError("The repaired wheel does not contain the OCCT STEP library")
        # Microsoft redist files must remain unmodified: exclude them from name mangling.
        for name in crt:
            matches = [entry for entry in dlls if Path(entry).name.lower() == name]
            if len(matches) != 1 or hashlib.sha256(archive.read(matches[0])).hexdigest() != sha256(runtime / name):
                raise RuntimeError(f"Redistributable was omitted or changed: {name}")
        for required in ("LICENSE_LGPL_21.txt", "OCCT_LGPL_EXCEPTION.txt", "OCCT_NOTICE.txt"):
            if not any(name.endswith("/" + required) for name in names):
                raise RuntimeError(f"The wheel omits {required}")
        return {"wheel": wheel.name, "sha256": sha256(wheel), "dlls": dlls,
                "unmodified_microsoft_runtime": crt}


def installed_smoke(source: Path, wheel: Path, build_dir: Path) -> None:
    environment = build_dir / "installed-test"
    venv.EnvBuilder(with_pip=True, clear=True).create(environment)
    python = environment / "Scripts/python.exe"
    run([str(python), "-m", "pip", "install", "--no-index", "--no-deps", str(wheel)], cwd=environment)
    fixture = environment / "cylinder.step"
    shutil.copyfile(source / "tests/fixtures/cad/cylinder.step", fixture)
    shutil.copyfile(source / "tests/fixtures/cad/assembly.step", environment / "assembly.step")
    clean = {key: value for key, value in os.environ.items()
             if key.upper() not in {"PYTHONPATH", "PYTHONHOME", "CASROOT", "OpenCASCADE_DIR".upper()}
             and not key.upper().startswith("CSF_")}
    # Windows runners have build tools installed; this proves no SDK/VC PATH is needed.
    windows = Path(os.environ.get("SystemRoot", "C:/Windows"))
    clean["PATH"] = os.pathsep.join([str(environment / "Scripts"), str(windows / "System32"), str(windows)])
    clean["PAACS_TEST_STEP_FIXTURE"] = str(fixture)
    clean["PAACS_EXPECT_INSTALLED"] = "1"
    clean["PYTHONNOUSERSITE"] = "1"
    run([str(python), "-I", "-c", "import paacs_mesh; print(paacs_mesh.__file__)"], cwd=environment, env=clean)
    result = run([str(python), "-I", "-m", "unittest", "discover", "-s", str(source / "python/tests"), "-v"],
                 cwd=environment, env=clean)
    if "Ran 0 tests" in result:
        raise RuntimeError("Installed-wheel qualification did not discover any tests")


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
    if (sys.platform != "win32" or platform.machine().lower() not in {"amd64", "x86_64"}
            or struct.calcsize("P") != 8):
        parser.error("The first supported wheel is Windows x64")
    if sys.implementation.name != "cpython" or sys.version_info[:2] != (3, 12):
        parser.error("Run this milestone's wheel builder with CPython 3.12 x64")
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
    env = dict(os.environ, CMAKE_GENERATOR="Visual Studio 17 2022", CMAKE_GENERATOR_PLATFORM="x64",
               CMAKE_BUILD_PARALLEL_LEVEL=str(args.jobs))
    run(command, cwd=source, env=env)
    wheels = list(raw.glob("*.whl"))
    if len(wheels) != 1:
        raise RuntimeError(f"Expected one raw wheel in {raw}; use a dedicated build directory")
    sdk = sdk_prefix(config_dir)
    runtime = runtime_directory(args.vc_runtime_dir)
    crt = crt_dependencies(wheels[0], sdk, runtime)
    paths = os.pathsep.join((str(sdk / "bin"), str(runtime)))
    run([sys.executable, "-m", "delvewheel", "show", str(wheels[0]), "--add-path", paths], cwd=source)
    repair = [sys.executable, "-m", "delvewheel", "repair", str(wheels[0]), "--add-path", paths,
              "--wheel-dir", str(output)]
    if crt:
        repair.extend(["--include", os.pathsep.join(crt), "--no-mangle", os.pathsep.join(crt)])
    run(repair, cwd=source)
    wheel = output / wheels[0].name
    receipt = validate_contents(wheel, runtime, crt)
    archive = source_archive(source, sdk, output)
    receipt["occt_source"] = {"file": archive.name, "sha256": sha256(archive)}
    installed_smoke(source, wheel, build_dir)
    receipt["installed_smoke"] = "passed with isolated Python, SDK-free PATH and cleared CAD variables"
    (output / "wheel-validation.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
