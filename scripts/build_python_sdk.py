#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Build or verify a toolchain-bound shared OCCT SDK, independent of Python ABI."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import tarfile

from build_python_wheel import cmake_generator, host_platform, occt_pin, run, sha256, validate_binary


SCHEMA = "paacs-python-occt-sdk-v1"
RECIPES = ("cmake/OCCTVersion.cmake", "cmake/step/CMakeLists.txt", "cmake/python/OCCT.cmake",
           "scripts/build_python_sdk.py")


def context(source: Path, work: Path) -> dict:
    """Probe the selected compiler, rather than keying only on an OS label."""
    probe_source, probe_build = work / "probe-source", work / "probe-build"
    probe_source.mkdir(parents=True, exist_ok=True)
    (probe_source / "CMakeLists.txt").write_text(
        'cmake_minimum_required(VERSION 3.30)\nproject(PaacsSDKProbe LANGUAGES C CXX)\n'
        'file(WRITE "${CMAKE_BINARY_DIR}/compiler.txt" '
        '"id=${CMAKE_CXX_COMPILER_ID}\\nversion=${CMAKE_CXX_COMPILER_VERSION}\\n'
        'compiler=${CMAKE_CXX_COMPILER}\\ngenerator=${CMAKE_GENERATOR}\\n'
        'c_id=${CMAKE_C_COMPILER_ID}\\nc_version=${CMAKE_C_COMPILER_VERSION}\\n'
        'c_compiler=${CMAKE_C_COMPILER}\\nlinker=${CMAKE_LINKER}\\n'
        'linker_flags=${CMAKE_SHARED_LINKER_FLAGS}\\n'
        'release_linker_flags=${CMAKE_SHARED_LINKER_FLAGS_RELEASE}\\n'
        'msvc_runtime=${CMAKE_MSVC_RUNTIME_LIBRARY}\\n'
        'platform=${CMAKE_GENERATOR_PLATFORM}\\ntoolset=${CMAKE_GENERATOR_TOOLSET}\\n'
        'msvc_toolset=${MSVC_TOOLSET_VERSION}\\nwindows_sdk=${CMAKE_VS_WINDOWS_TARGET_PLATFORM_VERSION}\\n'
        'system=${CMAKE_SYSTEM_NAME}\\nprocessor=${CMAKE_SYSTEM_PROCESSOR}\\n'
        'c_flags=${CMAKE_C_FLAGS}\\ncxx_flags=${CMAKE_CXX_FLAGS}\\n'
        'release_flags=${CMAKE_CXX_FLAGS_RELEASE}\\nc_release_flags=${CMAKE_C_FLAGS_RELEASE}\\n'
        'pointer_bytes=${CMAKE_SIZEOF_VOID_P}\\n")\n', encoding="utf-8")
    arguments = cmake_generator()
    command = ["cmake", "-S", str(probe_source), "-B", str(probe_build), *arguments]
    if platform.system() == "Darwin":
        command += ["-DCMAKE_OSX_ARCHITECTURES=arm64", "-DCMAKE_OSX_DEPLOYMENT_TARGET=15.0"]
    run(command, cwd=source)
    compiler = dict(line.split("=", 1) for line in
                    (probe_build / "compiler.txt").read_text(encoding="utf-8").splitlines())
    if compiler.get("pointer_bytes") != "8":
        raise RuntimeError("The reusable SDK requires a 64-bit compiler")
    sdk = subprocess.check_output(["xcrun", "--show-sdk-version"], text=True).strip() if platform.system() == "Darwin" else None
    os_release = (Path("/etc/os-release").read_text(encoding="utf-8") if platform.system() == "Linux"
                  else platform.mac_ver()[0] if platform.system() == "Darwin" else platform.win32_ver()[1])
    return {"schema": SCHEMA, "platform": host_platform(building=True), "occt": occt_pin(source),
            "recipes": {name: sha256(source / name) for name in RECIPES}, "compiler": compiler,
            "os_release": os_release, "configuration": "Release", "library_type": "Shared",
            "cmake": subprocess.check_output(["cmake", "--version"], text=True).splitlines()[0],
            "macos_deployment_target": "15.0" if platform.system() == "Darwin" else None,
            "macos_sdk": sdk, "manylinux_image": os.environ.get("PAACS_MANYLINUX_IMAGE")}


def cache_key(descriptor: dict) -> str:
    digest = hashlib.sha256(json.dumps(descriptor, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return f"paacs-occt-v1-{descriptor['platform']}-{digest}"


def inventory(root: Path) -> dict:
    records = {}
    resolved = root.resolve()
    for entry in sorted(root.rglob("*")):
        if not entry.resolve().is_relative_to(resolved):
            raise RuntimeError("An SDK entry escapes its verified root")
        if entry.is_symlink() and (not entry.exists() or entry.is_dir()):
            raise RuntimeError("SDK directory and broken symlinks are not supported")
        if entry.is_file():
            records[entry.relative_to(root).as_posix()] = {
                "sha256": sha256(entry), "symlink": os.readlink(entry) if entry.is_symlink() else None}
    return records


def verify_sdk(directory: Path, expected: dict | None = None) -> dict:
    receipt = json.loads((directory / "sdk-receipt.json").read_text(encoding="utf-8"))
    if receipt.get("schema") != SCHEMA or (expected is not None and receipt.get("context") != expected):
        raise RuntimeError("Cached SDK has another source recipe, architecture or toolchain")
    sdk = directory / "occt"
    files = inventory(sdk)
    if not files or files != receipt.get("files"):
        raise RuntimeError("SDK contents differ from their recorded hashes")
    if "lib/cmake/opencascade/OpenCASCADEConfig.cmake" not in files:
        raise RuntimeError("SDK has no relocatable OCCT configuration")
    names = {Path(name).name for name in files}
    for license_name in ("LICENSE_LGPL_21.txt", "OCCT_LGPL_EXCEPTION.txt"):
        if license_name not in names:
            raise RuntimeError(f"SDK omits {license_name}")
    binaries = [name for name in files if name.lower().endswith((".dll", ".dylib"))
                or re.search(r"\.so(?:\.|$)", name)]
    if not any(Path(name).name.lower().startswith(("tkernel", "libtkernel")) for name in binaries):
        raise RuntimeError("SDK omits the native OCCT kernel")
    if not any(Path(name).name.lower().startswith(("tkdestep", "libtkdestep")) for name in binaries):
        raise RuntimeError("SDK omits the native STEP library")
    for name in binaries:
        validate_binary((sdk / name).read_bytes(), name, receipt["context"]["platform"])
    pin = receipt["context"]["occt"]
    archive = directory / "downloads" / f"occt-{pin['VERSION']}.tar.gz"
    with archive.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha512").hexdigest()
    if digest.lower() != pin["SHA512"].lower():
        raise RuntimeError("Cached SDK source differs from the pinned OCCT source")
    return receipt


def build_sdk(source: Path, work: Path, directory: Path, descriptor: dict, jobs: int):
    if (directory / "sdk-receipt.json").is_file():
        verify_sdk(directory, descriptor)
        print("Reused and hash-verified OCCT SDK", flush=True)
        return
    if directory.exists() and any(directory.iterdir()):
        raise RuntimeError("Incomplete SDK cache; use an empty output directory")
    native = work / "occt-build"
    command = ["cmake", "-S", str(source / "cmake/step"), "-B", str(native), *cmake_generator(),
               "-DCMAKE_BUILD_TYPE=Release", "-DPAACS_STEP_BUILD_TOOLS=OFF", f"-DPAACS_STEP_JOBS={jobs}"]
    if platform.system() == "Darwin":
        command += ["-DCMAKE_OSX_ARCHITECTURES=arm64", "-DCMAKE_OSX_DEPLOYMENT_TARGET=15.0"]
    run(command, cwd=source)
    run(["cmake", "--build", str(native), "--config", "Release", "--target", "occt", "--parallel", str(jobs)], cwd=source)
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copytree(native / "occt", directory / "occt", symlinks=True)
    (directory / "downloads").mkdir()
    pin = descriptor["occt"]
    shutil.copy2(native / "downloads" / f"occt-{pin['VERSION']}.tar.gz", directory / "downloads")
    receipt = {"schema": SCHEMA, "context": descriptor, "files": inventory(directory / "occt")}
    (directory / "sdk-receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    verify_sdk(directory, descriptor)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("key", "build", "verify", "pack", "unpack"))
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--work", type=Path, default=Path("build/python-sdk-work"))
    parser.add_argument("--sdk", type=Path, default=Path("build/python-sdk"))
    parser.add_argument("--archive", type=Path, default=Path("build/python-sdk.tar.gz"))
    parser.add_argument("--github-output", type=Path)
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    source, work, directory = args.source.resolve(), args.work.resolve(), args.sdk.resolve()
    if args.operation in ("key", "build"):
        descriptor = context(source, work)
        if args.operation == "key":
            key = cache_key(descriptor)
            print(key)
            if args.github_output:
                with args.github_output.open("a", encoding="utf-8") as stream:
                    stream.write(f"cache-key={key}\n")
        else:
            build_sdk(source, work, directory, descriptor, args.jobs)
    elif args.operation == "unpack":
        if directory.exists() and any(directory.iterdir()):
            raise RuntimeError("SDK extraction requires an empty destination")
        directory.mkdir(parents=True, exist_ok=True)
        with tarfile.open(args.archive, "r:gz") as archive:
            archive.extractall(directory, filter="data")
        verify_sdk(directory, context(source, work))
    elif args.operation == "pack":
        verify_sdk(directory)
        args.archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(args.archive, "w:gz", dereference=False) as archive:
            for item in sorted(directory.iterdir()):
                archive.add(item, arcname=item.name)
    else:
        verify_sdk(directory)


if __name__ == "__main__":
    main()
