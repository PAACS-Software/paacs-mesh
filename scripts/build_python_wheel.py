#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
"""Build, repair and qualify a STEP wheel with its bundled CAD runtime.

Supports the Python/platform targets in the compatibility manifest. Build isolation
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
import sysconfig
import urllib.request
import venv
import zipfile

from python_compatibility import PLATFORM_TAGS, numpy_versions, python_tag, target


CRT_NAMES = {
    "msvcp140.dll", "msvcp140_1.dll", "msvcp140_2.dll",
    "msvcp140_atomic_wait.dll", "msvcp140_codecvt_ids.dll", "concrt140.dll",
    "vcruntime140.dll", "vcruntime140_1.dll", "vcruntime140_threads.dll",
}

def python_configuration() -> dict:
    """Accept only the regular, release CPython ABI named by a target wheel."""
    if (sys.implementation.name != "cpython" or struct.calcsize("P") != 8
            or sysconfig.get_config_var("Py_GIL_DISABLED")
            or sysconfig.get_config_var("Py_DEBUG") or getattr(sys, "abiflags", "")):
        raise RuntimeError("Wheel qualification requires regular-GIL, release 64-bit CPython")
    version = f"{sys.version_info[0]}.{sys.version_info[1]}"
    try:
        tag = python_tag(version)
    except ValueError as error:
        raise RuntimeError(str(error)) from error
    return {"python_version": version, "python_tag": tag}


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
    """Identify the native platform independently of the Python minor version."""
    if struct.calcsize("P") != 8:
        raise RuntimeError("Supported wheels require a 64-bit process")
    architecture = platform.machine().lower()
    if sys.platform == "win32":
        if architecture in {"amd64", "x86_64"}:
            return "windows-x86_64"
        if architecture in {"arm64", "aarch64"}:
            return "windows-arm64"
    if sys.platform == "darwin" and architecture in {"arm64", "aarch64"}:
        if int(platform.mac_ver()[0].split(".")[0]) < 15:
            raise RuntimeError("The supported macOS wheel requires macOS 15 or newer")
        return "macos-arm64"
    if sys.platform == "linux" and architecture in {"amd64", "x86_64", "arm64", "aarch64"}:
        libc, release = platform.libc_ver()
        try:
            compatible = libc == "glibc" and tuple(map(int, release.split(".")[:2])) >= (2, 28)
        except ValueError:
            compatible = False
        if not compatible:
            raise RuntimeError("The supported Linux wheel requires glibc 2.28 or newer")
        if building and release != "2.28":
            raise RuntimeError("Build Linux wheels inside the pinned manylinux 2.28 image")
        return "linux-arm64" if architecture in {"arm64", "aarch64"} else "linux-x86_64"
    raise RuntimeError("Supported wheels require Windows x64/ARM64, macOS 15+ ARM64 or Linux x64/ARM64")


def host_details() -> dict:
    result = {"system": platform.system(), "release": platform.release(),
              "version": platform.version(), "machine": platform.machine()}
    if sys.platform == "darwin":
        result["macos_version"] = platform.mac_ver()[0]
    elif sys.platform == "linux":
        result["libc"], result["libc_version"] = platform.libc_ver()
    return result


def target_identity() -> dict:
    configuration = python_configuration()
    current_platform = host_platform()
    selected = target(current_platform, configuration["python_version"])
    return {"target_id": selected["id"], **configuration,
            "platform": current_platform, "host": host_details()}


def wheel_platform(wheel: Path) -> str:
    configuration = python_configuration()
    fields = wheel.name.removesuffix(".whl").rsplit("-", 3)
    if len(fields) != 4 or fields[1:3] != [configuration["python_tag"]] * 2:
        raise RuntimeError("Expected a wheel matching the running CPython interpreter and native ABI")
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


def _visual_studio_installations() -> list[dict]:
    program_files = Path(os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)"))
    vswhere = program_files / "Microsoft Visual Studio/Installer/vswhere.exe"
    if not vswhere.is_file():
        return []
    return json.loads(subprocess.check_output(
        [str(vswhere), "-all", "-products", "*", "-format", "json", "-utf8"],
        encoding="utf-8"))


def _numeric_version(value: str) -> tuple[int, ...]:
    return tuple(map(int, value.split("."))) if re.fullmatch(r"\d+(?:\.\d+)*", value) else ()


def visual_studio_toolchain(selected: str | None = None) -> dict:
    """Choose one installed instance and its default target-capable C++ toolset."""
    selected = selected or host_platform()
    architecture = "arm64" if selected == "windows-arm64" else "x64"
    capabilities = json.loads(subprocess.check_output(["cmake", "-E", "capabilities"], text=True))
    generators = {}
    for item in capabilities.get("generators", []):
        match = re.match(r"Visual Studio (\d+) ", item["name"])
        if match:
            generators[int(match.group(1))] = item["name"]
    requested = os.environ.get("CMAKE_GENERATOR")
    requested_instance = os.environ.get("CMAKE_GENERATOR_INSTANCE")
    requested_toolset = os.environ.get("PAACS_MSVC_TOOLSET", "")
    if requested_toolset and not re.fullmatch(r"\d+\.\d+(?:\.\d+)?", requested_toolset):
        raise RuntimeError("PAACS_MSVC_TOOLSET must name an installed numeric toolset version, for example 14.44")
    candidates = []
    for installation in _visual_studio_installations():
        version = _numeric_version(installation.get("installationVersion", ""))
        if not version:
            continue
        name = generators.get(version[0])
        if not name or (requested and name != requested):
            continue
        root = Path(installation["installationPath"])
        if requested_instance and root.resolve() != Path(requested_instance).resolve():
            continue
        tools = root / "VC/Tools/MSVC"
        defaults = root / "VC/Auxiliary/Build"
        default_file = defaults / "Microsoft.VCToolsVersion.default.txt"
        default_version = default_file.read_text(encoding="utf-8-sig").strip() if default_file.is_file() else None
        toolsets = [tools / default_version] if default_version and not requested_toolset else sorted(
            (entry for entry in tools.glob("*") if _numeric_version(entry.name)),
            key=lambda entry: _numeric_version(entry.name), reverse=True)
        if requested_toolset:
            prefix = _numeric_version(requested_toolset)
            toolsets = [entry for entry in toolsets if _numeric_version(entry.name)[:len(prefix)] == prefix]
        compiler = None
        for toolset in toolsets:
            # Native ARM64 hosts and x64-hosted cross compilers are both valid.
            available = [entry for entry in toolset.glob(f"bin/Host*/{architecture}/cl.exe") if entry.is_file()]
            if available:
                native_host = "hostarm64" if platform.machine().lower() in {"arm64", "aarch64"} else "hostx64"
                compiler = min(available, key=lambda entry: (entry.parent.parent.name.lower() != native_host,
                                                            entry.as_posix().lower()))
                break
        if compiler is None:
            continue
        actual_tools_version = compiler.parents[3].name
        cmake_toolset = None
        if requested_toolset:
            families = []
            for entry in defaults.glob("Microsoft.VCToolsVersion.v*.default.*"):
                match = re.fullmatch(r"Microsoft\.VCToolsVersion\.(v\d+)\.default\.(?:txt|props)", entry.name)
                if not match:
                    continue
                contents = entry.read_text(encoding="utf-8-sig")
                versions = re.findall(r"\b\d+\.\d+(?:\.\d+)+\b", contents)
                if any(_numeric_version(value)[:2] == _numeric_version(actual_tools_version)[:2] for value in versions):
                    families.append(match.group(1))
            if not families:
                continue
            cmake_toolset = f"{max(set(families))},version={requested_toolset}"
        host_architecture = compiler.parent.parent.name[4:].lower()
        if requested_toolset or host_architecture != ("arm64" if native_host == "hostarm64" else "x64"):
            host_option = "host=" + ("ARM64" if host_architecture == "arm64" else host_architecture)
            cmake_toolset = ",".join(value for value in (cmake_toolset, host_option) if value)
        redist_file = defaults / "Microsoft.VCRedistVersion.default.txt"
        redist_version = redist_file.read_text(encoding="utf-8-sig").strip() if redist_file.is_file() else None
        if requested_toolset:
            compatible = sorted((entry.name for entry in (root / "VC/Redist/MSVC").glob("*")
                                 if _numeric_version(entry.name)[:2] == _numeric_version(actual_tools_version)[:2]),
                                key=_numeric_version, reverse=True)
            redist_version = compatible[0] if compatible else redist_version
            if redist_version and _numeric_version(redist_version)[:2] < _numeric_version(actual_tools_version)[:2]:
                redist_version = None
        record = {"installation": str(root.resolve()), "generator": name, "architecture": architecture,
                  "compiler": str(compiler.resolve()), "tools_version": actual_tools_version,
                  "host_architecture": host_architecture, "requested_toolset": requested_toolset or None,
                  "toolset": cmake_toolset, "redist_version": redist_version}
        candidates.append((version, str(root).lower(), record))
    if not candidates:
        raise RuntimeError("No installed Visual Studio instance has a supported generator and target C++ compiler")
    return max(candidates, key=lambda item: item[:2])[2]


def cmake_generator(platform: str | None = None) -> list[str]:
    """Pin the same installed toolchain used for compiler and CRT selection."""
    selected = platform or host_platform()
    if not selected.startswith("windows-"):
        return ["-G", "Ninja"]
    selected_toolchain = visual_studio_toolchain(selected)
    result = ["-G", selected_toolchain["generator"], "-A",
              "ARM64" if selected == "windows-arm64" else "x64",
              "-DCMAKE_GENERATOR_INSTANCE:INTERNAL=" + selected_toolchain["installation"]]
    if selected_toolchain["toolset"]:
        result.extend(["-T", selected_toolchain["toolset"]])
    return result


def binary_architecture(payload: bytes) -> tuple[str, str]:
    """Inspect binary headers without executing a bundled dependency."""
    if payload.startswith(b"MZ") and len(payload) >= 64:
        offset = struct.unpack_from("<I", payload, 0x3c)[0]
        if offset + 6 <= len(payload) and payload[offset:offset + 4] == b"PE\0\0":
            machine = struct.unpack_from("<H", payload, offset + 4)[0]
            return "pe", {0x8664: "x86_64", 0xaa64: "arm64"}.get(machine, "unsupported")
    if payload.startswith(b"\x7fELF") and len(payload) >= 20:
        if payload[4] != 2 or payload[5] not in (1, 2):
            raise RuntimeError("Expected a 64-bit ELF image with known byte order")
        machine = struct.unpack_from("<H" if payload[5] == 1 else ">H", payload, 18)[0]
        return "elf", {62: "x86_64", 183: "arm64"}.get(machine, "unsupported")
    endian = {b"\xcf\xfa\xed\xfe": "<", b"\xfe\xed\xfa\xcf": ">"}.get(payload[:4])
    if endian and len(payload) >= 12:
        machine = struct.unpack_from(endian + "I", payload, 4)[0]
        return "mach-o", {0x01000007: "x86_64", 0x0100000c: "arm64"}.get(machine, "unsupported")
    raise RuntimeError("Expected a PE, ELF64 or thin Mach-O64 binary")


def validate_binary(payload: bytes, name: str, selected: str | None = None) -> None:
    selected = selected or host_platform()
    expected_format = {"windows": "pe", "linux": "elf", "macos": "mach-o"}[selected.split("-")[0]]
    expected_arch = "x86_64" if selected.endswith("x86_64") else "arm64"
    try:
        actual = binary_architecture(payload)
    except RuntimeError as error:
        raise RuntimeError(f"Invalid bundled binary {name}: {error}") from error
    if actual != (expected_format, expected_arch):
        raise RuntimeError(f"Binary {name} has {actual}, expected {(expected_format, expected_arch)}")


def runtime_directory(explicit: Path | None) -> Path:
    selected = host_platform()
    architecture = "arm64" if selected == "windows-arm64" else "x64"
    candidates = [explicit] if explicit else []
    if not explicit:
        toolchain = visual_studio_toolchain(selected)
        redist = Path(toolchain["installation"]) / "VC/Redist/MSVC"
        default_version = toolchain["redist_version"]
        if default_version:
            candidates.extend(sorted((redist / default_version).glob(f"{architecture}/Microsoft.VC*.CRT")))
        else:
            toolset = _numeric_version(toolchain["tools_version"])
            versions = sorted((entry for entry in redist.glob("*")
                               if _numeric_version(entry.name)[:2] == toolset[:2]),
                              key=lambda entry: _numeric_version(entry.name), reverse=True)
            for version in versions:
                candidates.extend(sorted(version.glob(f"{architecture}/Microsoft.VC*.CRT")))
    for candidate in candidates:
        if candidate is not None and (candidate / "msvcp140.dll").is_file():
            validate_binary((candidate / "msvcp140.dll").read_bytes(), str(candidate / "msvcp140.dll"), selected)
            return candidate.resolve()
    raise RuntimeError(f"Visual C++ {architecture} redistributable files not found; set --vc-runtime-dir")


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


def _binary_crt_imports(payload: bytes, name: str) -> set[str]:
    import pefile  # Installed by delvewheel; inspect actual binary imports.
    image = pefile.PE(data=payload, fast_load=True)
    try:
        directories = ("IMAGE_DIRECTORY_ENTRY_IMPORT", "IMAGE_DIRECTORY_ENTRY_DELAY_IMPORT")
        image.parse_data_directories(directories=[pefile.DIRECTORY_ENTRY[key] for key in directories])
        imports = [entry.dll.decode("ascii").lower() for attribute in
                   ("DIRECTORY_ENTRY_IMPORT", "DIRECTORY_ENTRY_DELAY_IMPORT")
                   for entry in getattr(image, attribute, ())]
        return set(imports) & CRT_NAMES
    finally:
        image.close()


def crt_dependency_closure(payloads: list[tuple[str, bytes]], runtime: Path) -> list[str]:
    """Follow only imported CRT files, validating every required dependency's CPU."""
    required, pending = set(), list(payloads)
    while pending:
        name, payload = pending.pop()
        validate_binary(payload, name)
        for dependency in _binary_crt_imports(payload, name):
            if dependency in required:
                continue
            binary = runtime / dependency
            if not binary.is_file():
                raise RuntimeError(f"Required redistributable is missing: {binary}")
            required.add(dependency)
            pending.append((str(binary), binary.read_bytes()))
    return sorted(required)


def crt_dependencies(wheel: Path, sdk: Path, runtime: Path) -> list[str]:
    with zipfile.ZipFile(wheel) as archive:
        payloads = [(name, archive.read(name)) for name in archive.namelist() if name.endswith(".pyd")]
    payloads.extend((str(binary), binary.read_bytes()) for binary in (sdk / "bin").glob("*.dll"))
    return crt_dependency_closure(payloads, runtime)


def validate_contents(wheel: Path, runtime: Path | None = None,
                      crt: list[str] | tuple = ()) -> dict:
    current_platform = host_platform()
    if PLATFORM_TAGS[current_platform] not in wheel_platform(wheel).split("."):
        raise RuntimeError("The wheel platform tag does not match the supported host")
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        extension = ".pyd" if sys.platform == "win32" else ".so"
        extensions = [name for name in names if name.startswith("paacs_mesh/_native") and name.endswith(extension)]
        if len(extensions) != 1:
            raise RuntimeError("The wheel must contain exactly one native STEP binding")
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
        for name in [*extensions, *libraries]:
            validate_binary(archive.read(name), name, current_platform)
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
        result = {"wheel": wheel.name, "sha256": sha256(wheel), **target_identity(),
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
        command = [sys.executable, "-m", "auditwheel", "repair", "--plat", PLATFORM_TAGS[host_platform()],
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
    payloads = [(str(native_cli), native_cli.read_bytes())]
    payloads.extend((str(binary), binary.read_bytes()) for binary in (sdk / "bin").glob("*.dll"))
    for name in crt_dependency_closure(payloads, runtime):
        shutil.copyfile(runtime / name, native_cli.parent / name)


def installed_smoke(source: Path, wheel: Path, build_dir: Path) -> list[dict]:
    configuration = python_configuration()
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
        for version in numpy_versions(host_platform(), configuration["python_version"]):
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
    parser.add_argument("--vc-runtime-dir", type=Path, help="Licensed Microsoft CRT directory matching this architecture")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--no-isolation", action="store_true", help="Use already installed pinned build requirements")
    args = parser.parse_args()
    try:
        python_configuration()
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
        generator = cmake_generator(current_platform)
        env.update(CMAKE_GENERATOR=generator[1], CMAKE_GENERATOR_PLATFORM=generator[3])
        env.pop("CMAKE_GENERATOR_TOOLSET", None)
        settings = iter(generator[4:])
        for setting in settings:
            if setting == "-T":
                value = next(settings)
                env["CMAKE_GENERATOR_TOOLSET"] = value
                command.append(f"-Ccmake.define.CMAKE_GENERATOR_TOOLSET={value}")
            else:
                key, value = setting.removeprefix("-D").split("=", 1)
                command.append(f"-Ccmake.define.{key.split(':')[0]}={value}")
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
    if sys.platform == "win32":
        receipt["toolchain"] = visual_studio_toolchain(current_platform)
    build_info = build_dir / "cmake/build-info-Release.txt"
    if not build_info.is_file():
        raise RuntimeError("The native build omitted its compiler and configuration record")
    build_info_output = output / "native-build-info.txt"
    shutil.copyfile(build_info, build_info_output)
    receipt["build_info"] = {"file": build_info_output.name, "sha256": sha256(build_info_output)}
    if current_platform.startswith("linux-"):
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
