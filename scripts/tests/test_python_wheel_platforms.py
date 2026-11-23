"""Fail-closed interpreter, binary and independent-consumer qualification checks."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
import os
from pathlib import Path
import struct
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_python_wheel as build
import check_python_wheel as check
import python_compatibility as compatibility


def image(kind, architecture):
    if kind == "pe":
        payload = bytearray(128)
        payload[:2] = b"MZ"
        struct.pack_into("<I", payload, 0x3c, 64)
        payload[64:68] = b"PE\0\0"
        struct.pack_into("<H", payload, 68, 0x8664 if architecture == "x86_64" else 0xaa64)
    elif kind == "elf":
        payload = bytearray(64)
        payload[:6] = b"\x7fELF\x02\x01"
        struct.pack_into("<H", payload, 18, 62 if architecture == "x86_64" else 183)
    else:
        payload = bytearray(32)
        payload[:4] = b"\xcf\xfa\xed\xfe"
        struct.pack_into("<I", payload, 4, 0x01000007 if architecture == "x86_64" else 0x0100000c)
    return bytes(payload)


def identity(platform, version="3.12", *, macos="15.7", windows="10.0.26100"):
    selected = compatibility.target(platform, version)
    system = {"windows": "Windows", "macos": "Darwin", "linux": "Linux"}[platform.split("-")[0]]
    host = {"system": system, "release": "11" if system == "Windows" else "test",
            "version": windows, "machine": "x86_64" if platform.endswith("x86_64") else "arm64"}
    if system == "Darwin":
        host["macos_version"] = macos
    return {"target_id": selected["id"], "python_version": version,
            "python_tag": selected["python_tag"], "platform": platform, "host": host}


def installation(root, version, toolset, targets=("x64", "arm64")):
    defaults = root / "VC/Auxiliary/Build"
    defaults.mkdir(parents=True)
    for component in ("VCTools", "VCRedist"):
        (defaults / f"Microsoft.{component}Version.default.txt").write_text(toolset, encoding="utf-8")
    for target in targets:
        for host in ("Hostx64", "Hostarm64"):
            compiler = root / f"VC/Tools/MSVC/{toolset}/bin/{host}/{target}/cl.exe"
            compiler.parent.mkdir(parents=True)
            compiler.write_bytes(b"compiler")
    return {"installationVersion": version, "installationPath": str(root)}


class InterpreterPlatformTests(unittest.TestCase):
    def configuration(self, minor, flags=None, abiflags=""):
        with (patch.object(build.sys, "version_info", (3, minor, 0)),
              patch.object(build.sys, "abiflags", abiflags, create=True),
              patch.object(build.sysconfig, "get_config_var", side_effect=lambda key: (flags or {}).get(key, 0))):
            return build.python_configuration()

    def test_regular_release_interpreters_match_their_exact_native_abi(self):
        for minor in (12, 13, 14, 15):
            with self.subTest(minor=minor):
                self.assertEqual(self.configuration(minor),
                                 {"python_version": f"3.{minor}", "python_tag": f"cp3{minor}"})
        for flags, abi in (({"Py_GIL_DISABLED": 1}, ""), ({"Py_DEBUG": 1}, ""), ({}, "t"), ({}, "d")):
            with self.subTest(flags=flags, abi=abi), self.assertRaises(RuntimeError):
                self.configuration(13, flags, abi)
        with self.assertRaises(RuntimeError):
            self.configuration(11)
        with patch.object(build.sys, "implementation", SimpleNamespace(name="pypy")), self.assertRaises(RuntimeError):
            build.python_configuration()

    def test_platform_identification_does_not_depend_on_python_version(self):
        for machine, expected in (("AMD64", "windows-x86_64"), ("ARM64", "windows-arm64")):
            with (self.subTest(machine=machine), patch.object(build.sys, "platform", "win32"),
                  patch.object(build.platform, "machine", return_value=machine),
                  patch.object(build.sys, "version_info", (3, 11, 0))):
                self.assertEqual(build.host_platform(), expected)

    def test_linux_native_architecture_and_floor_are_preserved(self):
        for machine, expected in (("x86_64", "linux-x86_64"), ("aarch64", "linux-arm64")):
            with (self.subTest(machine=machine), patch.object(build.sys, "platform", "linux"),
                  patch.object(build.platform, "machine", return_value=machine),
                  patch.object(build.platform, "libc_ver", return_value=("glibc", "2.28"))):
                self.assertEqual(build.host_platform(building=True), expected)
        for libc, release in (("musl", "1.2"), ("glibc", "2.27"), ("glibc", "unknown")):
            with (self.subTest(libc=libc, release=release), patch.object(build.sys, "platform", "linux"),
                  patch.object(build.platform, "machine", return_value="aarch64"),
                  patch.object(build.platform, "libc_ver", return_value=(libc, release)), self.assertRaises(RuntimeError)):
                build.host_platform()
        with (patch.object(build.sys, "platform", "linux"), patch.object(build.platform, "machine", return_value="aarch64"),
              patch.object(build.platform, "libc_ver", return_value=("glibc", "2.40"))):
            self.assertEqual(build.host_platform(), "linux-arm64")
            with self.assertRaises(RuntimeError):
                build.host_platform(building=True)

    def test_macos_keeps_the_same_arm64_floor(self):
        for version in ("15.7", "26.1", "27.0"):
            with (self.subTest(version=version), patch.object(build.sys, "platform", "darwin"),
                  patch.object(build.platform, "machine", return_value="arm64"),
                  patch.object(build.platform, "mac_ver", return_value=(version, (), ""))):
                self.assertEqual(build.host_platform(), "macos-arm64")
        with (patch.object(build.sys, "platform", "darwin"), patch.object(build.platform, "machine", return_value="arm64"),
              patch.object(build.platform, "mac_ver", return_value=("14.7", (), "")), self.assertRaises(RuntimeError)):
            build.host_platform()

    def test_wheel_abi_must_match_running_interpreter(self):
        config = {"python_version": "3.14", "python_tag": "cp314"}
        with patch.object(build, "python_configuration", return_value=config):
            self.assertEqual(build.wheel_platform(Path("paacs_mesh-1.2.0-cp314-cp314-win_arm64.whl")), "win_arm64")
            for tag in ("cp312-cp312", "cp314-abi3", "cp314-cp314t"):
                with self.subTest(tag=tag), self.assertRaises(RuntimeError):
                    build.wheel_platform(Path(f"paacs_mesh-1.2.0-{tag}-win_arm64.whl"))


class BinaryQualificationTests(unittest.TestCase):
    def test_pe_elf_and_macho_headers_are_checked_for_every_native_architecture(self):
        for selected in compatibility.PLATFORM_TAGS:
            kind = {"windows": "pe", "macos": "mach-o", "linux": "elf"}[selected.split("-")[0]]
            arch = "x86_64" if selected.endswith("x86_64") else "arm64"
            with self.subTest(platform=selected):
                self.assertEqual(build.binary_architecture(image(kind, arch)), (kind, arch))
                build.validate_binary(image(kind, arch), "example", selected)
                with self.assertRaisesRegex(RuntimeError, "expected"):
                    build.validate_binary(image(kind, "arm64" if arch == "x86_64" else "x86_64"), "example", selected)
        for payload in (b"", b"MZ", b"\xca\xfe\xba\xbe" + b"\0" * 40,
                        b"\x7fELF\x01\x01" + b"\0" * 40):
            with self.subTest(payload=payload[:6]), self.assertRaises(RuntimeError):
                build.binary_architecture(payload)

    def synthetic_wheel(self, directory, row, *, wrong_arch=False):
        selected = row["platform"]
        is_windows = selected.startswith("windows-")
        kind = "pe" if is_windows else "mach-o" if selected == "macos-arm64" else "elf"
        arch = "x86_64" if selected.endswith("x86_64") else "arm64"
        wheel = directory / f"paacs_mesh-1.2.0-{row['python_tag']}-{row['python_tag']}-{row['platform_tag']}.whl"
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr("paacs_mesh/_native.test" + (".pyd" if is_windows else ".so"), image(kind, arch))
            for name in ("TKernel", "TKDESTEP"):
                library = name + ".dll" if is_windows else "lib" + name + (".dylib" if kind == "mach-o" else ".so.8.0")
                library_arch = ("arm64" if arch == "x86_64" else "x86_64") if wrong_arch else arch
                archive.writestr("paacs_mesh.libs/" + library, image(kind, library_arch))
            for name in ("LICENSE_LGPL_21.txt", "OCCT_LGPL_EXCEPTION.txt", "OCCT_NOTICE.txt", "pybind11-LICENSE.txt", "LICENSE", "NOTICE"):
                archive.writestr("paacs_mesh/_occt/" + name, "notice")
            for name in ("BOPAlgo", "SHMessage", "XSMessage", "XSTEPResource", "UnitsAPI", "StdResource"):
                archive.writestr(f"paacs_mesh/_occt/resources/{name}/example", "resource")
        return wheel

    def test_all_twenty_wheel_targets_validate_dependency_architecture_and_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            for row in compatibility.targets():
                selected = row["platform"]
                sys_platform = "win32" if selected.startswith("windows-") else "darwin" if selected == "macos-arm64" else "linux"
                config = {"python_version": row["python"], "python_tag": row["python_tag"]}
                with (self.subTest(target=row["id"]), patch.object(build, "host_platform", return_value=selected),
                      patch.object(build.sys, "platform", sys_platform),
                      patch.object(build, "python_configuration", return_value=config),
                      patch.object(build, "host_details", return_value={"machine": "test"})):
                    wheel = self.synthetic_wheel(Path(temporary), row)
                    receipt = build.validate_contents(wheel)
                    self.assertEqual(receipt["target_id"], row["id"])
                    self.assertEqual(receipt["python_version"], row["python"])
                    self.synthetic_wheel(Path(temporary), row, wrong_arch=True)
                    with self.assertRaisesRegex(RuntimeError, "Binary.*expected"):
                        build.validate_contents(wheel)


class ToolchainAndConsumerTests(unittest.TestCase):
    def test_visual_studio_generation_is_discovered_instead_of_hardcoded(self):
        capabilities = {"generators": [{"name": "Visual Studio 17 2022"}, {"name": "Visual Studio 18 2026"}]}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            installed = [installation(root / "vs17", "17.14.1", "14.44"),
                         installation(root / "vs18", "18.1.2", "14.50", targets=("x64",))]
            with (patch.object(build, "_visual_studio_installations", return_value=installed),
                  patch.object(build.subprocess, "check_output", return_value=json.dumps(capabilities)),
                  patch.dict(os.environ, {}, clear=True)):
                # The newest IDE lacks the ARM64 compiler; use the capable one.
                self.assertEqual(build.cmake_generator("windows-arm64")[:4], ["-G", "Visual Studio 17 2022", "-A", "ARM64"])
                selected = build.cmake_generator("windows-x86_64")
                self.assertEqual(selected[:4], ["-G", "Visual Studio 18 2026", "-A", "x64"])
                self.assertIn(str((root / "vs18").resolve()), selected[4])
                with patch.dict(os.environ, {"CMAKE_GENERATOR": "Visual Studio 17 2022"}):
                    self.assertEqual(build.cmake_generator("windows-arm64")[1], "Visual Studio 17 2022")
                with patch.object(build.platform, "machine", return_value="ARM64"):
                    self.assertIn("Hostarm64", build.visual_studio_toolchain("windows-arm64")["compiler"])
        self.assertEqual(build.cmake_generator("linux-arm64"), ["-G", "Ninja"])

    def test_crt_search_matches_architecture_without_fixing_a_vc_generation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            runtime = root / "VC/Redist/MSVC/14.50/arm64/Microsoft.VC145.CRT"
            runtime.mkdir(parents=True)
            (runtime / "msvcp140.dll").write_bytes(image("pe", "arm64"))
            with (patch.object(build, "host_platform", return_value="windows-arm64"),
                  patch.object(build, "visual_studio_toolchain", return_value={"installation": str(root),
                      "redist_version": None, "tools_version": "14.50"}),
                  patch.dict(os.environ, {}, clear=True)):
                self.assertEqual(build.runtime_directory(None), runtime.resolve())
                (runtime / "msvcp140.dll").write_bytes(image("pe", "x86_64"))
                with self.assertRaisesRegex(RuntimeError, "expected"):
                    build.runtime_directory(runtime)

    def test_generator_instance_and_crt_use_same_default_toolset_despite_vswhere_order(self):
        capabilities = {"generators": [{"name": "Visual Studio 17 2022"}]}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            installed = [installation(root / "old", "17.9", "14.9"), installation(root / "new", "17.10", "14.10")]
            runtime = root / "new/VC/Redist/MSVC/14.10/x64/Microsoft.VC145.CRT"
            runtime.mkdir(parents=True)
            (runtime / "msvcp140.dll").write_bytes(image("pe", "x86_64"))
            with (patch.object(build, "host_platform", return_value="windows-x86_64"),
                  patch.object(build, "_visual_studio_installations", return_value=installed),
                  patch.object(build.subprocess, "check_output", return_value=json.dumps(capabilities)),
                  patch.dict(os.environ, {"VCToolsRedistDir": str(root / "old/VC/Redist/MSVC/14.9")}, clear=True)):
                self.assertIn(str((root / "new").resolve()), build.cmake_generator()[4])
                self.assertEqual(build.runtime_directory(None), runtime.resolve())
                (root / "new/VC/Auxiliary/Build/Microsoft.VCToolsVersion.default.txt").write_text("14.99", encoding="utf-8")
                self.assertIn(str((root / "old").resolve()), build.cmake_generator()[4])

    def test_explicit_toolset_selects_installed_older_native_compiler_and_matching_crt(self):
        capabilities = {"generators": [{"name": "Visual Studio 17 2022"}, {"name": "Visual Studio 18 2026"}]}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            installed = [installation(root / "vs17", "17.14", "14.44.35207"),
                         installation(root / "vs18", "18.12", "14.52.37025")]
            selected_root = root / "vs18"
            older = selected_root / "VC/Tools/MSVC/14.44.35207/bin/Hostarm64/arm64/cl.exe"
            older.parent.mkdir(parents=True)
            older.write_bytes(b"older native compiler")
            (selected_root / "VC/Auxiliary/Build/Microsoft.VCToolsVersion.v143.default.txt").write_text("14.44.35207", encoding="utf-8")
            for version in ("14.44.35112", "14.52.37025"):
                runtime = selected_root / f"VC/Redist/MSVC/{version}/arm64/Microsoft.VC143.CRT"
                runtime.mkdir(parents=True)
                (runtime / "msvcp140.dll").write_bytes(image("pe", "arm64"))
            with (patch.object(build, "host_platform", return_value="windows-arm64"),
                  patch.object(build.platform, "machine", return_value="ARM64"),
                  patch.object(build, "_visual_studio_installations", return_value=installed),
                  patch.object(build.subprocess, "check_output", return_value=json.dumps(capabilities)),
                  patch.dict(os.environ, {}, clear=True)):
                self.assertEqual(build.visual_studio_toolchain()["tools_version"], "14.52.37025")
                with patch.dict(os.environ, {"PAACS_MSVC_TOOLSET": "14.44"}):
                    actual = build.visual_studio_toolchain()
                    self.assertEqual(actual["installation"], str(selected_root.resolve()))
                    self.assertEqual(actual["compiler"], str(older.resolve()))
                    self.assertEqual(actual["requested_toolset"], "14.44")
                    self.assertEqual(actual["host_architecture"], "arm64")
                    self.assertEqual(actual["redist_version"], "14.44.35112")
                    arguments = build.cmake_generator()
                    self.assertEqual(arguments[-2:], ["-T", "v143,version=14.44,host=ARM64"])
                    self.assertEqual(build.runtime_directory(None),
                        (selected_root / "VC/Redist/MSVC/14.44.35112/arm64/Microsoft.VC143.CRT").resolve())
                with patch.dict(os.environ, {"PAACS_MSVC_TOOLSET": "14.43"}), self.assertRaises(RuntimeError):
                    build.visual_studio_toolchain()
                with patch.dict(os.environ, {"PAACS_MSVC_TOOLSET": "invalid"}), self.assertRaisesRegex(RuntimeError, "numeric"):
                    build.visual_studio_toolchain()

    def test_requested_arm_toolset_can_use_matching_x64_host_cross_compiler(self):
        capabilities = {"generators": [{"name": "Visual Studio 18 2026"}]}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            installed = installation(root, "18.12", "14.44.35207")
            (root / "VC/Tools/MSVC/14.44.35207/bin/Hostarm64/arm64/cl.exe").unlink()
            (root / "VC/Auxiliary/Build/Microsoft.VCToolsVersion.v143.default.txt").write_text("14.44.35207", encoding="utf-8")
            with (patch.object(build.platform, "machine", return_value="ARM64"),
                  patch.object(build, "_visual_studio_installations", return_value=[installed]),
                  patch.object(build.subprocess, "check_output", return_value=json.dumps(capabilities)),
                  patch.dict(os.environ, {"PAACS_MSVC_TOOLSET": "14.44"}, clear=True)):
                self.assertEqual(build.cmake_generator("windows-arm64")[-1], "v143,version=14.44,host=x64")

    def test_wheel_backend_receives_the_same_cmake_toolset_as_the_sdk(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "sdk-config"
            config.mkdir()
            (config / "OpenCASCADEConfig.cmake").write_text("", encoding="utf-8")
            toolset = "v143,version=14.44,host=ARM64"
            generator = ["-G", "Visual Studio 18 2026", "-A", "ARM64",
                         "-DCMAKE_GENERATOR_INSTANCE:INTERNAL=" + str(root), "-T", toolset]
            arguments = ["builder", "--build-dir", str(root / "build"), "--out-dir", str(root / "dist"),
                         "--occt-dir", str(config)]
            with (patch.object(build.sys, "argv", arguments), patch.object(build.sys, "platform", "win32"),
                  patch.object(build, "host_platform", return_value="windows-arm64"),
                  patch.object(build, "python_configuration", return_value={
                      "python_version": "3.12", "python_tag": "cp312"}),
                  patch.object(build, "cmake_generator", return_value=generator),
                  patch.object(build, "run", side_effect=RuntimeError("stop before compilation")) as call):
                with self.assertRaisesRegex(RuntimeError, "stop before compilation"):
                    build.main()
                command = call.call_args.args[0]
                environment = call.call_args.kwargs["env"]
                self.assertEqual(environment["CMAKE_GENERATOR_TOOLSET"], toolset)
                self.assertIn("-Ccmake.define.CMAKE_GENERATOR_TOOLSET=" + toolset, command)
                self.assertIn("-Ccmake.define.CMAKE_GENERATOR_INSTANCE=" + str(root), command)

    def test_native_cli_copies_only_its_native_crt_dependency_closure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "cli/paacs-cad-mesh.exe"
            executable.parent.mkdir()
            executable.write_bytes(image("pe", "arm64"))
            sdk = root / "sdk"
            (sdk / "bin").mkdir(parents=True)
            (sdk / "bin/TKernel.dll").write_bytes(image("pe", "arm64"))
            runtime = root / "arm64/Microsoft.VC143.CRT"
            runtime.mkdir(parents=True)
            for name in ("msvcp140.dll", "vcruntime140.dll"):
                (runtime / name).write_bytes(image("pe", "arm64"))
            # Real ARM64 redistributable directories can also contain this x64
            # DLL. Its presence is harmless when no ARM64 binary imports it.
            (runtime / "vcruntime140_1.dll").write_bytes(image("pe", "x86_64"))
            graph = {"paacs-cad-mesh.exe": {"msvcp140.dll"}, "TKernel.dll": {"msvcp140.dll"},
                     "msvcp140.dll": {"vcruntime140.dll"}, "vcruntime140.dll": set()}
            with (patch.object(build.sys, "platform", "win32"),
                  patch.object(build, "host_platform", return_value="windows-arm64"),
                  patch.object(build, "_binary_crt_imports", side_effect=lambda payload, name: graph[Path(name).name])):
                build.prepare_native_runtime(executable, sdk, runtime)
                self.assertTrue((executable.parent / "msvcp140.dll").is_file())
                self.assertTrue((executable.parent / "vcruntime140.dll").is_file())
                self.assertFalse((executable.parent / "vcruntime140_1.dll").exists())
                # A transitive requirement is never skipped or accepted with
                # the wrong CPU just because the parent dependency is native.
                graph["msvcp140.dll"].add("vcruntime140_1.dll")
                with self.assertRaisesRegex(RuntimeError, "vcruntime140_1.*x86_64.*arm64"):
                    build.prepare_native_runtime(executable, sdk, runtime)
                self.assertFalse((executable.parent / "vcruntime140_1.dll").exists())
                (runtime / "vcruntime140_1.dll").write_bytes(image("pe", "arm64"))
                graph["vcruntime140_1.dll"] = {"vcruntime140.dll"}
                build.prepare_native_runtime(executable, sdk, runtime)
                self.assertEqual((executable.parent / "vcruntime140_1.dll").read_bytes(), image("pe", "arm64"))

    def test_each_macos_consumer_requires_its_actual_os_major_and_same_wheel(self):
        for major in (15, 26, 27):
            selected = f"cp314-macos-arm64-macos{major}"
            with (self.subTest(major=major), patch.object(check, "target_identity",
                  return_value=identity("macos-arm64", "3.14", macos=f"{major}.1"))):
                actual = check.consumer_identity(selected)
                self.assertEqual(actual["consumer_id"], selected)
                self.assertEqual(actual["target_id"], "cp314-macos-arm64")
            with patch.object(check, "target_identity", return_value=identity("macos-arm64", "3.14", macos="14.7")):
                with self.assertRaisesRegex(RuntimeError, "requires macOS"):
                    check.consumer_identity(selected)

    def test_consumer_cannot_claim_another_interpreter_cpu_or_old_windows(self):
        with patch.object(check, "target_identity", return_value=identity("windows-x86_64", "3.12")):
            for selected in ("cp313-windows-x86_64-windows2022", "cp312-windows-arm64-windows11arm"):
                with self.subTest(selected=selected), self.assertRaisesRegex(RuntimeError, "does not match"):
                    check.consumer_identity(selected)
        selected = "cp315-windows-arm64-windows11arm"
        with patch.object(check, "target_identity", return_value=identity("windows-arm64", "3.15")):
            self.assertEqual(check.consumer_identity(selected)["consumer_os"], "windows11arm")
        with patch.object(check, "target_identity", return_value=identity("windows-arm64", "3.15", windows="10.0.19045")):
            with self.assertRaisesRegex(RuntimeError, "Windows 11"):
                check.consumer_identity(selected)


if __name__ == "__main__":
    unittest.main()
