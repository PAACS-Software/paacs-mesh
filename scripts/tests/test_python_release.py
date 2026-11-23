"""Release gates require every Python/platform wheel and independent OS consumer."""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_python_release as release
import check_pypi_release as index_release


class PythonReleaseTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="paacs-release-gates-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.sha = "a" * 40
        (self.root / "VERSION.txt").write_text("1.2.0-rc.2\n")
        (self.root / "cmake").mkdir()
        source = b"test matching upstream source"
        (self.root / "cmake/OCCTVersion.cmake").write_text(
            'set(PAACS_OCCT_VERSION "8.0.1")\nset(PAACS_OCCT_URL "https://example.invalid/occt")\n'
            f'set(PAACS_OCCT_SHA512 "{hashlib.sha512(source).hexdigest()}")\n')
        self.builders, self.consumers = self.root / "builders", self.root / "consumers"
        self.output = self.root / "release"
        for identifier, row in release.TARGETS.items():
            platform, artifact = row["platform"], row["artifact"]
            directory = self.builders / artifact
            directory.mkdir(parents=True)
            wheel = directory / f"paacs_mesh-1.2.0rc2-{row['python_tag']}-{row['python_tag']}-{row['platform_tag']}.whl"
            with zipfile.ZipFile(wheel, "w") as archive:
                archive.writestr("paacs_mesh-1.2.0rc2.dist-info/METADATA",
                    "Metadata-Version: 2.4\nName: paacs-mesh\nVersion: 1.2.0rc2\nRequires-Python: " + release.REQUIRES_PYTHON + "\n")
            (directory / "occt-8.0.1.tar.gz").write_bytes(source)
            (directory / "dependency-audit.txt").write_text("test audited dependencies\n")
            release.write(directory / "native-parity.json", {"status": "passed", "cases": [{}, {}, {}]})
            record = {"source_revision": {"commit": self.sha, "clean": True}, "platform": platform,
                "target_id": identifier, "python_version": row["python"], "python_tag": row["python_tag"],
                "host": self.host(platform),
                "wheel": wheel.name, "sha256": release.sha256(wheel),
                "installed_smoke": "passed with isolated Python", "qualification": "passed with isolated Python",
                "installed_tests": {"discovered": 10, "passed": 9, "skipped": 1},
                "numpy_validation": [{"version": v, "status": "passed", "tests": 9}
                                     for v in row["numpy_versions"]],
                "occt_source": {"file": "occt-8.0.1.tar.gz", "sha256": release.sha256(directory / "occt-8.0.1.tar.gz")},
                "repair": {"audit": "dependency-audit.txt", "audit_sha256": release.sha256(directory / "dependency-audit.txt")},
                "native_parity": {"file": "native-parity.json", "sha256": release.sha256(directory / "native-parity.json")}}
            release.write(directory / "wheel-validation.json", record)
            for consumer_row in release.CONSUMERS.values():
                if consumer_row["target_id"] != identifier:
                    continue
                consumer = self.consumers / consumer_row["consumer_artifact"]
                consumer.mkdir(parents=True)
                release.write(consumer / "consumer-validation.json", {
                    **record, "consumer_id": consumer_row["id"], "consumer_os": consumer_row["consumer_os"],
                    "host": self.host(platform, consumer_row.get("required_os_major")),
                })

    def host(self, platform, macos_major=None):
        if platform.startswith("windows-"):
            return {"system": "Windows", "machine": "ARM64" if platform == "windows-arm64" else "AMD64",
                    "release": "11", "version": "10.0.26100"}
        if platform == "macos-arm64":
            return {"system": "Darwin", "machine": "arm64", "release": "24.0.0", "version": "Darwin",
                    "macos_version": f"{macos_major or 15}.0"}
        return {"system": "Linux", "machine": "aarch64" if platform == "linux-arm64" else "x86_64",
                "release": "6.0", "version": "Linux", "libc": "glibc", "libc_version": "2.28"}

    def builder_directory(self, platform, python="3.12"):
        return self.builders / release.ARTIFACTS["cp" + python.replace(".", "") + "-" + platform]

    def consumer_receipt(self, platform, python="3.12", os=None):
        identifier = "cp" + python.replace(".", "") + "-" + platform
        row = next(row for row in release.CONSUMERS.values()
                   if row["target_id"] == identifier and (os is None or row["consumer_os"] == os))
        return self.consumers / row["consumer_artifact"] / "consumer-validation.json"

    def assemble(self):
        return release.assemble(self.root, self.builders, self.consumers, self.sha, self.output)

    def change_builder(self, change):
        path = self.builders / next(iter(release.ARTIFACTS.values())) / "wheel-validation.json"
        record = release.load(path)
        change(record)
        release.write(path, record)

    def assert_gate_refuses(self):
        with self.assertRaises((RuntimeError, FileNotFoundError)):
            self.assemble()
        self.assertFalse(self.output.exists(), "A failed gate left a publication set")

    def use_qualified_dual_linux_wheel(self, platform="linux-x86_64", python="3.12"):
        directory = self.builder_directory(platform, python)
        original = next(directory.glob("*.whl"))
        architecture = "x86_64" if platform == "linux-x86_64" else "aarch64"
        tag = "cp" + python.replace(".", "")
        wheel = original.with_name(
            f"paacs_mesh-1.2.0rc2-{tag}-{tag}-manylinux_2_27_{architecture}.manylinux_2_28_{architecture}.whl")
        original.rename(wheel)
        with zipfile.ZipFile(wheel, "a") as archive:
            archive.writestr("paacs_mesh-1.2.0rc2.dist-info/WHEEL",
                "Wheel-Version: 1.0\nRoot-Is-Purelib: false\n"
                f"Tag: {tag}-{tag}-manylinux_2_27_{architecture}\n"
                f"Tag: {tag}-{tag}-manylinux_2_28_{architecture}\n")
        for path in (directory / "wheel-validation.json",
                     self.consumer_receipt(platform, python)):
            record = release.load(path)
            record.update(wheel=wheel.name, sha256=release.sha256(wheel))
            release.write(path, record)
        return wheel

    def test_actual_reviewed_dual_linux_filename_passes_full_assembly(self):
        wheel = self.use_qualified_dual_linux_wheel()
        manifest = self.assemble()
        self.assertEqual(manifest["wheels"][wheel.name], release.sha256(wheel))
        self.assertEqual(len(manifest["wheels"]), 20)
        self.assertEqual(release.sha256(self.output / "dist" / wheel.name), release.sha256(wheel))
        release.verify(manifest, self.output / "dist", self.sha)

    def test_linux_tags_require_expected_policy_and_only_reviewed_architecture(self):
        original = next(self.builder_directory("linux-x86_64").glob("*.whl"))
        for tags in ("manylinux_2_27_x86_64", "manylinux_2_28_aarch64",
                     "manylinux_2_28_x86_64.manylinux_2_29_x86_64",
                     "manylinux_2_26_x86_64.manylinux_2_28_x86_64",
                     "manylinux_2_28_x86_64.manylinux_2_28_aarch64",
                     "manylinux_2_28_x86_64.win_amd64"):
            with self.subTest(tags=tags):
                candidate = self.root / f"paacs_mesh-1.2.0rc2-cp312-cp312-{tags}.whl"
                candidate.write_bytes(original.read_bytes())
                with self.assertRaises(RuntimeError):
                    release.wheel_metadata(candidate, release.Version("1.2.0rc2"), "linux-x86_64")

    def test_dual_linux_tags_do_not_relax_python_abi_name_version_or_build(self):
        original = next(self.builder_directory("linux-x86_64").glob("*.whl"))
        dual = "manylinux_2_27_x86_64.manylinux_2_28_x86_64"
        for prefix in ("paacs_mesh-1.2.0rc2-cp311-cp312", "paacs_mesh-1.2.0rc2-cp312-cp311",
                       "paacs_mesh-1.2.0rc2-cp311.cp312-cp312", "paacs_mesh-1.2.0rc2-cp312-abi3.cp312",
                       "paacs_mesh-1.2.0rc3-cp312-cp312", "other_project-1.2.0rc2-cp312-cp312",
                       "paacs_mesh-1.2.0rc2-1-cp312-cp312"):
            with self.subTest(prefix=prefix):
                candidate = self.root / f"{prefix}-{dual}.whl"
                candidate.write_bytes(original.read_bytes())
                with self.assertRaises(RuntimeError):
                    release.wheel_metadata(candidate, release.Version("1.2.0rc2"), "linux-x86_64")

    def test_windows_and_macos_platform_tags_remain_exact(self):
        for platform in ("windows-x86_64", "macos-arm64"):
            original = next(self.builder_directory(platform).glob("*.whl"))
            candidate = self.root / original.name.replace(".whl", ".manylinux_2_28_x86_64.whl")
            candidate.write_bytes(original.read_bytes())
            with self.subTest(platform=platform), self.assertRaises(RuntimeError):
                release.wheel_metadata(candidate, release.Version("1.2.0rc2"), platform)

    def test_dual_linux_wheel_keeps_builder_and_consumer_hash_gates(self):
        wheel = self.use_qualified_dual_linux_wheel()
        original = wheel.read_bytes()
        wheel.write_bytes(original + b"changed after qualification")
        self.assert_gate_refuses()
        wheel.write_bytes(original)
        path = self.consumer_receipt("linux-x86_64")
        record = release.load(path)
        record["sha256"] = "0" * 64
        release.write(path, record)
        self.assert_gate_refuses()

    def test_complete_qualified_wheels_and_matching_source_are_staged_separately(self):
        manifest = self.assemble()
        self.assertEqual(len(list((self.output / "dist").iterdir())), 20)
        self.assertEqual(len(manifest["qualified_consumers"]), 28)
        self.assertTrue(all(p.suffix == ".whl" for p in (self.output / "dist").iterdir()))
        self.assertTrue((self.output / "release-assets/occt-8.0.1.tar.gz").is_file())
        release.verify(manifest, self.output / "dist", self.sha)
        with self.assertRaises(RuntimeError):
            self.assemble()

    def test_release_tag_and_workflow_revision_must_match(self):
        current = {"commit": self.sha, "tree": "c" * 40, "clean": True}
        with patch.object(release, "revision", return_value=current), \
                patch.object(release.subprocess, "check_output", return_value=self.sha):
            self.assertEqual(release.preflight(self.root, "v1.2.0-rc.2", self.sha)["commit"], self.sha)
            with self.assertRaises(RuntimeError):
                release.preflight(self.root, "v1.2.0-rc.1", self.sha)
            with self.assertRaises(RuntimeError):
                release.preflight(self.root, "v1.2.0-rc.2", "b" * 40)

    def test_an_unqualified_or_missing_consumer_blocks_release(self):
        path = self.consumer_receipt("linux-x86_64")
        path.unlink()
        self.assert_gate_refuses()

    def test_wheel_change_after_qualification_blocks_release(self):
        directory = self.builders / next(iter(release.ARTIFACTS.values()))
        next(directory.glob("*.whl")).write_bytes(b"modified wheel")
        self.assert_gate_refuses()

    def test_different_or_dirty_source_blocks_release(self):
        self.change_builder(lambda r: r["source_revision"].update(commit="b" * 40))
        self.assert_gate_refuses()
        self.change_builder(lambda r: r["source_revision"].update(commit=self.sha, clean=False))
        self.assert_gate_refuses()

    def test_missing_numpy_qualification_blocks_release(self):
        self.change_builder(lambda r: r["numpy_validation"].pop())
        self.assert_gate_refuses()

    def test_every_new_python_and_architecture_has_its_own_receipt(self):
        self.assertEqual(len(release.TARGETS), 20)
        self.assertEqual(len(release.CONSUMERS), 28)
        for platform, python in (("windows-arm64", "3.15"), ("linux-arm64", "3.14")):
            path = self.builder_directory(platform, python) / "wheel-validation.json"
            record = release.load(path)
            path.unlink()
            with self.subTest(platform=platform, python=python):
                self.assert_gate_refuses()
            release.write(path, record)

    def test_cross_python_abi_receipt_blocks_release(self):
        self.change_builder(lambda r: r.update(python_version="3.13", python_tag="cp313"))
        self.assert_gate_refuses()

    def test_consumer_receipt_cannot_substitute_another_python_or_target(self):
        path = self.consumer_receipt("linux-arm64", "3.15")
        original = release.load(path)
        for fields in ({"target_id": "cp314-linux-arm64"}, {"python_tag": "cp314"},
                       {"consumer_id": "cp315-linux-x86_64-linuxglibc228"}):
            release.write(path, {**original, **fields})
            with self.subTest(fields=fields):
                self.assert_gate_refuses()

    def test_macos_consumer_os_requires_actual_host_major_not_only_claimed_name(self):
        path = self.consumer_receipt("macos-arm64", os="macos27")
        record = release.load(path)
        record["host"]["macos_version"] = "26.1"
        release.write(path, record)
        self.assert_gate_refuses()

    def test_macos15_consumer_cannot_stand_in_for_newer_os(self):
        path = self.consumer_receipt("macos-arm64", "3.15", os="macos26")
        record = release.load(self.consumer_receipt("macos-arm64", "3.15", os="macos15"))
        release.write(path, record)
        self.assert_gate_refuses()

    def test_arm_hosts_and_linux_libc_are_verified(self):
        for platform, field, value in (("windows-arm64", "machine", "AMD64"),
                                       ("windows-arm64", "version", "10.0.19045"),
                                       ("linux-arm64", "libc_version", "2.27"),
                                       ("linux-arm64", "machine", "x86_64")):
            path = self.builder_directory(platform) / "wheel-validation.json"
            record = release.load(path)
            changed = json.loads(json.dumps(record))
            changed["host"][field] = value
            if platform == "windows-arm64" and field == "version":
                changed["host"]["release"] = "10"
            release.write(path, changed)
            with self.subTest(platform=platform, field=field):
                self.assert_gate_refuses()
            release.write(path, record)

    def test_aarch64_dual_linux_tags_remain_supported(self):
        wheel = self.use_qualified_dual_linux_wheel("linux-arm64", "3.15")
        manifest = self.assemble()
        self.assertEqual(manifest["wheels"][wheel.name], release.sha256(wheel))
        release.verify(manifest, self.output / "dist", self.sha)

    def test_unexpected_or_duplicate_artifact_directory_blocks_release(self):
        (self.builders / "paacs-mesh-cp312-win32").mkdir()
        self.assert_gate_refuses()

    def test_publisher_rechecks_target_and_consumer_completeness(self):
        manifest = self.assemble()
        for key in ("targets", "qualified_consumers"):
            changed = json.loads(json.dumps(manifest))
            if key == "targets":
                changed[key].pop(next(iter(changed[key])))
            else:
                changed[key].pop()
            with self.subTest(key=key), self.assertRaises(RuntimeError):
                release.verify(changed, self.output / "dist", self.sha)

    def test_publisher_rejects_mixed_python_target_even_when_file_hashes_match(self):
        manifest = self.assemble()
        first, second = "cp312-linux-arm64", "cp313-linux-arm64"
        manifest["targets"][first]["wheel"] = manifest["targets"][second]["wheel"]
        with self.assertRaises(RuntimeError):
            release.verify(manifest, self.output / "dist", self.sha)

    def test_index_install_refuses_partial_matrix_before_network_or_environment(self):
        manifest = self.assemble()
        manifest["qualified_consumers"].pop()
        release.write(self.output / "release-manifest.json", manifest)
        with patch.object(sys, "argv", ["check_pypi_release.py", "--manifest",
                str(self.output / "release-manifest.json"), "--fixture", str(self.root / "cylinder.step"),
                "--output", str(self.root / "index-install")]), \
                patch.object(index_release, "urlopen") as network, \
                patch.object(index_release.venv, "EnvBuilder") as environment:
            with self.assertRaises(RuntimeError):
                index_release.main()
            network.assert_not_called()
            environment.assert_not_called()

    def test_empty_installed_test_result_blocks_release(self):
        self.change_builder(lambda r: r.update(installed_tests={"discovered": 0, "passed": 0, "skipped": 0}))
        self.assert_gate_refuses()

    def test_unqualified_python_metadata_blocks_release_even_with_new_checksum(self):
        directory = self.builders / next(iter(release.ARTIFACTS.values()))
        wheel = next(directory.glob("*.whl"))
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr("paacs_mesh-1.2.0rc2.dist-info/METADATA",
                "Name: paacs-mesh\nVersion: 1.2.0rc2\nRequires-Python: >=3.12\n")
        self.change_builder(lambda r: r.update(sha256=release.sha256(wheel)))
        self.assert_gate_refuses()

    def test_source_archive_change_blocks_release(self):
        directory = self.builders / next(iter(release.ARTIFACTS.values()))
        (directory / "occt-8.0.1.tar.gz").write_bytes(b"wrong source")
        self.assert_gate_refuses()

    def test_publication_rechecks_the_exact_file_set(self):
        manifest = self.assemble()
        (self.output / "dist/extra.txt").write_text("must not upload")
        with self.assertRaises(RuntimeError):
            release.verify(manifest, self.output / "dist", self.sha)

    def test_existing_public_source_is_never_silently_replaced(self):
        self.assemble()
        assets = self.output / "release-assets"
        view = {"tagName": "v1.2.0-rc.2", "isDraft": False,
                "assets": [{"name": p.name} for p in assets.iterdir()]}
        def download(args, **kwargs):
            self.assertEqual(args[2], "download")
            target = Path(args[args.index("--dir") + 1]) / args[args.index("--pattern") + 1]
            target.write_bytes(b"different existing public asset")
        with patch.object(release.subprocess, "check_output", return_value=json.dumps(view)), \
                patch.object(release.subprocess, "run", side_effect=download) as calls:
            with self.assertRaises(RuntimeError):
                release.publish_sources(assets, "v1.2.0-rc.2", self.sha)
            self.assertEqual(calls.call_count, 1)

    def test_index_install_resolves_relative_output_before_changing_child_cwd(self):
        manifest = self.assemble()
        fixture = self.root / "cylinder.step"
        fixture.write_bytes(b"owned fixture; meshing subprocess is mocked")
        response = {"urls": [{"filename": name, "packagetype": "bdist_wheel",
                              "digests": {"sha256": digest}}
                             for name, digest in manifest["wheels"].items()]}
        relative_output = Path("build/index-install")
        expected_output = (self.root / relative_output).resolve()
        expected_environment = expected_output / "environment"
        previous_cwd = Path.cwd()
        try:
            os.chdir(self.root)
            with patch.object(sys, "argv", ["check_pypi_release.py", "--manifest",
                    str(self.output / "release-manifest.json"), "--fixture", str(fixture),
                    "--output", str(relative_output)]), \
                    patch.object(index_release, "urlopen", return_value=io.BytesIO(json.dumps(response).encode())), \
                    patch.object(index_release.venv, "EnvBuilder") as builder, \
                    patch.object(index_release.subprocess, "run") as processes, \
                    patch.object(sys, "stdout", new=io.StringIO()):
                builder.return_value.create.side_effect = lambda path: Path(path).mkdir(parents=True)
                index_release.main()
                builder.return_value.create.assert_called_once_with(expected_environment)
                self.assertEqual(processes.call_count, 3)
                for call in processes.call_args_list:
                    self.assertTrue(Path(call.args[0][0]).is_absolute())
                    self.assertTrue(Path(call.kwargs["cwd"]).is_absolute())
                    self.assertEqual(Path(call.kwargs["cwd"]), expected_environment)
                install = processes.call_args_list[0].args[0]
                requirement = Path(install[install.index("-r") + 1])
                self.assertTrue(requirement.is_absolute())
                self.assertEqual(requirement, expected_environment / "release.txt")
                self.assertEqual(requirement.read_text().count("--hash=sha256:"), 20)
                self.assertEqual((expected_environment / "cylinder.step").read_bytes(), fixture.read_bytes())
                self.assertTrue((expected_output / "index-validation.json").is_file())
        finally:
            os.chdir(previous_cwd)


if __name__ == "__main__":
    unittest.main()
