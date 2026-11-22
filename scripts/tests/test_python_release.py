"""Release gates use receipts and hashes from all three platform consumers."""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import prepare_python_release as release


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
        for platform, artifact in release.ARTIFACTS.items():
            directory = self.builders / artifact
            directory.mkdir(parents=True)
            wheel = directory / f"paacs_mesh-1.2.0rc2-cp312-cp312-{release.PLATFORM_TAGS[platform]}.whl"
            with zipfile.ZipFile(wheel, "w") as archive:
                archive.writestr("paacs_mesh-1.2.0rc2.dist-info/METADATA",
                    "Metadata-Version: 2.4\nName: paacs-mesh\nVersion: 1.2.0rc2\nRequires-Python: >=3.12,<3.13\n")
            (directory / "occt-8.0.1.tar.gz").write_bytes(source)
            (directory / "dependency-audit.txt").write_text("test audited dependencies\n")
            release.write(directory / "native-parity.json", {"status": "passed", "cases": [{}, {}, {}]})
            record = {"source_revision": {"commit": self.sha, "clean": True}, "platform": platform,
                "wheel": wheel.name, "sha256": release.sha256(wheel),
                "installed_smoke": "passed with isolated Python", "qualification": "passed with isolated Python",
                "installed_tests": {"discovered": 10, "passed": 9, "skipped": 1},
                "numpy_validation": [{"version": v, "status": "passed", "tests": 9}
                                     for v in release.NUMPY_TEST_VERSIONS],
                "occt_source": {"file": "occt-8.0.1.tar.gz", "sha256": release.sha256(directory / "occt-8.0.1.tar.gz")},
                "repair": {"audit": "dependency-audit.txt", "audit_sha256": release.sha256(directory / "dependency-audit.txt")},
                "native_parity": {"file": "native-parity.json", "sha256": release.sha256(directory / "native-parity.json")}}
            release.write(directory / "wheel-validation.json", record)
            consumer = self.consumers / ("paacs-mesh-installed-" + platform)
            consumer.mkdir(parents=True)
            release.write(consumer / "consumer-validation.json", record)

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

    def test_three_qualified_wheels_and_matching_source_are_staged_separately(self):
        manifest = self.assemble()
        self.assertEqual(len(list((self.output / "dist").iterdir())), 3)
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
        path = self.consumers / "paacs-mesh-installed-linux-x86_64/consumer-validation.json"
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


if __name__ == "__main__":
    unittest.main()
