"""Reusable CAD SDKs retain exact source, toolchain and installed content identity."""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import copy
import hashlib
import io
import json
from pathlib import Path
import sys
import struct
import tarfile
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_python_sdk as sdk


class SDKReuseTests(unittest.TestCase):
    @staticmethod
    def binary(architecture=0x8664):
        payload = bytearray(128)
        payload[:2] = b"MZ"
        struct.pack_into("<I", payload, 0x3c, 64)
        payload[64:68] = b"PE\0\0"
        struct.pack_into("<H", payload, 68, architecture)
        return bytes(payload)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        prefix = self.root / "occt"
        for name in ("lib/cmake/opencascade/OpenCASCADEConfig.cmake", "share/doc/LICENSE_LGPL_21.txt",
                     "share/doc/OCCT_LGPL_EXCEPTION.txt", "lib/TKernel.dll", "lib/TKDESTEP.dll"):
            target = prefix / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(self.binary() if name.endswith(".dll") else name.encode())
        archive = self.root / "downloads/occt-8.0.1.tar.gz"
        archive.parent.mkdir()
        archive.write_bytes(b"pinned source fixture")
        self.context = {"platform": "windows-x86_64", "occt": {"VERSION": "8.0.1",
            "SHA512": hashlib.sha512(archive.read_bytes()).hexdigest()},
            "compiler": {"id": "MSVC", "version": "19.50"}, "recipes": {"step": "a" * 64},
            "macos_deployment_target": None, "manylinux_image": None}
        self.receipt = {"schema": sdk.SCHEMA, "context": self.context, "files": sdk.inventory(prefix)}
        (self.root / "sdk-receipt.json").write_text(json.dumps(self.receipt), encoding="utf-8")

    def test_clean_sdk_can_be_reused_only_with_the_exact_toolchain_context(self):
        self.assertEqual(sdk.verify_sdk(self.root, self.context), self.receipt)
        changed = copy.deepcopy(self.context)
        changed["compiler"]["version"] = "19.51"
        with self.assertRaisesRegex(RuntimeError, "toolchain"):
            sdk.verify_sdk(self.root, changed)

    def test_modified_added_or_deleted_installed_files_are_refused(self):
        library = self.root / "occt/lib/TKernel.dll"
        library.write_bytes(b"changed")
        with self.assertRaisesRegex(RuntimeError, "hashes"):
            sdk.verify_sdk(self.root)
        library.write_bytes(self.binary())
        extra = self.root / "occt/lib/unreviewed.dll"
        extra.write_bytes(b"unexpected")
        with self.assertRaisesRegex(RuntimeError, "hashes"):
            sdk.verify_sdk(self.root)
        extra.unlink()
        library.unlink()
        with self.assertRaisesRegex(RuntimeError, "hashes"):
            sdk.verify_sdk(self.root)

    def test_source_archive_tampering_cannot_be_hidden_by_installed_hashes(self):
        (self.root / "downloads/occt-8.0.1.tar.gz").write_bytes(b"another source")
        with self.assertRaisesRegex(RuntimeError, "pinned OCCT"):
            sdk.verify_sdk(self.root)

    def test_cache_fingerprint_changes_for_each_binary_compatibility_input(self):
        original = sdk.cache_key(self.context)
        for key, value in (("platform", "windows-arm64"), ("compiler", {"version": "20"}),
                           ("recipes", {"step": "b" * 64}), ("macos_deployment_target", "26.0"),
                           ("manylinux_image", "another-digest"), ("occt", {"VERSION": "9"})):
            changed = copy.deepcopy(self.context)
            changed[key] = value
            self.assertNotEqual(sdk.cache_key(changed), original)
        self.assertEqual(sdk.cache_key(json.loads(json.dumps(self.context))), original)

    def test_wrong_architecture_is_rejected_even_with_matching_inventory(self):
        (self.root / "occt/lib/TKDESTEP.dll").write_bytes(self.binary(0xaa64))
        self.receipt["files"] = sdk.inventory(self.root / "occt")
        (self.root / "sdk-receipt.json").write_text(json.dumps(self.receipt), encoding="utf-8")
        with self.assertRaisesRegex(RuntimeError, "arm64.*x86_64"):
            sdk.verify_sdk(self.root, self.context)

    def test_verified_installed_prefix_survives_archive_relocation(self):
        archive_path = self.root / "sdk.tar.gz"
        with tarfile.open(archive_path, "w:gz", dereference=False) as archive:
            for item in (self.root / "occt", self.root / "downloads", self.root / "sdk-receipt.json"):
                archive.add(item, arcname=item.name)
        relocated = self.root / "relocated"
        with tarfile.open(archive_path, "r:gz") as archive:
            archive.extractall(relocated, filter="data")
        self.assertEqual(sdk.verify_sdk(relocated, self.context), self.receipt)
        (relocated / "occt/lib/TKernel.dll").write_bytes(b"tampered")
        with self.assertRaises(RuntimeError):
            sdk.verify_sdk(relocated, self.context)

    def test_unpack_refuses_archive_traversal_and_absolute_links(self):
        for name, link in (("../escaped.bin", None), ("unsafe-link", "/outside-sdk")):
            archive_path = self.root / "unsafe.tar.gz"
            with tarfile.open(archive_path, "w:gz") as archive:
                item = tarfile.TarInfo(name)
                if link:
                    item.type, item.linkname = tarfile.SYMTYPE, link
                    archive.addfile(item)
                else:
                    item.size = 4
                    archive.addfile(item, io.BytesIO(b"data"))
            destination = self.root / ("link-extraction" if link else "path-extraction")
            with patch.object(sys, "argv", ["sdk", "unpack", "--sdk", str(destination),
                                            "--archive", str(archive_path)]):
                with self.assertRaises(tarfile.FilterError):
                    sdk.main()
            self.assertFalse((self.root / "escaped.bin").exists())


if __name__ == "__main__":
    unittest.main()
