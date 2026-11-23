"""The release plan must qualify every promised ABI, architecture and consumer OS."""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from contextlib import redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import python_compatibility as compatibility


class PythonCompatibilityTests(unittest.TestCase):
    def setUp(self):
        self.manifest = compatibility.load_manifest()

    def test_complete_cartesian_build_matrix_and_unique_artifacts(self):
        rows = compatibility.targets()
        self.assertEqual(len(rows), 20)
        self.assertEqual({(row["python"], row["platform"]) for row in rows},
                         {(version, platform) for version in ("3.12", "3.13", "3.14", "3.15")
                          for platform in ("windows-x86_64", "macos-arm64", "linux-x86_64",
                                           "linux-arm64", "windows-arm64")})
        for field in ("id", "artifact", "job_name"):
            self.assertEqual(len({row[field] for row in rows}), 20)
        for row in rows:
            self.assertEqual(compatibility.target(row["platform"], row["python"]), row)
            self.assertEqual(compatibility.target_by_id(row["id"]), row)
            self.assertEqual(row["python_tag"], "cp" + row["python"].replace(".", ""))

    def test_independent_consumers_cover_each_wheel_and_all_macos_versions(self):
        rows = compatibility.consumers()
        self.assertEqual(len(rows), 28)
        for field in ("id", "consumer_artifact", "job_name"):
            self.assertEqual(len({row[field] for row in rows}), 28)
        self.assertEqual({row["target_id"] for row in rows},
                         {row["id"] for row in compatibility.targets()})
        for version in compatibility.PYTHON_VERSIONS:
            mac = [row for row in rows if row["platform"] == "macos-arm64" and row["python"] == version]
            self.assertEqual({(row["runner"], row["required_os_major"]) for row in mac},
                             {("macos-15", 15), ("macos-26", 26), ("xcode-27", 27)})
            self.assertEqual(len({row["artifact"] for row in mac}), 1,
                             "Newer macOS versions must consume the identical deployment-target-15 wheel")
            for row in mac:
                self.assertEqual(compatibility.consumer(row["id"]), row)
                self.assertEqual(row["platform_tag"], "macosx_15_0_arm64")
                self.assertEqual(row["preview"], row["required_os_major"] == 27)

    def test_only_compatible_numpy_versions_are_mandatory(self):
        for row in compatibility.targets():
            expected = (("1.26.4", "2.5.5") if row["python"] == "3.12"
                        and row["platform"] != "windows-arm64" else ("2.5.5",))
            self.assertEqual(compatibility.numpy_versions(row["platform"], row["python"]), expected)
            self.assertEqual(row["numpy_versions"], list(expected))

    def test_sdk_builds_are_python_independent_and_architecture_specific(self):
        rows = compatibility.platforms()
        self.assertEqual(len(rows), 5)
        self.assertEqual(len({row["sdk_artifact"] for row in rows}), 5)
        linux = [row for row in rows if row["id"].startswith("linux-")]
        self.assertEqual({row["image_architecture"] for row in linux}, {"amd64", "arm64"})
        self.assertTrue(all("@sha256:" in row["image"] for row in linux))
        self.assertNotEqual(linux[0]["image"], linux[1]["image"])

    def test_windows_arm64_toolset_pin_reaches_sdk_wheels_and_consumers_only(self):
        for rows in (compatibility.platforms(), compatibility.targets(), compatibility.consumers()):
            for row in rows:
                if row["platform"] == "windows-arm64":
                    self.assertEqual(row["msvc_toolset"], "14.44")
                else:
                    self.assertNotIn("msvc_toolset", row,
                                     "The ARM64 compiler workaround must not change other toolchains")

    def test_windows_arm64_toolset_pin_cannot_be_missing_or_changed(self):
        for selected in (None, "14.52", "v143", "14.44-arm64ec", 14.44, True):
            with self.subTest(selected=selected):
                value = deepcopy(self.manifest)
                row = next(row for row in value["platforms"] if row["id"] == "windows-arm64")
                if selected is None:
                    row.pop("msvc_toolset")
                else:
                    row["msvc_toolset"] = selected
                with self.assertRaises(RuntimeError):
                    compatibility.validate_manifest(value)

    def test_toolset_pin_is_rejected_on_other_platforms(self):
        for platform in ("windows-x86_64", "macos-arm64", "linux-x86_64", "linux-arm64"):
            with self.subTest(platform=platform):
                value = deepcopy(self.manifest)
                next(row for row in value["platforms"] if row["id"] == platform)["msvc_toolset"] = "14.44"
                with self.assertRaises(RuntimeError):
                    compatibility.validate_manifest(value)

    def test_manifest_rejects_missing_extra_and_duplicate_versions(self):
        for versions in (["3.12", "3.13", "3.14"], ["3.12", "3.13", "3.14", "3.15", "3.16"],
                         ["3.12", "3.13", "3.14", "3.14"]):
            with self.subTest(versions=versions):
                value = deepcopy(self.manifest)
                value["python_versions"] = versions
                with self.assertRaises(RuntimeError):
                    compatibility.validate_manifest(value)

    def test_manifest_rejects_missing_extra_and_duplicate_platforms(self):
        for mode in ("missing", "extra", "duplicate"):
            with self.subTest(mode=mode):
                value = deepcopy(self.manifest)
                if mode == "missing":
                    value["platforms"].pop()
                elif mode == "extra":
                    value["platforms"].append(deepcopy(value["platforms"][0]))
                else:
                    value["platforms"][-1] = deepcopy(value["platforms"][0])
                with self.assertRaises(RuntimeError):
                    compatibility.validate_manifest(value)

    def test_manifest_rejects_omitted_macos_consumer_and_unreviewed_labels(self):
        for mutation in (lambda row: row["consumers"].pop(),
                         lambda row: row["consumers"][2].update(runner="macos-27"),
                         lambda row: row["consumers"][2].update(required_os_major=26),
                         lambda row: row["consumers"][2].update(preview=False),
                         lambda row: row.update(minimum_macos="26.0")):
            value = deepcopy(self.manifest)
            mutation(next(row for row in value["platforms"] if row["id"] == "macos-arm64"))
            with self.assertRaises(RuntimeError):
                compatibility.validate_manifest(value)

    def test_manifest_rejects_wrong_or_floating_linux_architecture(self):
        for change in ({"image": "quay.io/pypa/manylinux_2_28_aarch64:latest"},
                       {"image_architecture": "amd64"}, {"minimum_glibc": "2.39"},
                       {"image_config_sha256": "not-a-digest"}):
            with self.subTest(change=change):
                value = deepcopy(self.manifest)
                next(row for row in value["platforms"] if row["id"] == "linux-arm64").update(change)
                with self.assertRaises(RuntimeError):
                    compatibility.validate_manifest(value)

    def test_unknown_settings_and_invalid_numpy_releases_fail_closed(self):
        for change in ({"allow_failure": True}, {"requires_python": ">=3.12"},
                       {"numpy_current": "latest"}, {"numpy_current": "1.26.4"}):
            value = deepcopy(self.manifest)
            value.update(change)
            with self.assertRaises(RuntimeError):
                compatibility.validate_manifest(value)

    def test_lookups_reject_unqualified_abis_and_platforms(self):
        for version in ("3.11", "3.16", "3.15t", "cp315", (3, 15), True):
            with self.subTest(version=version), self.assertRaises(RuntimeError):
                compatibility.python_tag(version)
        for identity in ("cp315t-windows-arm64", "cp312-macos-x86_64", "", None):
            with self.subTest(identity=identity), self.assertRaises(RuntimeError):
                compatibility.target_by_id(identity)
        with self.assertRaises(RuntimeError):
            compatibility.target("macos-x86_64", "3.12")
        with self.assertRaises(RuntimeError):
            compatibility.consumer("cp315-macos-arm64-macos28")
        with patch.object(compatibility.sys, "version_info") as current:
            current.major, current.minor = 3, 14
            self.assertEqual(compatibility.python_tag(), "cp314")

    def test_returned_rows_cannot_modify_subsequent_release_plan(self):
        row = compatibility.target("windows-x86_64", "3.12")
        row["numpy_versions"].clear()
        row["artifact"] = "untrusted"
        self.assertEqual(compatibility.target("windows-x86_64", "3.12")["numpy_versions"],
                         ["1.26.4", "2.5.5"])
        value = compatibility.platforms()
        value[1]["consumers"].clear()
        self.assertEqual(len(compatibility.consumers()), 28)

    def test_cli_outputs_named_complete_matrices_for_actions(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "github-output"
            for kind, count in (("targets", 20), ("consumers", 28), ("platforms", 5)):
                console = io.StringIO()
                with redirect_stdout(console):
                    compatibility.main(["--" + kind, "--github-output", str(output)])
                self.assertEqual(len(json.loads(console.getvalue())["include"]), count)
            lines = output.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(lines), 3)
            self.assertEqual([line.split("=", 1)[0] for line in lines], ["targets", "consumers", "platforms"])
            self.assertTrue(all(len(json.loads(line.split("=", 1)[1])["include"]) > 0 for line in lines))

    def test_manifest_size_and_encoding_are_bounded(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "manifest.json"
            path.write_bytes(b" " * 64_001)
            with self.assertRaises(RuntimeError):
                compatibility.load_manifest(path)
            path.write_bytes(b"\xff")
            with self.assertRaises(RuntimeError):
                compatibility.load_manifest(path)

    def test_manifest_rejects_duplicate_json_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "manifest.json"
            encoded = json.dumps(self.manifest)
            encoded = encoded.replace('"architecture": "x64"',
                                      '"architecture": "arm64", "architecture": "x64"', 1)
            path.write_text(encoded, encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "Duplicate compatibility manifest key"):
                compatibility.load_manifest(path)


if __name__ == "__main__":
    unittest.main()
