"""Pure stdlib/mock boundary-bundle admission tests: no PAACS/native imports."""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import struct
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch

SOURCE = Path(__file__).resolve().parents[2] / "python" / "paacs_mesh"
PACKAGE = "_paacs_experimental_boundary_unit"
package = types.ModuleType(PACKAGE)
package.__path__ = [str(SOURCE)]
sys.modules[PACKAGE] = package


class Model:
    kernel_version = "test-kernel"
    body_tags = [1]

    def __init__(self, source=b"owned source fixture", path=None):
        self._native = types.SimpleNamespace(_source_bytes=source)
        self.source_sha256 = hashlib.sha256(source).hexdigest()
        self.info = {"source_path": str(path), "file_bytes":len(source)}

    def _assert_identity(self):
        return None

    def inspect(self, kind, offset=0, limit=100, topology=False):
        values = {"bodies": [{"tag": 1, "bounds":{"lower":[0.,0.,0.],"upper":[10.,10.,10.]}, "face_uses": [{"face_tag": 2, "orientation": "forward", "sense": 1}]}],
                  "faces": [{"tag": 2, "surface": "plane", "owners": [{"body_tag": 1, "sense": 1, "outward_normal": [0, 0, 1]}], "wires": []}],
                  "edges": [], "vertices": []}[kind]
        return {"offset": offset, "total": len(values), "next_offset": None, "entities": values[offset:offset+limit]}


cad = types.ModuleType(PACKAGE + ".cad")
cad.CadModel, cad.read_step = Model, Mock()
sys.modules[cad.__name__] = cad
native = types.ModuleType(PACKAGE + "._native")
sys.modules[native.__name__] = native
spec = importlib.util.spec_from_file_location(PACKAGE + ".experimental", SOURCE / "experimental.py")
experimental = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = experimental
spec.loader.exec_module(experimental)


class Boundary:
    def __init__(self, options):
        self.options = options
        self.payload = b"mock archive; native authority deliberately mocked"
        self.report = {"schema": "paacs-experimental-boundary-report-v1", "prepared": True,
                       "size": options["size"], "effective_chord_error": options["max_chord_error"] or .01,
                       "maximum_relative_volume_error": options["max_relative_volume_error"],
                       "nodes": 3, "triangles": 1, "solid_tag": 1}
        self.to_numpy = Mock(return_value="owned optional arrays")

    def archive(self, maximum):
        if len(self.payload) > maximum:
            raise ValueError("byte allowance")
        return self.payload


class ExperimentalBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.mesh = types.SimpleNamespace(_experimental_set_python_timing=Mock())
        self.native = types.SimpleNamespace(engine_version="1.3.0-dev.0", source_revision="a"*40,
            __file__=str(Path(__file__).resolve()),
            source_dirty=False, _experimental_prepare_boundary=Mock(side_effect=lambda model, options: Boundary(options)),
            _experimental_rebind_boundary=Mock(side_effect=lambda model, archive, options: Boundary(options)),
            _experimental_mesh_boundary=Mock(return_value=self.mesh))
        self.patcher = patch.object(experimental, "_native", self.native)
        self.patcher.start()
        self.addCleanup(self.patcher.stop)

    def prepared(self):
        return experimental.prepare_boundary(Model(), size=.2)

    def read_model(self, path, **kwargs):
        self.assertTrue(kwargs["retain_source"])
        return Model(Path(path).read_bytes())

    def test_source_fingerprint_and_archive_identity_before_volume(self):
        boundary = self.prepared()
        self.native._experimental_mesh_boundary.assert_not_called()
        self.assertEqual(boundary.source_sha256, hashlib.sha256(boundary.source_bytes).hexdigest())
        self.assertEqual(boundary.boundary_sha256, hashlib.sha256(boundary._archive).hexdigest())
        report = boundary.report
        report["prepared"] = False
        self.assertTrue(boundary.report["prepared"])

    def test_dirty_native_builds_at_same_revision_have_distinct_binary_identity(self):
        with tempfile.TemporaryDirectory() as folder:
            one,two=Path(folder)/"one.pyd",Path(folder)/"two.pyd"
            one.write_bytes(b"first dirty binary");two.write_bytes(b"second dirty binary")
            self.native.source_dirty=True
            self.native.__file__=str(one);first=experimental._build()
            self.native.__file__=str(two);second=experimental._build()
            self.assertEqual(first["source_revision"],second["source_revision"])
            self.assertNotEqual(first["native_extension_sha256"],second["native_extension_sha256"])
            self.assertEqual(first["native_extension_sha256"],hashlib.sha256(one.read_bytes()).hexdigest())

    def test_none_and_empty_supplied_are_distinct_and_cloud_is_owned_ieee64(self):
        boundary = self.prepared()
        boundary.mesh()
        self.assertIsNone(self.native._experimental_mesh_boundary.call_args.args[-1])
        boundary.mesh([])
        self.assertEqual(self.native._experimental_mesh_boundary.call_args.args[-1], b"")
        cloud = [[.1, -.0, .3]]
        boundary.mesh(cloud)
        copied = self.native._experimental_mesh_boundary.call_args.args[-1]
        cloud[0][0] = 999
        self.assertEqual(copied, struct.pack("<3d", .1, -.0, .3))

    def test_invalid_clouds_refuse_before_native(self):
        boundary = self.prepared()
        for cloud in ("abc", b"", iter([(1,2,3)]), [(1,2)], [(True,2,3)], [(math.inf,2,3)], [(math.nan,2,3)]):
            with self.subTest(cloud=repr(cloud)), self.assertRaises((ValueError, TypeError)):
                boundary.mesh(cloud)
        self.native._experimental_mesh_boundary.assert_not_called()
        with self.assertRaises(ValueError):
            experimental._point_bytes([[1,2,3]]*2, {"max_nodes":1,"max_working_bytes":1024})

    def test_point_sequence_uses_bounded_indexing_not_untrusted_iteration(self):
        class Point:
            def __len__(self): return 3
            def __getitem__(self, index): return index + .5
            def __iter__(self): raise AssertionError("Iteration can ignore advertised length")
        self.assertEqual(experimental._point_bytes([Point()], {"max_nodes":1,"max_working_bytes":1024}), struct.pack("<3d", .5,1.5,2.5))

    def test_source_change_and_unsupported_controls_refuse_before_preparation(self):
        model = Model()
        model.source_sha256 = "0"*64
        with self.assertRaisesRegex(ValueError,"source"):
            experimental.prepare_boundary(model,size=.2)
        for kwargs in ({"face_sizes":[]},{"max_nodes":True},{"max_seconds":math.inf},{"max_seconds":86401}):
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                experimental.prepare_boundary(Model(),size=.2,**kwargs)
        self.native._experimental_prepare_boundary.assert_not_called()

    def test_nonretained_source_file_is_bounded_and_hash_checked(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"source.step";path.write_bytes(b"owned source fixture")
            model=Model(path=path);model._native._source_bytes=None
            self.assertEqual(experimental.prepare_boundary(model,size=.2).source_bytes,path.read_bytes())
            path.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError,"source"):
                experimental.prepare_boundary(model,size=.2)

    def test_save_load_reimports_and_rebinds_with_exact_bytes(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.prepared().save(Path(folder)/"bundle")
            before=(path/"manifest.json").read_bytes()
            with patch.object(experimental,"read_step",side_effect=self.read_model) as read:
                loaded=experimental.load_boundary(path)
            read.assert_called_once();self.native._experimental_rebind_boundary.assert_called_once()
            self.assertTrue(loaded.report["saved_preparation"]["prepared"])
            self.assertEqual(loaded.source_bytes,b"owned source fixture")
            with self.assertRaises(FileExistsError):self.prepared().save(path)
            self.assertEqual((path/"manifest.json").read_bytes(),before)

    def test_tightened_prepared_chord_rebinds_effective_value_preserving_request(self):
        def tightened(model,options):
            boundary=Boundary(options);boundary.report["effective_chord_error"]=.005;return boundary
        self.native._experimental_prepare_boundary.side_effect=tightened
        for requested in (None,.02):
            with self.subTest(requested=requested),tempfile.TemporaryDirectory() as folder:
                prepared=experimental.prepare_boundary(Model(),size=.2,max_chord_error=requested)
                path=prepared.save(Path(folder)/"boundary")
                with patch.object(experimental,"read_step",side_effect=self.read_model):loaded=experimental.load_boundary(path)
                self.assertEqual(self.native._experimental_rebind_boundary.call_args.args[2]["max_chord_error"],.005)
                loaded.mesh([])
                self.assertEqual(self.native._experimental_mesh_boundary.call_args.args[2]["max_chord_error"],requested)

    def test_replay_success_and_refusal_include_python_cost_snapshot(self):
        prepared=self.prepared()
        prepared.mesh([[.1,.2,.3]])
        total,packing,extension_hash=self.mesh._experimental_set_python_timing.call_args.args
        self.assertGreaterEqual(total,packing)
        self.assertEqual(extension_hash,hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        error=ValueError("Native refused")
        error.report={"reason":"geometry_error","experimental_boundary_volume":{"mode":"supplied"}}
        self.native._experimental_mesh_boundary.side_effect=error
        with self.assertRaises(ValueError) as failed:prepared.mesh([])
        timing=failed.exception.report["experimental_boundary_volume"]
        self.assertEqual(timing["mode"],"supplied")
        self.assertTrue(timing["native_binding_called"])
        self.assertEqual(timing["supplied_input_sha256"],hashlib.sha256(b"").hexdigest())
        self.assertGreaterEqual(timing["python_total_seconds"],timing["point_packing_seconds"])
        with self.assertRaises(ValueError) as invalid:prepared.mesh([(True,2,3)])
        self.assertFalse(invalid.exception.report["experimental_boundary_volume"]["native_binding_called"])

    def test_refusal_cloud_is_bounded_hash_checked_and_not_a_published_mesh(self):
        error=ValueError("Volume refused")
        error.initial_points_bytes=struct.pack("<3d",.1,.2,.3)
        error.report={"mesh_published":False,"experimental_boundary_volume":{
            "initial_points_available":True,"initial_points":1,
            "initial_points_sha256":hashlib.sha256(error.initial_points_bytes).hexdigest()}}
        self.assertEqual(experimental.initial_points(error),((.1,.2,.3),))
        self.assertFalse(error.report["mesh_published"])
        with self.assertRaisesRegex(ValueError,"max_bytes"):experimental.initial_points(error,max_bytes=65536)
        error.initial_points_bytes=struct.pack("<3d",.4,.5,.6)
        with self.assertRaisesRegex(ValueError,"scalar receipt"):experimental.initial_points(error)
        error.initial_points_bytes=None;error.initial_points_unavailable_reason="not_captured"
        with self.assertRaisesRegex(ValueError,"not_captured"):experimental.initial_points(error)

    def test_captured_empty_failure_cloud_differs_from_unavailable(self):
        error=ValueError("Boundary-only volume refused")
        error.initial_points_bytes=b""
        error.report={"experimental_boundary_volume":{"initial_points_available":True,"initial_points":0,
            "initial_points_sha256":hashlib.sha256(b"").hexdigest()}}
        self.assertEqual(experimental.initial_points(error),())

    def test_hash_tampering_rejected_before_import_or_rebind(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.prepared().save(Path(folder)/"bundle")
            (path/"source.step").write_bytes(b"changed")
            with patch.object(experimental,"read_step") as read,self.assertRaisesRegex(ValueError,"hash/length"):
                experimental.load_boundary(path)
            read.assert_not_called();self.native._experimental_rebind_boundary.assert_not_called()

    def test_manifest_paths_build_units_and_bool_counts_rejected_before_native(self):
        mutations=(lambda x:x["files"].update({"../outside":x["files"].pop("source.step")}),
                   lambda x:x["build"].update(source_revision="b"*40),lambda x:x.update(units="m"),
                   lambda x:x["files"]["source.step"].update(bytes=True))
        for mutate in mutations:
            with self.subTest(mutate=mutate),tempfile.TemporaryDirectory() as folder:
                path=self.prepared().save(Path(folder)/"bundle");manifest=json.loads((path/"manifest.json").read_text());mutate(manifest)
                (path/"manifest.json").write_text(json.dumps(manifest))
                with patch.object(experimental,"read_step") as read,self.assertRaises(ValueError):experimental.load_boundary(path)
                read.assert_not_called()

    def test_duplicate_and_nonfinite_json_refuse(self):
        for data in (b'{"schema":1,"schema":2}',b'{"x":NaN}',b'{"x":1e999}'):
            with self.assertRaises(ValueError):experimental._json(data)

    def test_matching_hashes_and_claimed_prepared_flag_never_bypass_native_rebind(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.prepared().save(Path(folder)/"bundle")
            self.native._experimental_rebind_boundary.side_effect=ValueError("Native CAD boundary rejected")
            with patch.object(experimental,"read_step",side_effect=self.read_model),self.assertRaisesRegex(ValueError,"Native CAD"):
                experimental.load_boundary(path)
            self.native._experimental_rebind_boundary.assert_called_once()

    def test_fresh_import_topology_and_kernel_must_match(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.prepared().save(Path(folder)/"bundle")
            def different(path,**kwargs):
                model=self.read_model(path,**kwargs);model.kernel_version="different kernel";return model
            with patch.object(experimental,"read_step",side_effect=different),self.assertRaisesRegex(ValueError,"identity"):
                experimental.load_boundary(path)
            self.native._experimental_rebind_boundary.assert_not_called()

    def test_optional_array_copy_budget_validation(self):
        boundary=self.prepared()
        self.assertEqual(boundary.to_numpy(max_bytes=12345),"owned optional arrays")
        boundary._boundary.to_numpy.assert_called_once_with(max_bytes=12345)
        with self.assertRaises(ValueError):boundary.to_numpy(max_bytes=True)


if __name__ == "__main__":unittest.main()
