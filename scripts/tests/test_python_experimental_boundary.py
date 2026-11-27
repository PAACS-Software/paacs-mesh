"""Pure stdlib/mock boundary-bundle admission tests: no PAACS/native imports."""
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0
import hashlib
import importlib.util
import json
import math
import os
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
            _experimental_mesh_boundary=Mock(return_value=self.mesh),
            _experimental_generate_initial_points=Mock(side_effect=self.generated))
        self.patcher = patch.object(experimental, "_native", self.native)
        self.patcher.start()
        self.addCleanup(self.patcher.stop)

    def prepared(self):
        return experimental.prepare_boundary(Model(), size=.2)

    def read_model(self, path, **kwargs):
        self.assertTrue(kwargs["retain_source"])
        return Model(Path(path).read_bytes())

    @staticmethod
    def generated(model, boundary, options):
        payload = struct.pack("<3d", .1, -.0, .3)
        return {"points": payload, "report": {
            "schema": "paacs-experimental-initial-points-report-v1", "generated": True,
            "reason": "none", "initial_points_available": True, "initial_points": 1,
            "initial_points_sha256": hashlib.sha256(payload).hexdigest(), "boundary_unchanged": True,
            "recovery_performed": False, "volume_construction_performed": False}}

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
        self.assertIsNone(self.native._experimental_mesh_boundary.call_args.args[3])
        boundary.mesh([])
        self.assertEqual(self.native._experimental_mesh_boundary.call_args.args[3], b"")
        cloud = [[.1, -.0, .3]]
        boundary.mesh(cloud)
        copied = self.native._experimental_mesh_boundary.call_args.args[3]
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

    def test_public_pipeline_facade_retains_experimental_class_identities(self):
        spec = importlib.util.spec_from_file_location(PACKAGE+".pipeline", SOURCE/"pipeline.py")
        facade = importlib.util.module_from_spec(spec);spec.loader.exec_module(facade)
        for name in facade.__all__:
            self.assertIs(getattr(facade, name), getattr(experimental, name))

    def test_capability_pipeline_discovery_is_additive_and_stdlib_only(self):
        spec = importlib.util.spec_from_file_location(PACKAGE+".diagnostics", SOURCE/"diagnostics.py")
        diagnostics = importlib.util.module_from_spec(spec);spec.loader.exec_module(diagnostics)
        with patch.object(package, "__version__", "1.2.1", create=True), \
                patch.object(package, "native_version", "1.2.1", create=True), \
                patch.object(cad, "_DEFAULT_RESOURCE_POLICY", "normal", create=True), \
                patch.object(cad, "_RESOURCE_POLICIES", ["normal", "bounded"], create=True):
            result = diagnostics.capabilities()
        self.assertEqual(result["schema"], "paacs-mesh-capabilities-v1")
        self.assertEqual(result["pipeline"]["module"], "paacs_mesh.pipeline")
        self.assertFalse(result["pipeline"]["numpy_required"])
        self.assertFalse(result["pipeline"]["initial_points"]["volume_certified"])
        self.assertFalse(result["pipeline"]["changes_default_meshing"])
        self.assertFalse(result["pipeline"]["includes_AI_or_training"])
        json.dumps(result, allow_nan=False)

    def test_generation_stops_before_any_volume_and_passes_owned_bytes_without_repacking(self):
        boundary = self.prepared();cloud = boundary.generate_initial_points()
        self.native._experimental_mesh_boundary.assert_not_called()
        self.assertEqual(cloud.points(), ((.1, -.0, .3),))
        self.assertEqual(experimental.initial_points(cloud), cloud.points())
        self.assertEqual(cloud.point_count, 1);self.assertEqual(cloud.units, "mm")
        self.assertFalse(cloud.report["recovery_performed"])
        self.assertFalse(cloud.report["volume_construction_performed"])
        self.assertNotIn("audit", cloud.report);self.assertFalse(cloud.report.get("mesh_published", False))
        with patch.object(experimental, "_point_bytes", side_effect=AssertionError("Do not repack an owned cloud")):
            boundary.mesh(cloud)
        self.assertIs(self.native._experimental_mesh_boundary.call_args.args[3], cloud._payload)
        with self.assertRaises(AttributeError):cloud._payload = b""
        with self.assertRaises(AttributeError):del cloud._payload
        copy = cloud.report;copy["generated"] = False
        self.assertTrue(cloud.report["generated"])
        options = cloud.options;options["size"] = 100
        self.assertEqual(cloud.options["size"], .2)

    def test_cloud_save_reload_is_hash_bound_and_does_not_import_cad_or_call_volume(self):
        boundary = self.prepared();cloud = boundary.generate_initial_points()
        with tempfile.TemporaryDirectory() as temporary:
            path = cloud.save(Path(temporary)/"points")
            with patch.object(experimental, "read_step") as read:
                loaded = experimental.load_initial_points(path)
            read.assert_not_called()
            self.assertEqual(loaded.points(), cloud.points())
            self.assertEqual(loaded.sha256, cloud.sha256)
            self.assertEqual(loaded.provenance, cloud.provenance)
            self.assertEqual(set(p.name for p in path.iterdir()), {"points.bin", "point-report.json", "manifest.json"})
            with self.assertRaises(FileExistsError):cloud.save(path)
        self.native._experimental_mesh_boundary.assert_not_called()

    def test_generated_empty_cloud_is_available_and_roundtrips_as_empty_supplied(self):
        def empty(*args):
            result = self.generated(*args);result["points"] = b""
            result["report"].update(initial_points=0, initial_points_sha256=hashlib.sha256(b"").hexdigest())
            return result
        self.native._experimental_generate_initial_points.side_effect = empty
        boundary = self.prepared();cloud = boundary.generate_initial_points()
        self.assertEqual(cloud.points(), ());self.assertEqual(cloud.point_count, 0)
        self.assertTrue(cloud.report["initial_points_available"])
        with tempfile.TemporaryDirectory() as temporary:
            loaded = experimental.load_initial_points(cloud.save(Path(temporary)/"empty"))
        boundary.mesh(loaded)
        self.assertEqual(self.native._experimental_mesh_boundary.call_args.args[3], b"")

    def test_cloud_source_boundary_build_and_options_identity_refuse_before_native(self):
        boundary = self.prepared();cloud = boundary.generate_initial_points()
        others = [experimental.prepare_boundary(Model(source=b"different CAD"), size=.2),
                  experimental.prepare_boundary(Model(), size=.3)]
        for other in others:
            with self.assertRaisesRegex(ValueError, "source/boundary/build/options"):
                other.mesh(cloud)
        self.native.source_revision = "b"*40
        with self.assertRaisesRegex(ValueError, "source/boundary/build/options"):
            boundary.mesh(cloud)
        self.native._experimental_mesh_boundary.assert_not_called()

    def test_cloud_build_identity_does_not_use_bool_integer_equality(self):
        boundary = self.prepared();cloud = boundary.generate_initial_points()
        metadata = cloud.provenance;metadata["build"]["source_dirty"] = 0
        report = cloud.report;report["provenance"] = metadata
        forged = experimental.InitialPointCloud(cloud._payload, cloud.options, metadata, report)
        with self.assertRaisesRegex(ValueError, "source/boundary/build/options"):
            boundary.mesh(forged)
        self.native._experimental_mesh_boundary.assert_not_called()

    def test_cloud_manifest_report_and_payload_tampering_refuse_without_native(self):
        mutations = (lambda m: m.update(point_count=True), lambda m: m.update(units="m"),
                     lambda m: m["build"].update(source_revision="b"*40),
                     lambda m: m["files"].update({"../points.bin": m["files"].pop("points.bin")}),
                     lambda m: m["files"]["points.bin"].update(bytes=10**12))
        cloud = self.prepared().generate_initial_points()
        for mutate in mutations:
            with tempfile.TemporaryDirectory() as temporary:
                path = cloud.save(Path(temporary)/"cloud");manifest = json.loads((path/"manifest.json").read_text())
                mutate(manifest);(path/"manifest.json").write_text(json.dumps(manifest))
                with self.assertRaises(ValueError):experimental.load_initial_points(path)
        with tempfile.TemporaryDirectory() as temporary:
            path = cloud.save(Path(temporary)/"cloud")
            (path/"points.bin").write_bytes(struct.pack("<3d", .4, .5, .6))
            with self.assertRaisesRegex(ValueError, "hash/length"):experimental.load_initial_points(path)
        self.native._experimental_mesh_boundary.assert_not_called()

    def test_generation_refusal_or_malformed_native_receipt_never_yields_partial_cloud(self):
        boundary = self.prepared()
        for field, value in (("generated", False), ("initial_points_available", False),
                             ("initial_points", True), ("recovery_performed", True),
                             ("volume_construction_performed", True), ("mesh_published", True),
                             ("units", "m"),
                             ("initial_points_sha256", "0"*64)):
            result = self.generated(None, None, None);result["report"][field] = value
            self.native._experimental_generate_initial_points.side_effect = None
            self.native._experimental_generate_initial_points.return_value = result
            with self.assertRaises(ValueError):boundary.generate_initial_points()
        failure = ValueError("Stage refused");failure.report = {"generated": False, "reason": "work_budget"}
        self.native._experimental_generate_initial_points.side_effect = failure
        with self.assertRaises(ValueError) as caught:boundary.generate_initial_points()
        self.assertFalse(caught.exception.report["generated"])
        self.assertTrue(caught.exception.report["native_binding_called"])
        self.assertGreaterEqual(caught.exception.report["python_generation_seconds"], 0.)
        self.native._experimental_mesh_boundary.assert_not_called()

    def test_cloud_copy_and_encoded_read_budgets_are_explicit(self):
        cloud = self.prepared().generate_initial_points()
        for value in (True, 0, 65536):
            with self.assertRaises(ValueError):cloud.points(max_bytes=value)
        with self.assertRaises(ValueError):cloud.to_numpy(max_bytes=65536)
        with tempfile.TemporaryDirectory() as temporary:
            path = cloud.save(Path(temporary)/"cloud")
            with self.assertRaises(ValueError):experimental.load_initial_points(path, max_bytes=23)
        payload = struct.pack("<3d", math.nan, 0., 0.)
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            experimental._point_payload(payload, cloud.options)

    def test_cad_incidence_records_are_bounded_copied_and_frozen_hash_checked(self):
        boundary = self.prepared();records = boundary.cad_topology()
        self.assertEqual(records["faces"][0]["owners"], [{"body_tag": 1, "sense": 1}])
        self.assertEqual(records["bodies"][0]["face_uses"][0]["face_tag"], 2)
        records["faces"][0]["owners"][0]["sense"] = -1
        self.assertEqual(boundary.cad_topology()["faces"][0]["owners"][0]["sense"], 1)
        with self.assertRaisesRegex(ValueError, "max_bytes"):boundary.cad_topology(max_bytes=1)

    def test_cloud_save_failure_rolls_back_before_manifest_publication(self):
        cloud = self.prepared().generate_initial_points();original = Path.open
        def fail_report(path, *args, **kwargs):
            if path.name == "point-report.json":
                raise OSError("Injected report write failure")
            return original(path, *args, **kwargs)
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)/"cloud"
            with patch.object(Path, "open", new=fail_report), self.assertRaises(OSError):
                cloud.save(destination)
            self.assertFalse(destination.exists())
            self.assertEqual(list(Path(temporary).iterdir()), [])

    def test_cloud_directory_publish_failure_cleans_staging_and_preserves_existing_destination(self):
        cloud = self.prepared().generate_initial_points()
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)/"cloud"
            with patch.object(experimental.os, "rename", side_effect=OSError("Injected rename failure")), self.assertRaises(OSError):
                cloud.save(destination)
            self.assertFalse(destination.exists())
            self.assertEqual(list(Path(temporary).iterdir()), [])
            destination.mkdir();marker = destination/"keep.txt";marker.write_bytes(b"existing user data")
            with self.assertRaises(FileExistsError):cloud.save(destination)
            self.assertEqual(marker.read_bytes(), b"existing user data")
            self.assertEqual(list(destination.iterdir()), [marker])
            self.assertEqual(list(Path(temporary).iterdir()), [destination])


    def test_loaded_boundary_protects_all_bundle_inputs_and_hardlinks_on_exports(self):
        spec = importlib.util.spec_from_file_location(PACKAGE+".pipeline_export", SOURCE/"pipeline_export.py")
        exported = importlib.util.module_from_spec(spec);spec.loader.exec_module(exported)
        with patch.dict(sys.modules, {PACKAGE+".pipeline_export": exported}), tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary);bundle = self.prepared().save(folder/"boundary")
            before = {name: (bundle/name).read_bytes() for name in (*experimental._FILES, "manifest.json")}
            with patch.object(experimental, "read_step", side_effect=self.read_model):
                loaded = experimental.load_boundary(bundle)
            expected = {bundle/name for name in (*experimental._FILES, "manifest.json")}
            self.assertTrue(expected <= set(loaded._protected_paths))
            calls = []
            def write(path, maximum):
                calls.append(Path(path));Path(path).write_bytes(b"mock stage export")
            loaded._boundary.write = write
            edge = types.SimpleNamespace(report={"complete": True}, write=write)
            with patch.object(self.native, "_experimental_cad_edges", return_value=edge, create=True), \
                    patch.object(exported.GeometryArtifact, "_snapshot", return_value=b"mock stage archive"):
                for name, artifact in (("surface", loaded.surface()), ("edges", loaded.edges())):
                    with self.assertRaises(ValueError):artifact.write(bundle/"boundary.paacsmesh")
                    alias = folder/(name+"-input-alias.paacsmesh")
                    alias.write_bytes(before["boundary.paacsmesh"])
                    original_samefile = os.path.samefile
                    def samefile(left, right):
                        if {Path(left).absolute(), Path(right).absolute()} == {alias, bundle/"boundary.paacsmesh"}:
                            return True
                        return original_samefile(left, right)
                    # Pure guard test only; real filesystem hardlinks are
                    # separately exercised in the installed Native suite.
                    with patch.object(os.path, "samefile", side_effect=samefile):
                        with self.assertRaises(ValueError):artifact.write(alias)
                    outside = folder/(name+"-outside.paacsmesh")
                    artifact.write(outside)
                    self.assertEqual(outside.read_bytes(), b"mock stage export")
            self.assertEqual(len(calls), 2)
            self.assertEqual(before, {name: (bundle/name).read_bytes() for name in before})
            loaded.mesh([])
            self.assertTrue(expected <= set(self.native._experimental_mesh_boundary.call_args.args[4]))

    def test_loaded_cloud_guards_are_joined_without_changing_bytes_or_metadata(self):
        boundary = self.prepared();cloud = boundary.generate_initial_points()
        self.assertEqual(cloud._protected_paths, ())
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary);boundary_bundle = boundary.save(folder/"boundary")
            point_bundle = cloud.save(folder/"points")
            with patch.object(experimental, "read_step", side_effect=self.read_model):
                loaded_boundary = experimental.load_boundary(boundary_bundle)
            loaded_cloud = experimental.load_initial_points(point_bundle)
            loaded_boundary.mesh(loaded_cloud)
            arguments = self.native._experimental_mesh_boundary.call_args.args
            self.assertEqual(arguments[3], cloud._payload)
            expected = {boundary_bundle/name for name in (*experimental._FILES, "manifest.json")}
            expected |= {point_bundle/name for name in (*experimental._POINT_FILES, "manifest.json")}
            self.assertTrue(expected <= set(arguments[4]))
            self.assertEqual(len(arguments[4]), len(set(arguments[4])))
            self.assertLessEqual(len(arguments[4]), 16)
            self.assertEqual(loaded_cloud.sha256, cloud.sha256)
            self.assertEqual(loaded_cloud.provenance, cloud.provenance)
            again = loaded_cloud.save(folder/"points-copy")
            for name in (*experimental._POINT_FILES, "manifest.json"):
                self.assertEqual((again/name).read_bytes(), (point_bundle/name).read_bytes())


if __name__ == "__main__":unittest.main()
