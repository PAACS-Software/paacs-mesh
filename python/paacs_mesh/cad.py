"""Bounded CAD inspection, explicit face selections and companion intent."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from dataclasses import dataclass
import hashlib
import json
import math
import numbers
from pathlib import Path

from . import _native
from ._cad_common import (
    MAX_JSON_BYTES, IntentError, SelectionError, boolean, bounds, close, integer, json_copy,
    json_text, metric, number, point, sequence, text,
)


_KINDS = {"faces": 2, "bodies": 3, "edges": 1, "vertices": 0}
_QUERY_FIELDS = {"tags", "name", "group", "surface", "body", "area", "radius",
                 "centroid_box", "normal", "adjacent_to", "expect", "tolerance"}
_MAX_FACE_CONTROLS = 10000
_DEFAULT_CONSTRUCTION_EFFORT = "legacy"
_DEFAULT_RESOURCE_POLICY = "normal"
_RESOURCE_POLICIES = ("normal", "bounded")
_RESOURCE_FIELDS = ("resource_policy",)
_FIDELITY_FIELDS = ("max_chord_error", "max_relative_volume_error")
_CURVATURE_FIELDS = ("curvature_points_per_turn",)
_CONSTRUCTION_LIMITS = ("max_surface_points", "max_surface_chord_points",
                        "max_interior_candidates", "max_interior_sites")
_CONSTRUCTION_FIELDS = ("construction_effort", *_CONSTRUCTION_LIMITS)


def _resource_policy_value(value):
    if not isinstance(value, str) or value not in _RESOURCE_POLICIES:
        raise ValueError("resource_policy must be normal or bounded")
    return value


def _fidelity_options(max_chord_error=None, max_relative_volume_error=None):
    result = {}
    for name, value in zip(_FIDELITY_FIELDS, (max_chord_error, max_relative_volume_error)):
        if value is not None:
            value = number(value, name, minimum=0)
            if value == 0 or (name == "max_relative_volume_error" and value >= 1):
                raise ValueError(f"{name} must be positive" + (" and less than one" if name == "max_relative_volume_error" else ""))
        result[name] = value
    return result


def _curvature_value(value):
    if value is not None:
        value = integer(value, "curvature_points_per_turn", minimum=0)
        if value > 64:
            raise ValueError("curvature_points_per_turn must be between zero and 64")
    return value


def _construction_options(construction_effort=_DEFAULT_CONSTRUCTION_EFFORT, *, max_surface_points=None,
                          max_surface_chord_points=None, max_interior_candidates=None,
                          max_interior_sites=None):
    if not isinstance(construction_effort, str) or construction_effort not in ("legacy", "auto", "thorough"):
        raise ValueError("construction_effort must be legacy, auto or thorough")
    result = {"construction_effort": construction_effort}
    for name, value in zip(_CONSTRUCTION_LIMITS, (max_surface_points, max_surface_chord_points,
                                                max_interior_candidates, max_interior_sites)):
        result[name] = None if value is None else integer(value, name)
    return result


def _tags(value, name):
    values = [value] if isinstance(value, numbers.Integral) else sequence(value, name, maximum=1000000)
    result = [integer(item, name) for item in values]
    if len(result) != len(set(result)):
        raise ValueError(f"{name} contains duplicate tags")
    return tuple(sorted(result))


def _expectation(value):
    if isinstance(value, numbers.Integral):
        count = integer(value, "expect")
        return [count, count]
    lo, hi = sequence(value, "expect", length=2)
    lo, hi = integer(lo, "expect minimum", minimum=0), integer(hi, "expect maximum", minimum=0)
    if lo > hi:
        raise ValueError("expect minimum must not exceed maximum")
    return [lo, hi]


def _query_values(values):
    if not isinstance(values, dict) or set(values) - _QUERY_FIELDS:
        raise ValueError("Selection query contains unsupported fields")
    query = {}
    for key in ("tags", "adjacent_to"):
        if values.get(key) is not None:
            query[key] = list(_tags(values[key], key))
    for key in ("name", "group", "surface"):
        if values.get(key) is not None:
            query[key] = text(values[key], key)
    if values.get("body") is not None:
        query["body"] = integer(values["body"], "body")
    for key in ("area", "radius"):
        if values.get(key) is not None:
            query[key] = metric(values[key], key)
    if values.get("centroid_box") is not None:
        query["centroid_box"] = bounds(values["centroid_box"], "centroid_box")
    if values.get("normal") is not None:
        direction = point(values["normal"], "normal")
        scale = max(abs(item) for item in direction)
        if scale == 0:
            raise ValueError("normal must have nonzero length")
        scaled = [item / scale for item in direction]
        length = math.hypot(*scaled)
        query["normal"] = [item / length for item in scaled]
    if values.get("expect") is not None:
        query["expect"] = _expectation(values["expect"])
    query["tolerance"] = number(values.get("tolerance", 1e-6), "tolerance", minimum=0)
    return query


def _required(record, key, *, nested=None):
    container = record if nested is None else record.get(nested, {})
    value = container.get(key) if isinstance(container, dict) else None
    if value is None:
        raise SelectionError("metric_unavailable", f"Face {record['tag']} has no usable {key}",
                             {"face_tag": record["tag"], "metric": key})
    return value


def _metric_matches(actual, requested, tolerance):
    actual = number(actual, "native metric", minimum=0)
    if isinstance(requested, list):
        lo, hi = requested
        return (actual >= lo - tolerance * max(1, abs(lo)) and
                actual <= hi + tolerance * max(1, abs(hi)))
    return math.isclose(actual, requested, rel_tol=tolerance, abs_tol=tolerance)


def _face_guard(record):
    guard = {"tag": integer(record["tag"], "face tag"), "surface": text(record["surface"], "surface")}
    owners = record.get("owners")
    if not isinstance(owners, list):
        raise ValueError("Native face ownership metadata is unavailable")
    guard["owners"] = sorted({integer(owner["body_tag"], "body tag") for owner in owners})
    guard["owner_uses"] = []
    for owner in owners:
        sense = integer(owner["sense"], "owner sense", minimum=-1)
        if sense not in (-1, 0, 1):
            raise ValueError("Native owner sense must be -1, 0 or 1")
        outward = owner.get("outward_normal")
        guard["owner_uses"].append({"body_tag": integer(owner["body_tag"], "body tag"), "sense": sense,
                                    "outward_normal": None if outward is None else point(outward, "outward_normal")})
    guard["owner_uses"].sort(key=lambda owner: (owner["body_tag"], owner["sense"], json_text(owner["outward_normal"])))
    if record.get("normal") is not None:
        guard["normal"] = point(record["normal"], "normal")
    for key in ("area", "tolerance"):
        if record.get(key) is not None:
            guard[key] = number(record[key], key, minimum=0)
    if record.get("area_centroid") is not None:
        guard["area_centroid"] = point(record["area_centroid"], "area_centroid")
    if record.get("bounds") is not None:
        box = record["bounds"]
        lo, hi = bounds([box["lower"], box["upper"]], "bounds")
        guard["bounds"] = {"lower": lo, "upper": hi}
    radius = record.get("geometry", {}).get("radius")
    if radius is not None:
        guard["radius"] = number(radius, "radius", minimum=0)
    return guard


@dataclass(frozen=True)
class FaceSelection:
    """Read-only face IDs owned by one imported model, with copied guards."""

    __slots__ = ("_owner", "_identity", "_tags", "_body", "_query_json", "_guards_json")

    _owner: "CadModel"
    _identity: int
    _tags: tuple
    _body: int | None
    _query_json: str
    _guards_json: str

    @property
    def tags(self):
        return self._tags

    @property
    def body(self):
        return self._body

    @property
    def model_identity(self):
        return self._identity

    @property
    def query(self):
        return json.loads(self._query_json)

    @property
    def guards(self):
        return json.loads(self._guards_json)

    def validate(self):
        self._owner._validate_selection(self)
        return True

    def __len__(self):
        return len(self._tags)

    def __iter__(self):
        return iter(self._tags)


def _selected_face_controls(model, values):
    """Validate model-owned selections before expanding bounded native controls."""
    result = []
    try:
        entries = sequence(() if values is None else values, "selected_face_sizes")
        for index, entry in enumerate(entries):
            try:
                entry = sequence(entry, f"selected_face_sizes[{index}]", maximum=3)
                if len(entry) not in (2, 3) or not isinstance(entry[0], FaceSelection):
                    raise ValueError("expected (FaceSelection, size) or (FaceSelection, size, growth)")
                selection = model._require_selection(entry[0])
                if len(result) + len(selection) > _MAX_FACE_CONTROLS:
                    raise SelectionError("control_budget", "Selected face sizing exceeds the control budget",
                                         {"control_limit": _MAX_FACE_CONTROLS,
                                          "requested_controls": len(result) + len(selection)}, stage="sizing")
                size = number(entry[1], "selected face size", minimum=0)
                if size == 0:
                    raise ValueError("selected face size must be positive")
                growth = entry[2] if len(entry) == 3 else None
                growth = None if growth is None else number(growth, "selected face growth", minimum=1)
                selection.validate()
                result.extend(_native.FaceSize(tag, size, growth=growth) for tag in selection.tags)
            except SelectionError:
                raise
            except (TypeError, ValueError, OverflowError) as error:
                raise SelectionError("invalid_argument", f"selected_face_sizes[{index}]: {error}",
                                     {"control_index": index}, stage="sizing") from error
        return result
    except SelectionError:
        raise
    except (TypeError, ValueError, OverflowError) as error:
        raise SelectionError("invalid_argument", str(error), stage="sizing") from error


class CadModel:
    """An imported native model with model-local selection and intent helpers."""

    __slots__ = ("_native", "_info_json", "_identity", "_sha256", "_kernel", "_tag_sets",
                 "_group_records", "_intent_revision")

    def __init__(self, native_model):
        self._native = native_model
        self._info_json = json_text(native_model.info)
        self._identity = integer(native_model.model_identity, "native model identity")
        self._sha256 = text(native_model.source_sha256, "source SHA-256", maximum=64)
        if len(self._sha256) != 64 or any(char not in "0123456789abcdef" for char in self._sha256):
            raise ValueError("Native source SHA-256 is invalid")
        self._kernel = text(native_model.kernel_version, "kernel version")
        self._tag_sets = {kind: _tags(getattr(native_model, attribute), kind)
                          for kind, attribute in (("faces", "face_tags"), ("bodies", "body_tags"),
                                                  ("edges", "edge_tags"), ("vertices", "vertex_tags"))}
        self._group_records = {}
        self._intent_revision = None

    @property
    def info(self):
        return json.loads(self._info_json)

    @property
    def source_sha256(self):
        return self._sha256

    @property
    def kernel_version(self):
        return self._kernel

    @property
    def model_identity(self):
        return self._identity

    @property
    def face_tags(self):
        return list(self._tag_sets["faces"])

    @property
    def body_tags(self):
        return list(self._tag_sets["bodies"])

    @property
    def edge_tags(self):
        return list(self._tag_sets["edges"])

    @property
    def vertex_tags(self):
        return list(self._tag_sets["vertices"])

    @property
    def groups(self):
        return {name: json.loads(record) for name, record in self._group_records.items()}

    def _assert_identity(self):
        if (self._native.model_identity != self._identity or self._native.source_sha256 != self._sha256 or
                self._native.kernel_version != self._kernel):
            raise SelectionError("model_changed", "The native model identity changed")

    def mesh(self, size, *, recovery="standard", assembly="preserve_parts",
             workers=1, validation_workers=1, max_nodes=None, max_cells=None,
             max_work=None, max_output_bytes=None, max_working_bytes=None,
             max_audit_bytes=None, max_seconds=None, face_sizes=None,
             selected_face_sizes=(), volume_controls=None, min_size=None, growth=1.3,
             native_sizing=None, warning_quality=None, minimum_quality=None,
             advancing_interior_front=None, construction_effort=_DEFAULT_CONSTRUCTION_EFFORT,
             max_surface_points=None, max_surface_chord_points=None,
             max_interior_candidates=None, max_interior_sites=None,
             resource_policy=_DEFAULT_RESOURCE_POLICY,
             max_chord_error=None, max_relative_volume_error=None, curvature_points_per_turn=None,
             initial_topology_quality=None):
        """Generate an owned, audited Tet4 mesh; lengths and sizes are in mm.

        ``face_sizes`` contains native ``FaceSize`` controls. Alternatively,
        ``selected_face_sizes`` accepts ``(FaceSelection, size)`` or
        ``(FaceSelection, size, growth)`` entries, with the same ownership and
        property guards as saved recipes. ``volume_controls`` contains BoxSize
        regions. Spatial controls require one solid and an explicit
        ``assembly="conforming"`` request; this method preserves the default
        ``preserve_parts`` assembly policy rather than changing it for controls.

        ``resource_policy`` selects normal (the default) or bounded admission.
        An omitted ``max_work`` uses the selected policy; an explicit positive
        limit always applies. Other omitted budgets retain native defaults.
        ``max_seconds`` is
        cooperative and starts after the native operation gate is acquired.
        Quality warnings do not replace validity checks. Native refusal reasons
        and hints remain available on ``MeshingError.report``.

        Construction effort controls bounded construction allowances. ``legacy``
        retains the established defaults. Optional positive integer caps override
        individual allowances within the existing global budgets; they do not
        select a recovery mode or weaken geometry and validity checks.
        """
        self._assert_identity()
        resource_policy = _resource_policy_value(resource_policy)
        fidelity = _fidelity_options(max_chord_error, max_relative_volume_error)
        curvature_points_per_turn = _curvature_value(curvature_points_per_turn)
        if initial_topology_quality is not None:
            boolean(initial_topology_quality, "initial_topology_quality")
        construction = _construction_options(construction_effort,
            max_surface_points=max_surface_points, max_surface_chord_points=max_surface_chord_points,
            max_interior_candidates=max_interior_candidates, max_interior_sites=max_interior_sites)
        selected = _selected_face_controls(self, selected_face_sizes)
        if selected:
            existing = sequence(() if face_sizes is None else face_sizes, "face_sizes")
            if len(existing) + len(selected) > _MAX_FACE_CONTROLS:
                raise SelectionError("control_budget", "Combined face sizing exceeds the control budget",
                                     {"control_limit": _MAX_FACE_CONTROLS,
                                      "requested_controls": len(existing) + len(selected)}, stage="sizing")
            face_sizes = existing + selected
        return self._native.mesh(size, recovery=recovery, assembly=assembly,
            workers=workers, validation_workers=validation_workers, max_nodes=max_nodes,
            max_cells=max_cells, max_work=max_work, max_output_bytes=max_output_bytes,
            max_working_bytes=max_working_bytes, max_audit_bytes=max_audit_bytes,
            max_seconds=max_seconds, face_sizes=face_sizes, volume_controls=volume_controls,
            min_size=min_size, growth=growth, native_sizing=native_sizing,
            warning_quality=warning_quality, minimum_quality=minimum_quality,
            advancing_interior_front=advancing_interior_front, resource_policy=resource_policy,
            curvature_points_per_turn=curvature_points_per_turn,
            initial_topology_quality=initial_topology_quality,
            **construction, **fidelity)

    def preview(self, selection=None, *, linear_deflection=.02, max_points=50000, max_triangles=100000):
        """Return owned, tagged display geometry for visual selection checks."""
        from .preview import preview_model
        return preview_model(self, selection, linear_deflection=linear_deflection,
                             max_points=max_points, max_triangles=max_triangles)

    def _records(self, kind, tags, *, topology=False):
        self._assert_identity()
        records = []
        payload_bytes = 2
        try:
            for start in range(0, len(tags), 1000):
                wanted = list(tags[start:start + 1000])
                received = self._native._inspect_entities(kind, wanted, topology=topology, max_uses=10000)
                if not isinstance(received, list):
                    raise ValueError("Native inspection did not return entity records")
                mapped = {integer(record["tag"], "native tag"): record for record in received}
                if len(mapped) != len(received) or set(mapped) != set(wanted):
                    raise ValueError("Native inspection omitted or duplicated requested entities")
                copied = json_copy([mapped[tag] for tag in wanted])
                payload_bytes += len(json_text(copied).encode("utf-8")) - 2 + bool(records and copied)
                if payload_bytes > MAX_JSON_BYTES:
                    raise SelectionError("inspection_budget", "Entity records exceed the 8 MiB metadata budget",
                                         {"kind": kind}, stage="inspection")
                records.extend(copied)
            return records
        except SelectionError:
            raise
        except (KeyError, TypeError, ValueError, OverflowError) as error:
            raise SelectionError("inspection_metadata", str(error), {"kind": kind}, stage="inspection") from error
        except (RuntimeError, OSError) as error:
            raise SelectionError("inspection_failed", str(error), {"kind": kind}, stage="inspection") from error

    def inspect(self, kind="faces", offset=0, limit=100, topology=False):
        """Return one bounded page of copied JSON metadata; IDs are model-local."""
        try:
            if kind not in _KINDS:
                raise ValueError("kind must be faces, bodies, edges or vertices")
            offset = integer(offset, "offset", minimum=0)
            limit = integer(limit, "limit")
            if limit > 1000:
                raise ValueError("limit must not exceed 1000")
            topology = boolean(topology, "topology")
            tags = self._tag_sets[kind]
            end = min(offset + limit, len(tags))
            result = {"schema": "paacs-cad-inspection-v1", "units": "mm",
                      "source_sha256": self._sha256, "kernel_version": self._kernel,
                      "kind": kind, "total": len(tags), "offset": offset,
                      "next_offset": end if end < len(tags) else None,
                      "entities": self._records(kind, tags[offset:end], topology=topology)}
            return json_copy(result)
        except SelectionError:
            raise
        except (TypeError, ValueError, OverflowError) as error:
            raise SelectionError("invalid_argument", str(error), stage="inspection") from error

    def _require_selection(self, selection):
        if not isinstance(selection, FaceSelection) or selection._owner is not self or selection._identity != self._identity:
            raise SelectionError("foreign_selection", "The selection belongs to a different imported model")
        self._assert_identity()
        return selection

    def _validate_selection(self, selection):
        self._require_selection(selection)
        known = set(self._tag_sets["faces"])
        if any(tag not in known for tag in selection.tags):
            raise SelectionError("unknown_face", "A selected face is absent from the model")
        if selection.body is not None and selection.body not in self._tag_sets["bodies"]:
            raise SelectionError("unknown_body", "The selected body is absent from the model")
        try:
            records = self._records("faces", selection.tags)
            if selection.body is not None and any(
                    not any(owner["body_tag"] == selection.body for owner in record["owners"])
                    for record in records):
                raise SelectionError("body_ownership", "A selected face does not belong to the selected body",
                                     {"body": selection.body, "face_tags": list(selection.tags)})
            actual = [_face_guard(record) for record in records]
        except SelectionError:
            raise
        except (KeyError, TypeError, ValueError, OverflowError) as error:
            raise SelectionError("guard_unavailable", "Selected face guards could not be evaluated",
                                 {"detail": str(error)}) from error
        if not close(actual, selection.guards):
            raise SelectionError("guard_mismatch", "Selected face properties changed")

    def _selection(self, tags, body, query, records=None, guards=None):
        if guards is None:
            if records is None:
                records = self._records("faces", tags)
            guards = [_face_guard(record) for record in records]
        return FaceSelection(self, self._identity, tuple(tags), body, json_text(query), json_text(guards))

    def _require_names(self):
        summary = self.info.get("representation_names", {})
        if not summary.get("requested", False):
            raise SelectionError("names_disabled", "STEP representation-name inspection was disabled")
        if not summary.get("complete", False):
            raise SelectionError("names_incomplete", "STEP representation-name mapping is incomplete",
                                 {"representation_names": summary})
        unmatched = summary.get("unmatched_by_dimension", {})
        face_gaps = unmatched.get("faces", unmatched.get("2", unmatched.get(2, 0)))
        if face_gaps:
            raise SelectionError("face_names_unmapped", "Some STEP face names could not be mapped to final placed faces",
                                 {"unmatched_face_names": face_gaps, "representation_names": summary})

    def _matches(self, record, query):
        tolerance = query["tolerance"]
        if "body" in query and not any(owner["body_tag"] == query["body"] for owner in record.get("owners", [])):
            return False
        if "surface" in query and record.get("surface") != query["surface"]:
            return False
        if "name" in query:
            status, names = record.get("name_status"), record.get("names", [])
            if status in ("disabled", "incomplete"):
                raise SelectionError("names_incomplete", "Face name metadata is incomplete", {"face_tag": record["tag"]})
            if query["name"] not in names:
                return False
            if status != "unique":
                raise SelectionError("ambiguous_name", "A matching face has conflicting representation names",
                                     {"face_tag": record["tag"], "names": names})
        for key in ("area", "radius"):
            if key in query:
                # A single analytic radius is exposed for cylinders and
                # spheres. Cones and tori have different radius parameters.
                if key == "radius" and record.get("surface") not in ("cylinder", "sphere"):
                    return False
                if not _metric_matches(_required(record, key, nested="geometry" if key == "radius" else None),
                                       query[key], tolerance):
                    return False
        if "centroid_box" in query:
            centroid = point(_required(record, "area_centroid"), "area_centroid")
            lo, hi = query["centroid_box"]
            if not all(a - tolerance <= value <= b + tolerance for a, value, b in zip(lo, centroid, hi)):
                return False
        if "adjacent_to" in query:
            adjacent = _required(record, "adjacent_faces")
            if not set(adjacent).intersection(query["adjacent_to"]):
                return False
        if "normal" in query:
            # A constant outward normal is only defined here for planar faces.
            if record.get("surface") != "plane":
                return False
            owners = record.get("owners", [])
            if "body" in query:
                owners = [owner for owner in owners if owner["body_tag"] == query["body"]]
            if len(owners) != 1:
                raise SelectionError("body_required", "A shared face needs a body scope for outward-normal selection",
                                     {"face_tag": record["tag"], "owners": [owner["body_tag"] for owner in owners]})
            outward = owners[0].get("outward_normal")
            if outward is None:
                raise SelectionError("metric_unavailable", "The face has no usable outward plane normal",
                                     {"face_tag": record["tag"]})
            outward = point(outward, "outward_normal")
            scale = max(abs(item) for item in outward)
            if scale == 0:
                raise SelectionError("metric_unavailable", "The outward plane normal has zero length")
            scaled = [item / scale for item in outward]
            length = math.hypot(*scaled)
            if sum(a * b / length for a, b in zip(scaled, query["normal"])) < 1 - tolerance:
                return False
        return True

    def select_faces(self, tags=None, name=None, group=None, surface=None, body=None,
                     area=None, radius=None, centroid_box=None, normal=None, adjacent_to=None,
                     expect=None, tolerance=1e-6, max_faces=10000):
        """Combine explicit predicates and enforce their expected cardinality."""
        try:
            if isinstance(tags, FaceSelection):
                tags = self._require_selection(tags).tags
            if isinstance(adjacent_to, FaceSelection):
                adjacent_to = self._require_selection(adjacent_to).tags
            query = _query_values(dict(tags=tags, name=name, group=group, surface=surface, body=body,
                                       area=area, radius=radius, centroid_box=centroid_box, normal=normal,
                                       adjacent_to=adjacent_to, expect=expect, tolerance=tolerance))
            if not set(query) - {"tolerance", "expect"}:
                raise ValueError("Provide face tags, a group, a name or a geometric predicate")
            maximum = integer(max_faces, "max_faces")
            candidates = self._tag_sets["faces"]
            for key in ("tags", "adjacent_to"):
                unknown = set(query.get(key, [])).difference(candidates)
                if unknown:
                    raise SelectionError("unknown_face", "The query contains unknown face IDs", {"tags": sorted(unknown)[:20]})
            if "body" in query and query["body"] not in self._tag_sets["bodies"]:
                raise SelectionError("unknown_body", "The query contains an unknown body ID", {"body": query["body"]})
            if "tags" in query:
                candidates = tuple(query["tags"])
            group_record = None
            if "group" in query:
                if query["group"] not in self._group_records:
                    raise SelectionError("unknown_group", "The query references an undefined companion group",
                                         {"group": query["group"]})
                group_record = json.loads(self._group_records[query["group"]])
                self._validate_group_record(group_record)
                saved_body = group_record["body"]
                if saved_body is not None:
                    if query.get("body", saved_body) != saved_body:
                        raise SelectionError("group_body_conflict", "The query contradicts the saved group's body",
                                             {"group": query["group"], "body": query["body"], "saved_body": saved_body})
                    query.setdefault("body", saved_body)
                group_tags = set(group_record["tags"])
                candidates = tuple(tag for tag in candidates if tag in group_tags)
            if len(candidates) > maximum:
                raise SelectionError("scan_budget", "The complete face scan exceeds max_faces",
                                     {"faces_to_scan": len(candidates), "max_faces": maximum})
            if "name" in query:
                self._require_names()
            selected, guards, guard_bytes = [], [], 2
            # Inspect the complete candidate set in bounded pages. A large
            # model can yield a small selection without retaining every face.
            for start in range(0, len(candidates), 1000):
                for record in self._records("faces", candidates[start:start + 1000]):
                    if self._matches(record, query):
                        guard = _face_guard(record)
                        guard_bytes += len(json_text(guard).encode("utf-8")) + bool(selected)
                        if guard_bytes > MAX_JSON_BYTES:
                            raise SelectionError("selection_budget", "Selected guards exceed the 8 MiB metadata budget")
                        selected.append(record["tag"])
                        guards.append(guard)
            if "expect" not in query:
                geometric = set(query) - {"tags", "group", "tolerance"}
                count = len(group_record["tags"]) if group_record is not None else 1 if geometric else len(candidates)
                query["expect"] = [count, count]
            lo, hi = query["expect"]
            if not lo <= len(selected) <= hi:
                raise SelectionError("cardinality", "Face selection did not meet its expected count",
                                     {"expected": [lo, hi], "matched": len(selected),
                                      "candidate_tags": selected[:20], "query": query})
            return self._selection(selected, query.get("body"), query, guards=guards)
        except SelectionError:
            raise
        except (KeyError, TypeError, ValueError, OverflowError) as error:
            raise SelectionError("invalid_argument", str(error)) from error

    def _validate_group_record(self, record):
        selection = FaceSelection(self, self._identity, tuple(record["tags"]), record["body"],
                                  json_text(record["query"]), json_text(record["guards"]))
        self._validate_selection(selection)

    def define_group(self, name, selection, group_id=None):
        from ._cad_intent import define_group
        return define_group(self, name, selection, group_id)

    def save_intent(self, path, revision=None):
        from ._cad_intent import save_intent
        return save_intent(self, path, revision)

    def load_intent(self, path, rebind=False):
        from ._cad_intent import load_intent
        return load_intent(self, path, rebind)

    def _source_snapshot(self):
        """Return exact import-time bytes or refuse a changed/deleted source."""
        self._assert_identity()
        info = self.info
        expected = integer(info["file_bytes"], "source file bytes", minimum=0)
        payload = self._native._source_bytes
        if payload is None:
            try:
                with Path(info["source_path"]).open("rb") as stream:
                    payload = stream.read(expected + 1)
            except (OSError, OverflowError) as error:
                raise IntentError("source_unavailable", "The imported STEP bytes are unavailable",
                                  {"detail": str(error)}) from error
        if (not isinstance(payload, bytes) or len(payload) != expected or
                hashlib.sha256(payload).hexdigest() != self._sha256):
            raise IntentError("source_changed", "The STEP source no longer matches its import-time hash")
        return payload


def read_step(path, *, max_file_bytes=None, intent=None, names=True, retain_source=False):
    """Import STEP and optionally apply a strictly matching companion manifest."""
    names, retain_source = boolean(names, "names"), boolean(retain_source, "retain_source")
    model = CadModel(_native.read_step(path, max_file_bytes=max_file_bytes,
                                      names=names, retain_source=retain_source))
    if intent is not None:
        model.load_intent(intent)
    return model
