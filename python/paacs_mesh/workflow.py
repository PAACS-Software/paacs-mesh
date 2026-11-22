"""Reproducible, bounded meshing attempts and separate local publication.

Workflow deadlines are cooperative. They include Python validation, native
operation queueing and every retry; OCCT calls and queue waits can overrun them.
Publication is separate and reports its own time. No validity check is disabled.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from dataclasses import dataclass
import errno
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time

from . import _native
from ._cad_common import (
    MAX_JSON_BYTES, INTENT_SCHEMA, TAG_SCHEME, IntentError, SelectionError, atomic_json, boolean,
    bounds, integer, json_copy, json_text, load_json, number, sequence, text,
)
from ._cad_intent import _group
from .cad import CadModel, FaceSelection, read_step


RECIPE_SCHEMA = "paacs-mesh-recipe-v1"
DEFAULT_SOURCE_BYTES = 64 * 1024 * 1024
_BUDGETS = ("max_nodes", "max_cells", "max_work", "max_output_bytes",
            "max_working_bytes", "max_audit_bytes")
_QUALITY_FIELDS = ("warning_quality", "minimum_quality")
_PLACEMENT_FIELDS = ("advancing_interior_front",)
_OPTION_FIELDS = {"size", "assembly", "workers", "validation_workers", "face_sizes",
                  "volume_controls", "min_size", "growth", "native_sizing", *_QUALITY_FIELDS, *_PLACEMENT_FIELDS, *_BUDGETS}
_FATAL_REASONS = {"invalid_request", "unsupported_body", "unresolved_body_interface",
                  "unsupported_control", "budget_exceeded", "cancelled", "allocation_failed",
                  "internal_error", "publication_checks_failed"}
_ADVANCED_DIAGNOSTICS = {"spaced-surface-constraint-seed-refused", "native-plc-seed-candidate-refused",
                         "native-plc-cavity-material-refusal", "native-curvature-material-refused"}
_FINER_DIAGNOSTICS = _ADVANCED_DIAGNOSTICS | {
    "native-boundary-volume-error", "native-plc-source-volume-error", "native-plc-final-volume-error",
    "native-plc-adaptive-facet-error", "native-plc-final-centroid-outside",
    "native-plc-final-interior-node-outside", "native_tet_centroid_outside",
    "native-plc-final-native-validation", "spaced-surface-quality-requirement",
}
_FINER_SURFACES = {"missing_triangulation", "nonconforming_edge", "invalid_surface"}
_REPORT_RESERVE = 65536


def _small(value):
    if isinstance(value, str):
        return value[:256].encode("utf-8", errors="backslashreplace").decode("utf-8")
    if value is None or isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value if -(2 ** 63) <= value < 2 ** 63 else None
    if isinstance(value, float) and math.isfinite(value):
        return value
    return None


def _native_summary(report):
    if not isinstance(report, dict):
        return {"available": False}
    result = {key: _small(report.get(key)) for key in ("status", "reason", "mesh_published")}
    for stage in ("audit", "quality"):
        data = report.get(stage)
        if isinstance(data, dict):
            result[stage] = {key: _small(data.get(key)) for key in ("complete", "valid")}
    bodies = report.get("bodies")
    if isinstance(bodies, list):
        result["body_count"] = len(bodies)
        result["bodies"] = [{key: _small(body.get(key)) for key in
                              ("solid_tag", "reason", "surface_reason", "volume_reason")}
                             for body in bodies[:8] if isinstance(body, dict)]
        result["body_summaries_complete"] = len(bodies) <= 8
    return result


def _error_summary(report, original_reason):
    """Last-resort bounded fallback; never present summaries as complete data."""
    result = {"schema": "paacs-mesh-workflow-v1", "status": "failed", "diagnostics_complete": False,
              "original_reason": _small(original_reason), "report_omission": "serialization_or_metadata_budget",
              "attempts": []}
    if not isinstance(report, dict):
        return result
    result["seconds"] = _small(report.get("seconds"))
    if "native_report" in report:
        result["native_report_summary"] = _native_summary(report["native_report"])
    attempts = report.get("attempts", [])
    if isinstance(attempts, list):
        result["attempt_count"] = len(attempts)
        for attempt in attempts[:8]:
            if not isinstance(attempt, dict):
                continue
            item = {key: _small(attempt.get(key)) for key in
                    ("number", "status", "reason", "stage", "started_seconds", "seconds")}
            item.update(options_complete=False, native_report_complete=False,
                        native_report_summary=_native_summary(attempt.get("native_report")))
            options = attempt.get("options", {})
            if isinstance(options, dict):
                item["options_summary"] = {key: _small(options.get(key)) for key in
                                            ("size", "recovery", "assembly", *_QUALITY_FIELDS, *_BUDGETS)}
            result["attempts"].append(item)
    return result


class WorkflowError(RuntimeError):
    """A recipe, bounded run or publication failed without publishing a mesh."""

    def __init__(self, reason, message, report=None, *, stage="workflow"):
        super().__init__(f"{stage}: {message}")
        self.stage, self.reason = stage, reason
        try:
            self.report = {} if report is None else json_copy(report)
        except (TypeError, ValueError, OverflowError) as error:
            self.reason = "report_budget" if "budget" in str(error) else "invalid_report"
            self.report = json_copy(_error_summary(report, reason))
            self.args = (f"{stage}: {message}; diagnostic report omitted ({self.reason})",)


@dataclass(frozen=True)
class WorkflowResult:
    """An owned successful mesh with independent recipe and attempt reports."""

    # Explicit slots keep the original class identity used by the frozen
    # dataclass setter on CPython 3.12, including assignments to properties.
    __slots__ = ("_model", "_mesh", "_recipe_json", "_report_json")

    _model: CadModel
    _mesh: object
    _recipe_json: str
    _report_json: str

    @property
    def mesh(self):
        return self._mesh

    @property
    def source_sha256(self):
        return self.recipe["source"]["sha256"]

    @property
    def recipe(self):
        return json.loads(self._recipe_json)

    @property
    def report(self):
        return json.loads(self._report_json)


def _report_room(attempts, decisions, candidate):
    try:
        payload = json_text({"attempts": [*attempts, candidate], "decisions": decisions})
        return len(payload.encode("utf-8")) <= MAX_JSON_BYTES - _REPORT_RESERVE
    except (TypeError, ValueError, OverflowError):
        return False


def _report_failure(attempts, decisions, started, source_hash, *, detail):
    return {"schema": "paacs-mesh-workflow-v1", "status": "failed", "diagnostics_complete": False,
            "attempts": attempts, "decisions": decisions, "seconds": time.perf_counter() - started,
            "source_sha256": source_hash, "report_budget_bytes": MAX_JSON_BYTES,
            "report_reserved_bytes": _REPORT_RESERVE, "cooperative_deadline": True,
            "budgets_apply_per_attempt": True,
            "report_budget_detail": detail, "publication_included": False}


def _record_native_attempt(attempts, decisions, raw_report, *, actual, trigger, status,
                           reason, stage, began, started, source_hash):
    entry = {"number": len(attempts) + 1, "status": status, "reason": reason, "stage": stage,
             "options": actual, "options_complete": True, "trigger": trigger,
             "started_seconds": began - started, "native_report": None,
             "native_report_complete": False}
    size, issue = None, None
    try:
        encoded = json.dumps(raw_report, ensure_ascii=False, allow_nan=False,
                             separators=(",", ":"), sort_keys=True)
        size = len(encoded.encode("utf-8"))
        if size > MAX_JSON_BYTES:
            issue = "single_native_report_exceeds_budget"
        else:
            full = {**entry, "native_report": json.loads(encoded), "native_report_complete": True,
                    "native_report_bytes": size, "seconds": time.perf_counter() - began}
            if _report_room(attempts, decisions, full):
                full["seconds"] = time.perf_counter() - began
                attempts.append(full)
                return full["native_report"]
            issue = "aggregate_report_exceeds_budget"
    except (TypeError, ValueError, OverflowError, RecursionError):
        issue = "native_report_is_not_finite_json"
    entry.update(native_report_bytes=size, native_report_omission=issue,
                 native_report_summary=_native_summary(raw_report), seconds=time.perf_counter() - began)
    attempts.append(entry)
    reason = "invalid_report" if issue == "native_report_is_not_finite_json" else "report_budget"
    raise WorkflowError(reason, "Native attempt recorded with an explicitly omitted diagnostic report",
                        _report_failure(attempts, decisions, started, source_hash, detail=issue), stage="diagnostics")


def _positive(value, name):
    value = number(value, name, minimum=0)
    if value == 0:
        raise ValueError(f"{name} must be positive")
    return value


def _source(model):
    return {"sha256": model.source_sha256, "bytes": model.info["file_bytes"],
            "name": Path(model.info["source_path"]).name, "units": "mm"}


def _intent(model):
    return {"schema": INTENT_SCHEMA, "source": _source(model), "kernel_version": model.kernel_version,
            "tag_scheme": TAG_SCHEME, "revision": model._intent_revision,
            "groups": [{"name": name, **record} for name, record in sorted(model.groups.items())]}


def _options(values):
    # Recipe v1 predates optional quality and initial-placement choices.
    # Missing optional choices retain automatic defaults; unknown fields fail.
    required = _OPTION_FIELDS.difference(_QUALITY_FIELDS, _PLACEMENT_FIELDS)
    if not isinstance(values, dict) or not required <= values.keys() or values.keys() - _OPTION_FIELDS:
        raise ValueError("Recipe meshing options have unsupported or missing fields")
    result = {"size": _positive(values["size"], "size"), "assembly": values["assembly"],
              "workers": integer(values["workers"], "workers"),
              "validation_workers": integer(values["validation_workers"], "validation_workers"),
              "growth": number(values["growth"], "growth", minimum=1),
              "native_sizing": values["native_sizing"], "min_size": values["min_size"]}
    result["advancing_interior_front"] = values.get("advancing_interior_front")
    if result["advancing_interior_front"] is not None:
        boolean(result["advancing_interior_front"], "advancing_interior_front")
    if result["assembly"] not in ("preserve_parts", "conforming"):
        raise ValueError("assembly must be preserve_parts or conforming")
    if result["workers"] > 1024 or result["validation_workers"] > 8:
        raise ValueError("workers must be 1..1024 and validation_workers 1..8")
    if result["native_sizing"] is not None:
        boolean(result["native_sizing"], "native_sizing")
    if result["min_size"] is not None:
        result["min_size"] = number(result["min_size"], "min_size", minimum=0)
        if result["min_size"] > result["size"]:
            raise ValueError("min_size must not exceed global size")
    for key in _QUALITY_FIELDS:
        value = values.get(key)
        result[key] = None if value is None else number(value, key, minimum=0)
        if result[key] is not None and result[key] > 1:
            raise ValueError(f"{key} must be between zero and one")
    for key in _BUDGETS:
        result[key] = None if values[key] is None else integer(values[key], key)
    result["face_sizes"], result["volume_controls"] = [], []
    seen = set()
    for control in sequence(values["face_sizes"], "face sizes"):
        if not isinstance(control, dict) or set(control) != {"tag", "size", "growth"}:
            raise ValueError("Invalid recipe face control")
        tag = integer(control["tag"], "face tag")
        if tag in seen:
            raise ValueError("Recipe contains duplicate face sizing controls")
        seen.add(tag)
        result["face_sizes"].append({"tag": tag, "size": _positive(control["size"], "face size"),
                                     "growth": None if control["growth"] is None else number(control["growth"], "face growth", minimum=1)})
    for control in sequence(values["volume_controls"], "volume controls"):
        if not isinstance(control, dict) or set(control) != {"lower", "upper", "size", "growth"}:
            raise ValueError("Invalid recipe box control")
        lo, hi = bounds([control["lower"], control["upper"]], "box")
        if any(a >= b for a, b in zip(lo, hi)):
            raise ValueError("Box lower coordinates must be strictly below upper coordinates")
        result["volume_controls"].append({"lower": lo, "upper": hi, "size": _positive(control["size"], "box size"),
                                          "growth": None if control["growth"] is None else number(control["growth"], "box growth", minimum=1)})
    if result["volume_controls"] and result["native_sizing"] is False:
        raise ValueError("Box controls require native spatial sizing")
    return result


def _policy(values):
    if not isinstance(values, dict) or set(values) != {"max_attempts", "max_seconds", "attempt_seconds", "allow_advanced", "finer_factors"}:
        raise ValueError("Invalid workflow retry policy")
    result = {"max_attempts": integer(values["max_attempts"], "max_attempts"),
              "max_seconds": _positive(values["max_seconds"], "max_seconds"),
              "attempt_seconds": values["attempt_seconds"],
              "allow_advanced": boolean(values["allow_advanced"], "allow_advanced"), "finer_factors": []}
    if result["max_attempts"] > 8 or result["max_seconds"] > 86400:
        raise ValueError("A workflow supports at most 8 attempts and 86400 cooperative seconds")
    if result["attempt_seconds"] is not None:
        result["attempt_seconds"] = _positive(result["attempt_seconds"], "attempt_seconds")
        if result["attempt_seconds"] > result["max_seconds"]:
            raise ValueError("attempt_seconds must not exceed the total workflow budget")
    previous = 1.0
    for factor in sequence(values["finer_factors"], "finer factors", maximum=7):
        factor = _positive(factor, "finer factor")
        if factor >= previous:
            raise ValueError("finer_factors must be strictly decreasing and below 1")
        result["finer_factors"].append(factor)
        previous = factor
    return result


def create_recipe(model, size, *, face_sizes=(), selected_face_sizes=(), volume_controls=(),
                  min_size=None, growth=1.3, native_sizing=None, assembly="preserve_parts",
                  warning_quality=None, minimum_quality=None, advancing_interior_front=None,
                  workers=1, validation_workers=1, max_nodes=None, max_cells=None,
                  max_work=None, max_output_bytes=None, max_working_bytes=None, max_audit_bytes=None,
                  allow_advanced=False, finer_factors=(), max_attempts=1, max_seconds=90,
                  attempt_seconds=None):
    """Capture controls, source identity and intent without changing model policy.

    selected_face_sizes contains (FaceSelection, size) or (FaceSelection, size,
    growth) entries. Retry sizes are explicit fractions of the original global
    size; existing face/box/minimum controls and resource budgets are unchanged.
    Optional warning and strict quality thresholds remain fixed across retries.
    """
    try:
        if not isinstance(model, CadModel):
            raise ValueError("model must be a CadModel")
        model._assert_identity()
        faces, boxes = [], []
        for control in sequence(face_sizes, "face_sizes"):
            if not isinstance(control, _native.FaceSize):
                raise ValueError("face_sizes must contain FaceSize controls")
            faces.append({"tag": control.tag, "size": control.size, "growth": control.growth})
        for entry in sequence(selected_face_sizes, "selected_face_sizes"):
            entry = sequence(entry, "selected face size", maximum=3)
            if len(entry) not in (2, 3) or not isinstance(entry[0], FaceSelection):
                raise ValueError("selected_face_sizes entries must contain a FaceSelection and size")
            selection = model._require_selection(entry[0])
            selection.validate()
            for tag in selection.tags:
                faces.append({"tag": tag, "size": entry[1], "growth": entry[2] if len(entry) == 3 else None})
        for control in sequence(volume_controls, "volume_controls"):
            if not isinstance(control, _native.BoxSize):
                raise ValueError("volume_controls must contain BoxSize controls")
            boxes.append({"lower": list(control.lower), "upper": list(control.upper),
                          "size": control.size, "growth": control.growth})
        options = _options(dict(size=size, assembly=assembly, workers=workers,
            validation_workers=validation_workers, face_sizes=faces, volume_controls=boxes,
            min_size=min_size, growth=growth, native_sizing=native_sizing,
            warning_quality=warning_quality, minimum_quality=minimum_quality,
            advancing_interior_front=advancing_interior_front,
            max_nodes=max_nodes, max_cells=max_cells, max_work=max_work,
            max_output_bytes=max_output_bytes, max_working_bytes=max_working_bytes, max_audit_bytes=max_audit_bytes))
        policy = _policy(dict(max_attempts=max_attempts, max_seconds=max_seconds, attempt_seconds=attempt_seconds,
                              allow_advanced=allow_advanced, finer_factors=list(sequence(finer_factors, "finer_factors", maximum=7))))
        tags = sorted(control["tag"] for control in options["face_sizes"])
        selected = model.select_faces(tags=tags, expect=[len(tags), len(tags)]) if tags else None
        face_guards = None if selected is None else {"tags": tags, "body": None,
                      "query": selected.query, "guards": selected.guards}
        recipe = {"schema": RECIPE_SCHEMA, "source": _source(model),
                  "kernel_version": model.kernel_version, "tag_scheme": TAG_SCHEME,
                  "engine_version": str(_native.engine_version), "intent": _intent(model),
                  "face_guards": face_guards, "options": options, "retry": policy}
        return _validate_recipe(model, recipe)
    except WorkflowError:
        raise
    except (SelectionError, IntentError) as error:
        raise WorkflowError(error.reason, str(error), error.report, stage="recipe") from error
    except (TypeError, ValueError, OverflowError) as error:
        raise WorkflowError("invalid_recipe", str(error), stage="recipe") from error


def _validate_recipe(model, recipe):
    if not isinstance(model, CadModel):
        raise ValueError("model must be a CadModel")
    recipe = json_copy(recipe)
    fields = {"schema", "source", "kernel_version", "tag_scheme", "engine_version", "intent", "face_guards", "options", "retry"}
    if not isinstance(recipe, dict) or set(recipe) != fields:
        raise ValueError("Invalid mesh recipe fields")
    if recipe["schema"] != RECIPE_SCHEMA:
        raise WorkflowError("unknown_schema", "Unsupported mesh recipe schema", stage="recipe")
    model._assert_identity()
    expected = _source(model)
    source = recipe["source"]
    if not isinstance(source, dict) or set(source) != set(expected):
        raise ValueError("Invalid recipe source metadata")
    text(source["name"], "source name", maximum=1024)
    integer(source["bytes"], "source bytes", minimum=0)
    if any(source[key] != expected[key] for key in ("sha256", "bytes", "units")):
        raise WorkflowError("source_mismatch", "Recipe source does not match the imported STEP", stage="recipe")
    if recipe["kernel_version"] != model.kernel_version:
        raise WorkflowError("kernel_mismatch", "Recipe kernel does not match this import", stage="recipe")
    if recipe["tag_scheme"] != TAG_SCHEME:
        raise WorkflowError("tag_scheme_mismatch", "Recipe tag scheme does not match this import", stage="recipe")
    if recipe["engine_version"] != str(_native.engine_version):
        raise WorkflowError("engine_mismatch", "Recipe native engine version does not match", stage="recipe")
    recipe["options"], recipe["retry"] = _options(recipe["options"]), _policy(recipe["retry"])
    intent = recipe["intent"]
    if (not isinstance(intent, dict) or set(intent) != {"schema", "source", "kernel_version", "tag_scheme", "revision", "groups"}
            or intent["schema"] != INTENT_SCHEMA or intent["source"] != source
            or intent["kernel_version"] != model.kernel_version or intent["tag_scheme"] != TAG_SCHEME):
        raise ValueError("Recipe contains inconsistent CAD intent")
    if intent["revision"] is not None:
        text(intent["revision"], "intent revision")
    names, identifiers = set(), set()
    for value in sequence(intent["groups"], "intent groups"):
        name, record = _group(value)
        if name in names or record["id"] in identifiers:
            raise ValueError("Recipe contains duplicate companion groups")
        names.add(name)
        identifiers.add(record["id"])
        model._validate_group_record(record)
    face_tags = sorted(control["tag"] for control in recipe["options"]["face_sizes"])
    guards = recipe["face_guards"]
    if face_tags:
        if (not isinstance(guards, dict) or set(guards) != {"tags", "body", "query", "guards"}
                or guards["tags"] != face_tags or guards["body"] is not None):
            raise ValueError("Recipe face controls lack matching property guards")
        # The manifest parser validates all guard fields without mutating model groups.
        _, normalized = _group({"name": "workflow controls", "id": "workflow controls", **guards})
        model._validate_group_record(normalized)
    elif guards is not None:
        raise ValueError("Recipe has unexpected face-control guards")
    return recipe


def save_recipe(model, recipe, path):
    """Atomically save a recipe, protecting the original STEP path."""
    try:
        recipe = _validate_recipe(model, recipe)
        return atomic_json(path, recipe, model.info["source_path"])
    except WorkflowError:
        raise
    except (SelectionError, IntentError) as error:
        raise WorkflowError(error.reason, str(error), error.report, stage="recipe") from error
    except (OSError, TypeError, ValueError, OverflowError) as error:
        raise WorkflowError("recipe_save_failed", str(error), stage="recipe") from error


def load_recipe(model, path):
    """Read bounded JSON and validate it against an already imported model."""
    try:
        return _validate_recipe(model, load_json(path))
    except WorkflowError:
        raise
    except (SelectionError, IntentError) as error:
        raise WorkflowError(error.reason, str(error), error.report, stage="recipe") from error
    except (OSError, TypeError, ValueError, OverflowError) as error:
        raise WorkflowError("invalid_recipe", str(error), stage="recipe") from error


def _diagnostics(report):
    bodies = report.get("bodies", [])
    return ({body.get("volume_reason") for body in bodies}, {body.get("surface_reason") for body in bodies})


def _retry_allowed(reason, report, route):
    # The aggregate reason is the failure authority. Body diagnostics only
    # qualify whether an explicitly opted-in alternative is useful.
    if reason in _FATAL_REASONS:
        return False
    volume, surface = _diagnostics(report)
    if route == "advanced":
        return reason in ("meshing_failed", "geometry_error") and bool(volume & _ADVANCED_DIAGNOSTICS)
    return (reason == "size_requirement_not_met" or
            reason in ("meshing_failed", "geometry_error", "quality_requirement_not_met") and
            bool(volume & _FINER_DIAGNOSTICS or surface & _FINER_SURFACES))


def _native_options(options, size, recovery, seconds):
    result = {key: options[key] for key in ("assembly", "workers", "validation_workers", "min_size", "growth", "native_sizing", *_QUALITY_FIELDS, *_PLACEMENT_FIELDS, *_BUDGETS)}
    result.update(size=size, recovery=recovery, max_seconds=seconds,
                  face_sizes=[_native.FaceSize(control["tag"], control["size"], growth=control["growth"])
                              for control in options["face_sizes"]],
                  volume_controls=[_native.BoxSize(control["lower"], control["upper"], control["size"], growth=control["growth"])
                                   for control in options["volume_controls"]])
    return result


def _attempt_options(options, size, recovery, seconds):
    result = json_copy(options)
    result.update(size=size, recovery=recovery, max_seconds=seconds)
    return result


def _run(model, recipe, started, import_seconds=0):
    attempts, decisions = [], []
    try:
        recipe = _validate_recipe(model, recipe)
        options, policy = recipe["options"], recipe["retry"]
        size, recovery, finer_index = options["size"], "standard", 0
        advanced_used = False
        trigger = {"route": "standard", "reason": "initial"}
        while len(attempts) < policy["max_attempts"]:
            remaining = policy["max_seconds"] - (time.perf_counter() - started)
            if remaining <= 0:
                break
            seconds = min(remaining, policy["attempt_seconds"] or remaining)
            actual = _attempt_options(options, size, recovery, seconds)
            call_options = _native_options(options, size, recovery, seconds)
            # Control conversion and report preparation consume the same wall
            # budget. Do not begin another native call after that work expires it.
            remaining = policy["max_seconds"] - (time.perf_counter() - started)
            if remaining <= 0:
                break
            seconds = min(seconds, remaining)
            actual["max_seconds"] = call_options["max_seconds"] = seconds
            pending = {"number": len(attempts) + 1, "options": actual, "trigger": trigger,
                       "native_report": None, "native_report_complete": False, "seconds": 0}
            if not _report_room(attempts, decisions, pending):
                report = _report_failure(attempts, decisions, started, model.source_sha256,
                                         detail="next_attempt_options_exceed_budget")
                report["next_attempt_not_started"] = len(attempts) + 1
                raise WorkflowError("report_budget", "Another attempt cannot be fully recorded within the report budget",
                                    report, stage="diagnostics")
            remaining = policy["max_seconds"] - (time.perf_counter() - started)
            if remaining <= 0:
                break
            actual["max_seconds"] = call_options["max_seconds"] = min(seconds, remaining)
            began = time.perf_counter()
            try:
                mesh = model.mesh(**call_options)
                raw_report = mesh.report
                # A workflow never upgrades a refused or incomplete native result.
                valid = (isinstance(raw_report, dict) and raw_report.get("mesh_published") is True and raw_report.get("reason") == "none"
                         and raw_report.get("audit", {}).get("complete") is True
                         and raw_report.get("audit", {}).get("valid") is True
                         and raw_report.get("quality", {}).get("complete") is True
                         and raw_report.get("quality", {}).get("valid") is True)
                native_report = _record_native_attempt(attempts, decisions, raw_report, actual=actual,
                    trigger=trigger, status="succeeded" if valid else "failed",
                    reason="none" if valid else "publication_checks_failed", stage="mesh", began=began,
                    started=started, source_hash=model.source_sha256)
                if not valid:
                    raise WorkflowError("publication_checks_failed", "Native mesh did not retain complete validity checks",
                                        {"attempts": attempts, "seconds": time.perf_counter() - started})
            except _native.MeshingError as error:
                native_report = _record_native_attempt(attempts, decisions, error.report, actual=actual,
                    trigger=trigger, status="failed", reason=error.reason, stage=error.stage, began=began,
                    started=started, source_hash=model.source_sha256)
                next_route = None
                if policy["allow_advanced"] and not advanced_used and _retry_allowed(error.reason, native_report, "advanced"):
                    advanced_used, recovery = True, "advanced"
                    next_route = "advanced"
                elif _retry_allowed(error.reason, native_report, "finer"):
                    while finer_index < len(policy["finer_factors"]):
                        factor = policy["finer_factors"][finer_index]
                        finer_index += 1
                        candidate = options["size"] * factor
                        if candidate <= 0:
                            decisions.append({"route": "finer", "factor": factor, "status": "skipped_unrepresentable_size"})
                            continue
                        if candidate < (options["min_size"] or 0):
                            decisions.append({"route": "finer", "factor": factor, "status": "skipped_minimum_size"})
                            continue
                        size, next_route = candidate, "finer"
                        break
                if next_route is None:
                    break
                trigger = {"route": next_route, "reason": error.reason,
                           "volume_diagnostics": sorted(_diagnostics(native_report)[0] & _FINER_DIAGNOSTICS)}
                continue
            except (MemoryError, RuntimeError, OSError, TypeError, ValueError, OverflowError) as error:
                if isinstance(error, WorkflowError):
                    raise
                reason = "allocation_failed" if isinstance(error, MemoryError) else "native_exception"
                attempts.append({"number": len(attempts) + 1, "status": "failed", "reason": reason,
                                 "stage": "mesh", "options": actual, "trigger": trigger,
                                 "started_seconds": began - started, "seconds": time.perf_counter() - began,
                                 "native_report": None, "native_report_complete": False,
                                 "native_report_omission": "native_exception", "detail": _small(str(error))})
                raise WorkflowError(reason, "A native call did not return its usual diagnostic report",
                                    {"attempts": attempts, "seconds": time.perf_counter() - started}) from error
            total = time.perf_counter() - started
            report = {"schema": "paacs-mesh-workflow-v1", "status": "succeeded", "attempts": attempts,
                      "diagnostics_complete": True,
                      "decisions": decisions, "source_sha256": model.source_sha256,
                      "seconds": total, "import_seconds": import_seconds,
                      "overhead_seconds": max(0, total - sum(attempt["seconds"] for attempt in attempts) - import_seconds),
                      "cooperative_deadline": True, "time_budget_overrun": total > policy["max_seconds"],
                      "budgets_apply_per_attempt": True, "publication_included": False}
            return WorkflowResult(model, mesh, json_text(recipe), json_text(report))
        total = time.perf_counter() - started
        reason = "time_budget" if total >= policy["max_seconds"] else attempts[-1]["reason"] if attempts else "no_attempt"
        raise WorkflowError(reason, "The bounded workflow did not produce a valid mesh",
                            {"schema": "paacs-mesh-workflow-v1", "status": "failed", "attempts": attempts,
                             "diagnostics_complete": True,
                             "decisions": decisions, "seconds": total, "import_seconds": import_seconds,
                             "overhead_seconds": max(0, total - sum(attempt["seconds"] for attempt in attempts) - import_seconds),
                             "cooperative_deadline": True, "budgets_apply_per_attempt": True,
                             "publication_included": False,
                             "stopped_by": "time_budget" if total >= policy["max_seconds"] else
                                           "attempt_limit" if len(attempts) >= policy["max_attempts"] else "no_eligible_retry",
                             "source_sha256": model.source_sha256})
    except WorkflowError:
        raise
    except (SelectionError, IntentError) as error:
        raise WorkflowError(error.reason, str(error), {"attempts": attempts, "detail": error.report}) from error
    except (TypeError, ValueError, OverflowError) as error:
        raise WorkflowError("invalid_recipe", str(error), {"attempts": attempts}) from error


def run_recipe(model, recipe):
    """Run an existing model with total cooperative time and attempt limits."""
    return _run(model, recipe, time.perf_counter())


def run_saved_recipe(path, source_path, *, max_file_bytes=DEFAULT_SOURCE_BYTES, names=True, retain_source=False):
    """Import an explicitly supplied local source and run its saved recipe.

    The manifest's display name is never used to locate or open a source file.
    Import and recipe reading count toward total time but cannot be interrupted.
    """
    started = time.perf_counter()
    try:
        recipe = load_json(path)
        before_import = time.perf_counter()
        model = read_step(source_path, max_file_bytes=max_file_bytes, names=names, retain_source=retain_source)
        return _run(model, recipe, started, time.perf_counter() - before_import)
    except WorkflowError:
        raise
    except _native.CadImportError as error:
        raise WorkflowError(error.reason, str(error), {"attempts": [], "seconds": time.perf_counter() - started,
                                                      "native_report": error.report}, stage="import") from error
    except (OSError, TypeError, ValueError, OverflowError) as error:
        raise WorkflowError("recipe_read_failed", str(error), {"attempts": []}, stage="recipe") from error


def _rename_new_directory(source, destination):
    """Atomically rename without replacing even an existing empty directory."""
    if sys.platform == "win32":
        os.rename(source, destination)  # Windows refuses an existing destination.
        return
    import ctypes
    library = ctypes.CDLL(None, use_errno=True)
    if sys.platform.startswith("linux"):
        function = getattr(library, "renameat2", None)
        if function is None:
            raise OSError(errno.ENOTSUP, "Exclusive atomic directory rename is unavailable")
        function.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        function.restype = ctypes.c_int
        status = function(-100, os.fsencode(source), -100, os.fsencode(destination), 1)  # RENAME_NOREPLACE
    elif sys.platform == "darwin":
        function = getattr(library, "renamex_np", None)
        if function is None:
            raise OSError(errno.ENOTSUP, "Exclusive atomic directory rename is unavailable")
        function.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        function.restype = ctypes.c_int
        status = function(os.fsencode(source), os.fsencode(destination), 4)  # RENAME_EXCL
    else:
        raise OSError(errno.ENOTSUP, "Exclusive atomic directory rename is unsupported on this platform")
    if status:
        code = ctypes.get_errno()
        raise OSError(code, os.strerror(code), str(destination))


def publish_result(result, directory, *, max_source_bytes=DEFAULT_SOURCE_BYTES):
    """Publish a new local directory after all exports and metadata succeed.

    Existing directories, files and symlinks are refused. A failed publication
    leaves the owned result usable. Native export validity checks remain active.
    Publication uses a sibling staging directory and exclusive atomic rename;
    filesystems lacking that operation fail without a non-atomic fallback.
    """
    started, staging, parent = time.perf_counter(), None, None
    try:
        if not isinstance(result, WorkflowResult):
            raise ValueError("result must be a WorkflowResult")
        maximum = integer(max_source_bytes, "max_source_bytes")
        if not str(directory) or "\0" in str(directory):
            raise ValueError("Publication directory must be nonempty and contain no NUL")
        destination = Path(directory).absolute()
        if os.path.lexists(destination):
            raise WorkflowError("output_exists", "Publication never replaces an existing output", stage="publication")
        parent = destination.parent.resolve(strict=True)
        destination = parent / destination.name
        if not parent.is_dir():
            raise ValueError("Publication parent must be an existing directory")
        model, recipe = result._model, result.recipe
        _validate_recipe(model, recipe)
        if model.info["file_bytes"] > maximum:
            raise WorkflowError("source_budget", "Source embedding exceeds max_source_bytes", stage="publication")
        source = model._source_snapshot()
        if len(source) > maximum:
            raise WorkflowError("source_budget", "Source embedding exceeds max_source_bytes", stage="publication")
        staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.paacs-", dir=parent)).resolve()
        if staging.parent != parent:
            raise WorkflowError("staging_path", "Staging directory escaped its publication parent", stage="publication")
        result.mesh.write_archive(staging / "mesh.paacsmesh")
        result.mesh.write_vtk(staging / "mesh.vtk")
        (staging / "source.step").write_bytes(source)
        for name, value in (("recipe.json", recipe), ("intent.json", recipe["intent"])):
            (staging / name).write_text(json_text(value) + "\n", encoding="utf-8")
        report = result.report
        report["publication"] = {"source_sha256": hashlib.sha256(source).hexdigest(),
                                 "source_bytes": len(source), "staging_seconds": time.perf_counter() - started,
                                 "atomic_directory": True, "replaced_existing": False}
        (staging / "report.json").write_text(json_text(report) + "\n", encoding="utf-8")
        try:
            _rename_new_directory(staging, destination)
        except FileExistsError as error:
            raise WorkflowError("output_exists", "Another writer created the publication output",
                                stage="publication") from error
        except OSError as error:
            if error.errno in (errno.ENOTSUP, errno.ENOSYS, errno.EINVAL):
                raise WorkflowError("atomic_publication_unsupported", "The filesystem cannot publish this directory exclusively",
                                    {"detail": str(error)}, stage="publication") from error
            raise
        staging = None
        return {"path": str(destination), "seconds": time.perf_counter() - started,
                "source_sha256": result.source_sha256, "files": ["mesh.paacsmesh", "mesh.vtk", "source.step",
                                                                    "recipe.json", "intent.json", "report.json"]}
    except WorkflowError:
        raise
    except (SelectionError, IntentError) as error:
        raise WorkflowError(error.reason, str(error), error.report, stage="publication") from error
    except _native.ExportError as error:
        raise WorkflowError(error.reason, str(error), {"native_report": error.report,
                            "seconds": time.perf_counter() - started}, stage="publication") from error
    except (OSError, TypeError, ValueError, OverflowError) as error:
        raise WorkflowError("publication_failed", str(error), {"seconds": time.perf_counter() - started}, stage="publication") from error
    finally:
        if staging is not None and parent is not None and staging.parent == parent and staging.name.startswith("."):
            original_error = sys.exception()
            try:
                if staging.parent.resolve() != parent:
                    raise OSError("Refusing cleanup outside the verified publication parent")
                shutil.rmtree(staging)
            except OSError as cleanup_error:
                if isinstance(original_error, WorkflowError):
                    original_error.report["staging_cleanup_error"] = str(cleanup_error)
                else:
                    raise
