"""Atomic, source-bound companion groups for CAD selection workflows."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
from pathlib import Path
from uuid import uuid4

from ._cad_common import (
    INTENT_SCHEMA, TAG_SCHEME, IntentError, SelectionError, atomic_json, boolean,
    bounds, close, integer, json_copy, json_text, load_json, number, point, sequence, text,
)
from .cad import _expectation, _query_values, _tags


def _sha256(value):
    value = text(value, "source SHA-256", maximum=64)
    if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
        raise ValueError("Source SHA-256 must contain 64 lowercase hexadecimal characters")
    return value


def _guard(value):
    required = {"tag", "surface", "owners", "owner_uses"}
    allowed = required | {"area", "tolerance", "area_centroid", "bounds", "radius", "normal"}
    if not isinstance(value, dict) or not required <= value.keys() or value.keys() - allowed:
        raise ValueError("Invalid face property guard fields")
    result = {"tag": integer(value["tag"], "guard tag"), "surface": text(value["surface"], "guard surface"),
              "owners": list(_tags(sequence(value["owners"], "guard owners"), "guard owners"))}
    uses = []
    for use in sequence(value["owner_uses"], "owner uses"):
        if not isinstance(use, dict) or set(use) != {"body_tag", "sense", "outward_normal"}:
            raise ValueError("Invalid owner use guard")
        sense = integer(use["sense"], "owner sense", minimum=-1)
        if sense not in (-1, 0, 1):
            raise ValueError("Owner sense must be -1, 0 or 1")
        outward = use["outward_normal"]
        uses.append({"body_tag": integer(use["body_tag"], "owner body"), "sense": sense,
                     "outward_normal": None if outward is None else point(outward, "outward_normal")})
    uses.sort(key=lambda use: (use["body_tag"], use["sense"], json_text(use["outward_normal"])))
    if sorted({use["body_tag"] for use in uses}) != result["owners"]:
        raise ValueError("Owner use guards contradict the owner tags")
    result["owner_uses"] = uses
    if "normal" in value:
        result["normal"] = point(value["normal"], "normal")
    for key in ("area", "tolerance", "radius"):
        if key in value:
            result[key] = number(value[key], key, minimum=0)
    if "area_centroid" in value:
        result["area_centroid"] = point(value["area_centroid"], "area_centroid")
    if "bounds" in value:
        box = value["bounds"]
        if not isinstance(box, dict) or set(box) != {"lower", "upper"}:
            raise ValueError("Invalid guard bounds")
        lo, hi = bounds([box["lower"], box["upper"]], "guard bounds")
        result["bounds"] = {"lower": lo, "upper": hi}
    return result


def _group(value):
    fields = {"name", "id", "tags", "body", "query", "guards"}
    if not isinstance(value, dict) or set(value) != fields:
        raise ValueError("Invalid companion group fields")
    name = text(value["name"], "group name")
    identifier = text(value["id"], "group ID", maximum=128)
    tags = list(_tags(sequence(value["tags"], "group tags"), "group tags"))
    body = None if value["body"] is None else integer(value["body"], "group body")
    query = _query_values(value["query"])
    if "expect" not in query:
        raise ValueError("Saved selection query must include an expected count")
    lo, hi = query["expect"]
    if not lo <= len(tags) <= hi:
        raise ValueError("Group size contradicts its saved query")
    if query.get("body") != body:
        raise ValueError("Group body contradicts its saved query")
    guards = [_guard(item) for item in sequence(value["guards"], "group guards")]
    guards.sort(key=lambda item: item["tag"])
    if [item["tag"] for item in guards] != tags:
        raise ValueError("Group guards must describe every selected tag exactly once")
    return name, {"id": identifier, "tags": tags, "body": body, "query": query, "guards": guards}


def define_group(model, name, selection, group_id):
    try:
        name = text(name, "group name")
        model._require_selection(selection)
        selection.validate()
        previous = json.loads(model._group_records[name]) if name in model._group_records else None
        identifier = previous["id"] if group_id is None and previous else str(uuid4()) if group_id is None else group_id
        identifier = text(identifier, "group ID", maximum=128)
        for other_name, stored in model._group_records.items():
            if other_name != name and json.loads(stored)["id"] == identifier:
                raise SelectionError("duplicate_group_id", "A different companion group already has this ID",
                                     {"id": identifier, "group": other_name})
        record = {"id": identifier, "tags": list(selection.tags), "body": selection.body,
                  "query": selection.query, "guards": selection.guards}
        candidate = dict(model._group_records)
        candidate[name] = json_text(record)
        json_text({key: json.loads(value) for key, value in candidate.items()})
        model._group_records = candidate
        return json_copy({"name": name, **record})
    except SelectionError:
        raise
    except (TypeError, ValueError, OverflowError) as error:
        raise SelectionError("invalid_group", str(error)) from error


def save_intent(model, path, revision):
    try:
        model._assert_identity()
        revision = model._intent_revision if revision is None else text(revision, "revision")
        groups = []
        for name, stored in sorted(model._group_records.items()):
            record = json.loads(stored)
            model._validate_group_record(record)
            groups.append({"name": name, **record})
        info = model.info
        manifest = {"schema": INTENT_SCHEMA,
                    "source": {"sha256": model.source_sha256, "bytes": info["file_bytes"],
                               "name": Path(info["source_path"]).name, "units": "mm"},
                    "kernel_version": model.kernel_version, "tag_scheme": TAG_SCHEME,
                    "revision": revision, "groups": groups}
        destination = atomic_json(path, manifest, info["source_path"])
        model._intent_revision = revision
        return {"path": str(destination), "schema": INTENT_SCHEMA, "groups": len(groups),
                "source_sha256": model.source_sha256, "revision": revision}
    except IntentError:
        raise
    except SelectionError as error:
        raise IntentError(error.reason, "A companion group failed validation", error.report) from error
    except (OSError, TypeError, ValueError, OverflowError) as error:
        raise IntentError("save_failed", str(error)) from error


def load_intent(model, path, rebind):
    try:
        rebind = boolean(rebind, "rebind")
        model._assert_identity()
        manifest = load_json(path)
        required = {"schema", "source", "kernel_version", "tag_scheme", "groups"}
        if not isinstance(manifest, dict) or not required <= manifest.keys() or manifest.keys() - (required | {"revision"}):
            raise ValueError("Invalid intent manifest fields")
        if manifest["schema"] != INTENT_SCHEMA:
            raise IntentError("unknown_schema", "Unsupported CAD intent schema", {"schema": manifest["schema"]})
        if manifest["tag_scheme"] != TAG_SCHEME:
            raise IntentError("tag_scheme_mismatch", "Unsupported CAD tag scheme")
        source = manifest["source"]
        if not isinstance(source, dict) or set(source) != {"sha256", "bytes", "name", "units"}:
            raise ValueError("Invalid intent source metadata")
        source_hash = _sha256(source["sha256"])
        source_bytes = integer(source["bytes"], "source bytes", minimum=0)
        text(source["name"], "source name", maximum=1024)
        if source["units"] != "mm":
            raise IntentError("units_mismatch", "CAD intent must use millimetres")
        kernel = text(manifest["kernel_version"], "kernel version")
        revision = manifest.get("revision")
        if revision is not None:
            revision = text(revision, "revision")
        source_matches = source_hash == model.source_sha256 and source_bytes == model.info["file_bytes"]
        identity_matches = source_matches and kernel == model.kernel_version
        if not identity_matches and not rebind:
            reason = "source_mismatch" if not source_matches else "kernel_mismatch"
            raise IntentError(reason, "Intent does not match this import; remapping requires explicit rebind=True",
                              {"saved_source_sha256": source_hash, "source_sha256": model.source_sha256,
                               "saved_kernel_version": kernel, "kernel_version": model.kernel_version})
        candidate, identifiers, changes = {}, set(), []
        for value in sequence(manifest["groups"], "groups"):
            name, record = _group(value)
            if name in candidate or record["id"] in identifiers:
                raise IntentError("duplicate_group", "Intent contains duplicate group names or IDs")
            identifiers.add(record["id"])
            if identity_matches:
                if record["body"] is not None and record["body"] not in model.body_tags:
                    raise IntentError("unknown_body", "Intent references an unknown body")
                model._validate_group_record(record)
            else:
                query = dict(record["query"])
                forbidden = {"tags", "body", "adjacent_to", "group"}.intersection(query)
                if forbidden or not set(query) - {"expect", "tolerance"}:
                    raise IntentError("unsafe_rebind", "This group has no independent semantic/geometric query",
                                      {"group": name, "model_local_fields": sorted(forbidden)})
                # Count is fixed to the saved selection, even if its original query used a range.
                count = len(record["tags"])
                query["expect"] = [count, count]
                selection = model.select_faces(**query)
                before = sorted(guard["surface"] for guard in record["guards"])
                after = sorted(guard["surface"] for guard in selection.guards)
                if before != after:
                    raise IntentError("surface_kind_changed", "Rebinding changed the selected surface kinds", {"group": name})
                changes.append({"group": name, "id": record["id"], "old_tags": record["tags"],
                                "new_tags": list(selection.tags), "old_body": record["body"], "new_body": selection.body})
                record = {"id": record["id"], "tags": list(selection.tags), "body": selection.body,
                          "query": selection.query, "guards": selection.guards}
            candidate[name] = json_text(record)
        json_text({key: json.loads(value) for key, value in candidate.items()})
        # Every record has passed before replacing any existing group or revision.
        model._group_records, model._intent_revision = candidate, revision
        return {"schema": INTENT_SCHEMA, "source_matched": source_matches,
                "rebound": not identity_matches, "groups": len(candidate), "changes": changes,
                "source_sha256": model.source_sha256, "kernel_version": model.kernel_version,
                "revision": revision}
    except IntentError:
        raise
    except SelectionError as error:
        raise IntentError(error.reason, "An intent selection failed validation", error.report) from error
    except (OSError, UnicodeError, TypeError, ValueError, OverflowError) as error:
        raise IntentError("invalid_manifest", str(error)) from error
