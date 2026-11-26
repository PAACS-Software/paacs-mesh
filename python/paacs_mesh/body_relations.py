"""Bounded, opt-in geometric relations between explicitly selected CAD bodies."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from ._cad_common import SelectionError, integer, json_copy, number, sequence


def inspect_body_relations(model, pairs, *, distance_tolerance=1e-7, max_work=100000,
                           max_working_bytes=64 * 1024 * 1024, max_seconds=30):
    """Query selected body pairs without constructing or changing a mesh.

    A contact classification describes native geometry, not a conforming mesh
    interface or a choice of bonded/contact physics. Near contacts that cannot
    be established by the native intersection remain unresolved. The caller
    supplies at most 64 pairs per query; no quadratic all-body scan is implicit.
    """
    try:
        try:
            model._assert_identity()
        except SelectionError as error:
            raise SelectionError(error.reason, "The imported model identity changed", error.report,
                                 stage="body_relations") from error
        requested = sequence(pairs, "pairs", maximum=64)
        if not requested:
            raise ValueError("pairs must contain at least one pair")
        known = set(model.body_tags)
        selected, seen = [], set()
        for item in requested:
            first, second = sequence(item, "body pair", length=2)
            first, second = integer(first, "first body tag"), integer(second, "second body tag")
            if first == second:
                raise ValueError("A body pair must contain two different body tags")
            if first not in known or second not in known:
                raise SelectionError("unknown_body", "A requested body is absent from the model",
                                     {"body_tags": [first, second]}, stage="body_relations")
            key = tuple(sorted((first, second)))
            if key in seen:
                raise ValueError("pairs contains a duplicate unordered body pair")
            seen.add(key)
            selected.append(list(key))
        tolerance = number(distance_tolerance, "distance_tolerance", minimum=0)
        work = integer(max_work, "max_work")
        memory = integer(max_working_bytes, "max_working_bytes")
        seconds = number(max_seconds, "max_seconds", minimum=0)
        if seconds == 0:
            raise ValueError("max_seconds must be positive")
        # Keep the native conversion and returned diagnostic metadata bounded.
        if work > 10000000 or memory > 1024 * 1024 * 1024:
            raise ValueError("Body relation caps are 10000000 work units and 1 GiB tracked storage")
        result = model._native._body_relations(selected, distance_tolerance=tolerance,
            max_work=work, max_working_bytes=memory, max_seconds=seconds)
        result = json_copy(result)
        if (result.get("source_sha256") != model.source_sha256 or
                result.get("model_identity") != model.model_identity):
            raise SelectionError("model_changed", "Body relations do not match the imported model",
                                 stage="body_relations")
        if not result.get("completed", False):
            raise SelectionError(result.get("reason", "query_failed"),
                                 "Body relation query did not complete", result, stage="body_relations")
        return result
    except SelectionError:
        raise
    except (TypeError, ValueError, OverflowError) as error:
        raise SelectionError("invalid_argument", str(error), stage="body_relations") from error
    except (RuntimeError, OSError) as error:
        raise SelectionError("query_failed", str(error), stage="body_relations") from error
