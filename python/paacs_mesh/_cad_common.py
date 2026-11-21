"""Validation and bounded JSON helpers for CAD selections and intent."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import json
import math
import numbers
import operator
import os
from pathlib import Path
import tempfile


MAX_JSON_BYTES = 8 * 1024 * 1024
TAG_SCHEME = "paacs-occt-tags-v1"
INTENT_SCHEMA = "paacs-cad-intent-v1"


class SelectionError(ValueError):
    """A selection could not be resolved or validated unambiguously."""

    def __init__(self, reason, message, report=None, *, stage="selection"):
        super().__init__(f"{stage}: {message}")
        self.stage, self.reason = stage, reason
        self.report = {} if report is None else json_copy(report)


class IntentError(ValueError):
    """A CAD intent manifest could not be saved or applied safely."""

    def __init__(self, reason, message, report=None):
        super().__init__(f"intent: {message}")
        self.stage, self.reason = "intent", reason
        self.report = {} if report is None else json_copy(report)


def json_text(value):
    try:
        text = json.dumps(value, ensure_ascii=False, allow_nan=False,
                          separators=(",", ":"), sort_keys=True)
    except RecursionError as error:
        raise ValueError("Metadata JSON nesting is too deep") from error
    if len(text.encode("utf-8")) > MAX_JSON_BYTES:
        raise ValueError("JSON exceeds the 8 MiB metadata budget")
    return text


def json_copy(value):
    return json.loads(json_text(value))


def integer(value, name, *, minimum=1):
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an integer, not bool")
    result = operator.index(value)
    if result < minimum or result > 2 ** 63 - 1:
        raise ValueError(f"{name} is outside its supported range")
    return result


def number(value, name, *, minimum=None):
    if isinstance(value, (bool, str, bytes)):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result) or (minimum is not None and result < minimum):
        requirement = "finite" if minimum is None else f"finite and at least {minimum}"
        raise ValueError(f"{name} must be {requirement}")
    return result


def boolean(value, name):
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be bool")
    return value


def text(value, name, *, maximum=256):
    if (not isinstance(value, str) or not value.strip() or "\0" in value or
            len(value) > maximum):
        raise ValueError(f"{name} must be a nonempty string of at most {maximum} characters")
    return value


def sequence(value, name, *, length=None, maximum=10000):
    if isinstance(value, (str, bytes, bytearray, dict)):
        raise ValueError(f"{name} must be a bounded sequence")
    try:
        count = len(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a bounded sequence") from error
    if count > maximum or (length is not None and count != length):
        raise ValueError(f"{name} has an unsupported length")
    try:
        return [value[index] for index in range(count)]
    except (TypeError, KeyError, IndexError) as error:
        raise ValueError(f"{name} must be an indexable sequence") from error


def point(value, name):
    return [number(item, name) for item in sequence(value, name, length=3)]


def bounds(value, name):
    lo, hi = sequence(value, name, length=2)
    lo, hi = point(lo, name), point(hi, name)
    if any(a > b for a, b in zip(lo, hi)):
        raise ValueError(f"{name} lower coordinates must not exceed upper coordinates")
    return [lo, hi]


def metric(value, name):
    if isinstance(value, numbers.Real) and not isinstance(value, bool):
        result = number(value, name, minimum=0)
        return result
    lo, hi = sequence(value, name, length=2)
    lo, hi = number(lo, name, minimum=0), number(hi, name, minimum=0)
    if lo > hi:
        raise ValueError(f"{name} lower bound must not exceed its upper bound")
    return [lo, hi]


def close(first, second, tolerance=1e-8):
    if isinstance(first, bool) or isinstance(second, bool):
        return first is second
    # IDs, counts and senses are exact even when their magnitude is large.
    if isinstance(first, int) or isinstance(second, int):
        return first == second
    if isinstance(first, (int, float)) and isinstance(second, (int, float)):
        return math.isclose(first, second, rel_tol=tolerance, abs_tol=tolerance)
    if isinstance(first, dict) and isinstance(second, dict):
        return first.keys() == second.keys() and all(close(first[key], second[key], tolerance)
                                                     for key in first)
    if isinstance(first, list) and isinstance(second, list):
        return len(first) == len(second) and all(close(a, b, tolerance) for a, b in zip(first, second))
    return first == second


def load_json(path):
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON field {key!r}")
            result[key] = value
        return result

    with Path(path).open("rb") as stream:
        payload = stream.read(MAX_JSON_BYTES + 1)
    if len(payload) > MAX_JSON_BYTES:
        raise ValueError("Intent exceeds the 8 MiB metadata budget")
    def reject_constant(value):
        raise ValueError(f"Nonfinite JSON number {value}")
    try:
        result = json.loads(payload.decode("utf-8"), object_pairs_hook=unique_pairs,
                            parse_constant=reject_constant)
    except RecursionError as error:
        raise ValueError("Intent JSON nesting is too deep") from error
    # json.loads can parse an overflowed exponent as inf without using parse_constant.
    json_text(result)
    return result


def protected_destination(path, source):
    destination = Path(path)
    if not str(destination) or "\0" in str(destination):
        raise ValueError("Intent output path must be nonempty and contain no NUL")
    original = Path(source)
    if destination.resolve() == original.resolve():
        raise ValueError("Intent output cannot overwrite the STEP source")
    if destination.exists() and original.exists() and os.path.samefile(destination, original):
        raise ValueError("Intent output cannot overwrite the STEP source")
    return destination


def atomic_json(path, value, source):
    destination = protected_destination(path, source)
    payload = (json_text(value) + "\n").encode("utf-8")
    if len(payload) > MAX_JSON_BYTES:
        raise ValueError("Intent exceeds the 8 MiB metadata budget")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="wb", dir=destination.parent,
                                         prefix=f".{destination.name}.", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return destination
