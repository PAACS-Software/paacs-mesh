"""Independent, scale-aware diagnostics for linear four-node tetrahedra.

These measurements describe an element, not its CAD fidelity or mesh validity.
The signed scaled Jacobian uses the common-denominator convention documented in
Sandia's Verdict implementation:
https://github.com/sandialabs/verdict/blob/master/V_TetMetric.cpp
Metric definitions are also described in:
https://coreform.com/cubit_help/latest/mesh_generation/mesh_quality_assessment/tetrahedral_metrics.htm
No third-party implementation is required or incorporated here.
"""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import math
from itertools import islice


METRIC_CONVENTIONS = {
    "version": "paacs-mesh-tet4-quality-v1",
    "orientation": "sign of dot(cross(p1-p0,p2-p0),p3-p0); zero is degenerate",
    "mean_ratio": "12*(3*V)^(2/3)/sum(six squared edges) for positive orientation; otherwise zero",
    "minimum_dihedral_deg": "minimum of six interior dihedral angles, in degrees; null for degeneracy",
    "scaled_jacobian": "sqrt(2)*signed determinant / maximum of four incident edge-length products",
    "signed_volume": "signed determinant / 6, in input length units cubed",
    "edge_median": "arithmetic mean of the third and fourth of six sorted edge lengths",
    "volume_equivalent_edge": "edge of a regular tetrahedron of equal absolute volume",
    "range_limited": "a nonzero measurement could not be represented as a finite nonzero binary64 value",
    "scope": "Tet4 shape only; no CAD, connectivity, overlap, or simulation certificate",
}

_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
_FACES = ((1, 3, 2), (0, 2, 3), (0, 3, 1), (0, 1, 2))
_INCIDENT = ((0, 1, 2), (0, 3, 4), (1, 3, 5), (2, 4, 5))
_ROOT_TWO = math.sqrt(2.0)


def _subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _physical(value, exponent, limited):
    """Rescale without confusing true zeros with numerical range loss."""
    try:
        result = math.ldexp(value, exponent)
    except OverflowError:
        result = math.inf
    if not math.isfinite(result) or (value != 0 and result == 0):
        limited[0] = True
        return None
    return result


def _integer_parts(value):
    """A positive integer as a bounded mantissa and a base-two exponent."""
    if not value:
        return 0.0, 0
    exponent = value.bit_length()
    # Integer true division rounds the bounded quotient without first converting
    # either possibly very large integer to float.
    mantissa = value / (1 << exponent)
    if mantissa == 1.0:
        return 0.5, exponent + 1
    return mantissa, exponent


def _sqrt_parts(value):
    mantissa, exponent = _integer_parts(value)
    quotient, remainder = divmod(exponent, 2)
    return math.sqrt(math.ldexp(mantissa, remainder)), quotient


def _cube_root_parts(mantissa, exponent):
    quotient, remainder = divmod(exponent, 3)
    return math.cbrt(math.ldexp(mantissa, remainder)), quotient


def _result(determinant, volume, mean_ratio, jacobian, angle, lengths,
            median, rms, equivalent, limited):
    return {
        "mean_ratio": mean_ratio,
        "minimum_dihedral_deg": angle,
        "scaled_jacobian": jacobian,
        "signed_volume": volume,
        "orientation": "positive" if determinant > 0 else "inverted" if determinant < 0 else "degenerate",
        "edge_min": lengths[0], "edge_median": median,
        "edge_max": lengths[-1], "edge_rms": rms,
        "volume_equivalent_edge": equivalent,
        "range_limited": limited[0],
    }


def _fast_metrics(points):
    edges = [_subtract(points[j], points[i]) for i, j in _PAIRS]
    largest = max(abs(value) for edge in edges for value in edge)
    if not math.isfinite(largest) or largest == 0:
        return None
    exponent = math.frexp(largest)[1]
    normalized = [tuple(math.ldexp(value, -exponent) for value in edge) for edge in edges]
    if any(value != 0 and scaled == 0
           for edge, scaled_edge in zip(edges, normalized)
           for value, scaled in zip(edge, scaled_edge)):
        return None
    squared = [math.fsum(value * value for value in edge) for edge in normalized]
    if any(value == 0 for value in squared):
        return None
    a, b, c = normalized[:3]
    terms = (a[0] * b[1] * c[2], a[1] * b[2] * c[0], a[2] * b[0] * c[1],
             -a[2] * b[1] * c[0], -a[1] * b[0] * c[2], -a[0] * b[2] * c[1])
    determinant = math.fsum(terms)
    permanent = math.fsum(abs(value) for value in terms)
    # The conservative filter includes subtraction and product roundoff. The
    # exact dyadic path below also handles subnormal products and cancellation.
    if abs(determinant) <= 1e-14 * permanent or abs(determinant) < 1e-280:
        return None
    lengths = [math.hypot(*edge) for edge in normalized]
    denominator = max(math.prod(lengths[i] for i in incident) for incident in _INCIDENT)
    if denominator == 0:
        return None
    normalized_points = ((0.0, 0.0, 0.0), a, b, c)
    normals = []
    for i, j, k in _FACES:
        normal = _cross(_subtract(normalized_points[j], normalized_points[i]),
                        _subtract(normalized_points[k], normalized_points[i]))
        size = math.hypot(*normal)
        if size < 1e-140:
            return None
        normals.append(tuple(value / size for value in normal))
    angles = [math.degrees(math.atan2(math.hypot(*_cross(normals[i], normals[j])),
                                     -math.fsum(x * y for x, y in zip(normals[i], normals[j]))))
              for i, j in _PAIRS]
    # An underflowed cross product must not manufacture a zero angle.
    if min(angles) == 0:
        return None
    limited = [False]
    ordered = sorted(lengths)
    physical_lengths = [_physical(value, exponent, limited) for value in ordered]
    equivalent, equivalent_exponent = _cube_root_parts(abs(determinant) * _ROOT_TWO, 3 * exponent)
    mean_ratio = (12 * (determinant / 2) ** (2 / 3) / math.fsum(squared)
                  if determinant > 0 else 0.0)
    return _result(determinant, _physical(determinant / 6, 3 * exponent, limited),
                   min(1.0, mean_ratio), max(-1.0, min(1.0, _ROOT_TWO * determinant / denominator)),
                   min(angles), physical_lengths,
                   _physical((ordered[2] + ordered[3]) / 2, exponent, limited),
                   _physical(math.sqrt(math.fsum(squared) / 6), exponent, limited),
                   _physical(equivalent, equivalent_exponent, limited), limited)


def _exact_points(points):
    """Represent finite binary64 coordinates as integers times one power of two."""
    parts = []
    for point in points:
        for value in point:
            numerator, denominator = value.as_integer_ratio()
            if numerator:
                trailing = (abs(numerator) & -abs(numerator)).bit_length() - 1
                parts.append((numerator >> trailing, trailing - (denominator.bit_length() - 1)))
            else:
                parts.append((0, 0))
    exponent = min((power for value, power in parts if value), default=0)
    integers = [value << (power - exponent) if value else 0 for value, power in parts]
    return [tuple(integers[i:i + 3]) for i in range(0, 12, 3)], exponent


def _exact_metrics(points):
    points, exponent = _exact_points(points)
    edges = [_subtract(points[j], points[i]) for i, j in _PAIRS]
    squared = [_dot(edge, edge) for edge in edges]
    determinant = _dot(_cross(edges[0], edges[1]), edges[2])
    determinant_mantissa, determinant_exponent = _integer_parts(abs(determinant))
    total_mantissa, total_exponent = _integer_parts(sum(squared))
    limited = [False]
    ordered = sorted(squared)
    length_parts = [_sqrt_parts(value) for value in ordered]
    physical_lengths = [_physical(value, power + exponent, limited) for value, power in length_parts]
    # Align before adding: either central length may individually exceed the
    # binary64 range although their mean divided by two is representable.
    median_exponent = max(length_parts[2][1], length_parts[3][1])
    median = sum(math.ldexp(value, power - median_exponent)
                 for value, power in length_parts[2:4]) / 2
    rms, rms_exponent = _sqrt_parts(sum(squared))
    equivalent, equivalent_exponent = _cube_root_parts(
        determinant_mantissa * _ROOT_TWO, determinant_exponent + 3 * exponent)
    orientation = (determinant > 0) - (determinant < 0)
    mean_ratio = 0.0
    jacobian = 0.0
    angle = None
    if determinant:
        if orientation > 0:
            quotient, remainder = divmod(2 * determinant_exponent, 3)
            value = 12 * (determinant_mantissa / 2) ** (2 / 3) / total_mantissa
            value *= 2 ** (remainder / 3)
            mean_ratio = _physical(value, quotient - total_exponent, limited)
            mean_ratio = 0.0 if mean_ratio is None else min(1.0, mean_ratio)
        product = max(math.prod(squared[i] for i in incident) for incident in _INCIDENT)
        denominator, denominator_exponent = _sqrt_parts(product)
        jacobian = _physical(orientation * _ROOT_TWO * determinant_mantissa / denominator,
                             determinant_exponent - denominator_exponent, limited)
        jacobian = 0.0 if jacobian is None else max(-1.0, min(1.0, jacobian))
        normals = [_cross(_subtract(points[j], points[i]), _subtract(points[k], points[i]))
                   for i, j, k in _FACES]
        angles = []
        for i, j in _PAIRS:
            cross = _cross(normals[i], normals[j])
            dot = -_dot(normals[i], normals[j])
            scale = max(abs(value).bit_length() for value in (*cross, dot))
            divisor = 1 << scale
            sine = math.hypot(*(value / divisor for value in cross))
            if sine == 0 and any(cross):
                limited[0] = True
            angles.append(math.degrees(math.atan2(sine, dot / divisor)))
        angle = min(angles)
    return _result(determinant, _physical(orientation * determinant_mantissa / 6,
                                        determinant_exponent + 3 * exponent, limited),
                   mean_ratio, jacobian, angle, physical_lengths,
                   _physical(median, median_exponent + exponent, limited),
                   _physical(rms / math.sqrt(6), rms_exponent + exponent, limited),
                   _physical(equivalent, equivalent_exponent, limited), limited)


def tet4_metrics(points):
    """Measure four finite XYZ coordinates without NumPy or engine state.

    Node order controls orientation and signed volume/Jacobian. Mean ratio is
    zero for nonpositive orientation, matching the native mixed-cell report.
    Inverted cells still have geometric dihedral and edge measurements.
    Undefined angles and unrepresentable dimensional values are ``None``.
    ``range_limited`` distinguishes numerical range loss from a true zero.
    Coordinates are interpreted as their supplied binary64 values; sub-ULP
    features lost before this call cannot be recovered. Most cells use bounded
    floating arithmetic. Uncertain orientation or extreme ranges use exact
    dyadic integer arithmetic with bounded-size finite binary64 inputs.
    """
    try:
        rows = tuple(islice(iter(points), 5))
        if len(rows) != 4:
            raise ValueError("Tet4 metrics require four finite XYZ coordinates")
        values = tuple(tuple(float(value) for value in islice(iter(point), 4)) for point in rows)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Tet4 metrics require four finite XYZ coordinates") from exc
    if (len(values) != 4 or any(len(point) != 3 for point in values)
            or not all(math.isfinite(value) for point in values for value in point)):
        raise ValueError("Tet4 metrics require four finite XYZ coordinates")
    return _fast_metrics(values) or _exact_metrics(values)
