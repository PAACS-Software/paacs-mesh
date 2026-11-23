"""Analytic and independent high-precision references for Tet4 measurements."""

# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import json
import math
import sys
import unittest

from paacs_mesh.quality_metrics import METRIC_CONVENTIONS, tet4_metrics


RIGHT = ((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.))
REGULAR = ((0., 0., 0.), (1., 0., 0.), (.5, math.sqrt(3) / 2, 0.),
           (.5, math.sqrt(3) / 6, math.sqrt(2 / 3)))


def exact_determinant(points):
    """Leibniz determinant on exact supplied floats, independent of production."""
    points = [[Fraction(value) for value in point] for point in points]
    rows = [[points[i][j] - points[0][j] for j in range(3)] for i in (1, 2, 3)]
    result = Fraction(0)
    for permutation in itertools.permutations(range(3)):
        inversions = sum(permutation[i] > permutation[j] for i in range(3) for j in range(i + 1, 3))
        result += (-1) ** inversions * math.prod(rows[i][permutation[i]] for i in range(3))
    return result


def high_precision_reference(points):
    """Weighted-reference Jacobian and edge Gram geometry, not native helpers."""
    with localcontext() as context:
        context.prec = 800
        p = [[Decimal.from_float(float(value)) for value in point] for point in points]
        edges = [[p[i][j] - p[0][j] for j in range(3)] for i in (1, 2, 3)]
        determinant = exact_determinant(points)
        determinant = Decimal(determinant.numerator) / Decimal(determinant.denominator)
        weighted = (edges[0],
                    [(-edges[0][j] + 2 * edges[1][j]) / Decimal(3).sqrt() for j in range(3)],
                    [(-edges[0][j] - edges[1][j] + 3 * edges[2][j]) / Decimal(6).sqrt()
                     for j in range(3)])
        frobenius_squared = sum(value * value for vector in weighted for value in vector)
        mean_ratio = (3 * (Decimal(2).sqrt() * determinant) ** (Decimal(2) / 3) / frobenius_squared
                      if determinant > 0 else Decimal(0))
        squared = {(i, j): sum((p[i][axis] - p[j][axis]) ** 2 for axis in range(3))
                   for i, j in itertools.combinations(range(4), 2)}
        incident_products = [math.prod(squared[tuple(sorted((i, j)))] for j in range(4) if i != j)
                             for i in range(4)]
        jacobian = Decimal(2).sqrt() * determinant / max(incident_products).sqrt()
        angles = []
        for i, j in itertools.combinations(range(4), 2):
            k, l = [index for index in range(4) if index not in (i, j)]
            u, v, w = [[p[index][axis] - p[i][axis] for axis in range(3)] for index in (j, k, l)]
            dot = lambda a, b: sum(x * y for x, y in zip(a, b))
            # Projection onto a plane perpendicular to the common edge gives
            # the interior angle directly, without constructing face normals.
            cosine = dot(u, u) * dot(v, w) - dot(u, v) * dot(u, w)
            sine = dot(u, u).sqrt() * abs(determinant)
            scale = max(abs(cosine), sine)
            angles.append(math.degrees(math.atan2(float(sine / scale), float(cosine / scale))))
        return float(mean_ratio), float(jacobian), min(angles), determinant


class Tet4MetricTests(unittest.TestCase):
    def assertClose(self, actual, expected, relative=3e-13, absolute=0.0):
        self.assertTrue(math.isclose(actual, expected, rel_tol=relative, abs_tol=absolute),
                        f"{actual!r} != {expected!r}")

    def test_regular_tetrahedron_analytic_values_and_definitions(self):
        metrics = tet4_metrics(REGULAR)
        self.assertEqual(METRIC_CONVENTIONS["version"], "paacs-mesh-tet4-quality-v1")
        self.assertEqual(metrics["orientation"], "positive")
        self.assertClose(metrics["mean_ratio"], 1)
        self.assertClose(metrics["scaled_jacobian"], 1)
        self.assertClose(metrics["minimum_dihedral_deg"], math.degrees(math.acos(1 / 3)))
        self.assertClose(metrics["signed_volume"], 1 / (6 * math.sqrt(2)))
        for name in ("edge_min", "edge_median", "edge_max", "edge_rms", "volume_equivalent_edge"):
            self.assertClose(metrics[name], 1)
        self.assertFalse(metrics["range_limited"])
        json.dumps(metrics, allow_nan=False)

    def test_right_tetrahedron_analytic_values(self):
        metrics = tet4_metrics(RIGHT)
        self.assertClose(metrics["mean_ratio"], (4 / 3) * 2 ** (-2 / 3))
        self.assertClose(metrics["scaled_jacobian"], 1 / math.sqrt(2))
        self.assertClose(metrics["minimum_dihedral_deg"], math.degrees(math.acos(1 / math.sqrt(3))))
        self.assertClose(metrics["signed_volume"], 1 / 6)
        self.assertClose(metrics["edge_min"], 1)
        self.assertClose(metrics["edge_max"], math.sqrt(2))
        self.assertClose(metrics["edge_median"], (1 + math.sqrt(2)) / 2)
        self.assertClose(metrics["edge_rms"], math.sqrt(1.5))
        self.assertClose(metrics["volume_equivalent_edge"], 2 ** (1 / 6))

    def test_all_permutations_keep_geometry_and_record_orientation(self):
        original = tet4_metrics(RIGHT)
        for permutation in itertools.permutations(range(4)):
            points = [RIGHT[index] for index in permutation]
            metrics = tet4_metrics(points)
            sign = 1 if exact_determinant(points) > 0 else -1
            self.assertEqual(metrics["orientation"], "positive" if sign > 0 else "inverted")
            self.assertClose(metrics["scaled_jacobian"], sign / math.sqrt(2))
            self.assertClose(metrics["signed_volume"], sign / 6)
            self.assertClose(metrics["mean_ratio"], original["mean_ratio"] if sign > 0 else 0)
            self.assertClose(metrics["minimum_dihedral_deg"], original["minimum_dihedral_deg"])

    def test_reflection_changes_signed_values_without_changing_shape_angles(self):
        reflected = tet4_metrics([(-x, y, z) for x, y, z in RIGHT])
        self.assertEqual(reflected["orientation"], "inverted")
        self.assertEqual(reflected["mean_ratio"], 0)
        # This catches taking the minimum of four *negative* corner ratios:
        # that other convention would incorrectly produce -sqrt(2) here.
        self.assertClose(reflected["scaled_jacobian"], -1 / math.sqrt(2))
        self.assertClose(reflected["minimum_dihedral_deg"], tet4_metrics(RIGHT)["minimum_dihedral_deg"])

    def test_coplanar_repeated_and_coincident_nodes_have_no_dihedral(self):
        cases = (((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (1., 1., 0.)),
                 (RIGHT[0], RIGHT[1], RIGHT[2], RIGHT[2]), ((0., 0., 0.),) * 4)
        for points in cases:
            metrics = tet4_metrics(points)
            self.assertEqual(metrics["orientation"], "degenerate")
            for name in ("mean_ratio", "scaled_jacobian", "signed_volume", "volume_equivalent_edge"):
                self.assertEqual(metrics[name], 0)
            self.assertIsNone(metrics["minimum_dihedral_deg"])
            self.assertFalse(metrics["range_limited"])
            json.dumps(metrics, allow_nan=False)

    def test_uniform_scale_translation_and_rotation_preserve_shape(self):
        original = tet4_metrics(RIGHT)
        for exponent in (-400, -100, 0, 100, 400):
            # An exact axis rotation and power-of-two translation preserve the
            # represented geometry, including at very small and large scales.
            scale = math.ldexp(1.0, exponent)
            points = [(scale * (4 - y), scale * (8 + x), scale * (16 + z)) for x, y, z in RIGHT]
            metrics = tet4_metrics(points)
            for name in ("mean_ratio", "scaled_jacobian", "minimum_dihedral_deg"):
                self.assertClose(metrics[name], original[name])
            for name in ("edge_min", "edge_median", "edge_max", "edge_rms", "volume_equivalent_edge"):
                self.assertClose(metrics[name], original[name] * scale)
            if abs(exponent) == 400:
                self.assertIsNone(metrics["signed_volume"])
                self.assertTrue(metrics["range_limited"])
            else:
                self.assertClose(metrics["signed_volume"], math.ldexp(1 / 6, 3 * exponent))

    def test_independent_weighted_jacobian_and_gram_references(self):
        cases = (((.25, -.5, .125), (2., .25, -.5), (-.25, 1.5, .75), (.5, .25, 2.)),
                 ((0., 0., 0.), (1., 1., 1.), (1., 1. + 2 ** -40, 1.), (1., 1., 1. + 2 ** -40)),
                 ((0., 0., 0.), (1., 0., 0.), (.5, 1., 0.), (.3, .3, 2 ** -500)))
        for points in cases:
            metrics = tet4_metrics(points)
            mean_ratio, jacobian, angle, determinant = high_precision_reference(points)
            self.assertEqual(metrics["orientation"], "positive" if determinant > 0 else "inverted")
            self.assertClose(metrics["mean_ratio"], mean_ratio)
            self.assertClose(metrics["scaled_jacobian"], jacobian)
            self.assertClose(metrics["minimum_dihedral_deg"], angle, relative=2e-12)
            self.assertClose(metrics["signed_volume"], float(determinant / 6))

    def test_nonzero_extremely_thin_tetrahedron_is_not_flat_or_fake_ninety_degrees(self):
        points = (RIGHT[0], RIGHT[1], RIGHT[2], (0., 0., 2 ** -1000))
        metrics = tet4_metrics(points)
        reference = high_precision_reference(points)
        self.assertEqual(metrics["orientation"], "positive")
        self.assertGreater(metrics["mean_ratio"], 0)
        self.assertGreater(metrics["minimum_dihedral_deg"], 0)
        self.assertClose(metrics["mean_ratio"], reference[0])
        self.assertClose(metrics["scaled_jacobian"], reference[1])
        self.assertClose(metrics["minimum_dihedral_deg"], reference[2])
        self.assertFalse(metrics["range_limited"])

    def test_full_finite_coordinate_range_keeps_exact_orientation_and_json(self):
        tiny = math.ulp(0.)
        huge = sys.float_info.max
        points = ((-huge, 0., 0.), (huge, 0., 0.), (0., tiny, 0.), (0., 0., tiny))
        metrics = tet4_metrics(points)
        self.assertGreater(exact_determinant(points), 0)
        self.assertEqual(metrics["orientation"], "positive")
        self.assertIsNone(metrics["signed_volume"])
        self.assertIsNone(metrics["edge_max"])
        self.assertIsNotNone(metrics["edge_min"])
        self.assertIsNotNone(metrics["volume_equivalent_edge"])
        self.assertTrue(metrics["range_limited"])
        json.dumps(metrics, allow_nan=False)

    def test_only_some_overflowing_dimensional_metrics_are_null(self):
        scale = 1e308
        metrics = tet4_metrics(tuple(tuple(value * scale for value in point) for point in RIGHT))
        self.assertIsNone(metrics["signed_volume"])
        self.assertTrue(metrics["range_limited"])
        self.assertClose(metrics["edge_min"], scale)
        self.assertClose(metrics["edge_max"], math.sqrt(2) * scale)
        self.assertClose(metrics["volume_equivalent_edge"], 2 ** (1 / 6) * scale)
        self.assertClose(metrics["mean_ratio"], tet4_metrics(RIGHT)["mean_ratio"])

    def test_invalid_input_is_rejected_without_nonfinite_outputs(self):
        cases = (None, (), RIGHT[:3], (*RIGHT, RIGHT[0]), ((0., 0.),) * 4,
                 ((0., 0., 0., 0.),) * 4, ((math.nan, 0., 0.),) * 4,
                 ((math.inf, 0., 0.),) * 4, ((-math.inf, 0., 0.),) * 4,
                 ((10 ** 1000, 0., 0.),) * 4, ((object(), 0., 0.),) * 4)
        for points in cases:
            with self.assertRaises(ValueError):
                tet4_metrics(points)

    def test_input_iterators_are_consumed_with_fixed_coordinate_bounds(self):
        rows = []

        def endless_points():
            while True:
                rows.append(1)
                yield (0., 0., 0.)

        with self.assertRaises(ValueError):
            tet4_metrics(endless_points())
        self.assertEqual(len(rows), 5)
        coordinates = []

        def endless_coordinates():
            while True:
                coordinates.append(1)
                yield 0.

        with self.assertRaises(ValueError):
            tet4_metrics((endless_coordinates(), *RIGHT[1:]))
        self.assertEqual(len(coordinates), 4)


if __name__ == "__main__":
    unittest.main()
