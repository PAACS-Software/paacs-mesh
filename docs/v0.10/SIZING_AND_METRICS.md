# V0.10 sizing and metric foundation

## Existing scalar controls retained

The existing `SizeField`, compiled spatial fields, face/body controls, growth
limits and structured-engine planners are unchanged. They continue to feed
edge, surface, Tet, Hex and layer paths through `MeshSpec` and compiled control
objects.

## SPD metric foundation

V0.10 adds:

```cpp
MetricTensor3
MetricField
ConstantMetricField
AdditiveMetricField
MetricSizeFieldAdapter
DirectionalMetricControl
```

For a displacement `dx`, metric length is:

```text
sqrt(dx^T M dx)
```

Validation includes:

- finite tensor entries;
- symmetric positive definiteness;
- deterministic eigenvalue bounds;
- scalar size envelope checks.

`directionalMetric()` creates a metric with one desired size along a supplied
axis and a transverse size in the perpendicular plane.

The conservative scalar adapter returns:

```text
1 / sqrt(maximum metric eigenvalue)
```

so current scalar meshers can consume a directional control safely without
claiming metric-Delaunay behavior.

## Scope

V0.10 provides a stable representation and composition boundary for future
curvature, proximity, thickness and solution fields. It does not implement a
complete anisotropic 3D Delaunay refinement theory.
