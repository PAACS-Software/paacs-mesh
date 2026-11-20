# Geometry-provider boundary for OCCT

## Purpose

The meshers consume the existing PAACS geometry model and small specialized
providers. OCCT should translate into these capabilities rather than leaking
OCCT classes into algorithm code.

## Curve capabilities

Current `GEdge` capabilities include:

- persistent tag/identity;
- start/end vertex tags;
- parameter bounds;
- point and tangent evaluation;
- length;
- curvature;
- periodicity;
- closest parameter and closest point.

## Surface capabilities

Current `GFace` capabilities include:

- persistent tag/identity;
- parameter bounds;
- point and normal evaluation;
- inverse parameter mapping / projection;
- curvature;
- metric scale;
- trimmed-domain containment;
- U/V periodicity;
- singularity query;
- ordered boundary edge tags for the simple model interface.

Prepared multi-wire surface requests carry richer loop/chain/orientation data.

## Topology and ownership

`ModelBoundaryOwnershipProvider` translates:

```text
solid.boundaryFaceTags()
        ↓
face persistent identity
        ↓
region owner + directed sense
```

The built-in implementation determines sense from face normal and solid
inside/outside evaluation. A real OCCT adapter should supply authoritative
oriented face uses directly when available.

## Tet10 placement

`ModelTet10PlacementProvider` maps classified source edges to the corresponding
abstract curve/surface and returns a projected midside coordinate with a
provider stamp. The existing T4-B engine handles identity, budgets,
certification and refusal.

## Intended OCCT adapter

```text
OCCT import/healing
       ↓
persistent PAACS vertices/edges/wires/faces/shells/bodies
       ↓
GEdge/GFace evaluation + projection implementations
       ↓
authoritative ownership/orientation provider
       ↓
existing V0.10 curve/surface/Tet/layer/hybrid algorithms
```
