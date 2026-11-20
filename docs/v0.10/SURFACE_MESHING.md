# V0.10 prepared surface and curve meshing

## Canonical edge discretization

`CanonicalCurveSet` is the new explicit shared-edge artifact.

- Every model edge is meshed once using `meshEdges()`.
- Records retain the natural endpoint-to-endpoint node order.
- Consumers request either natural or reversed orientation.
- Classifications and exact node identity are retained.
- The artifact carries a deterministic hash.
- Atomic append to another `mesh::Mesh` is available.

The box integration test produces twelve records and verifies that every
forward sequence equals the reverse of its reversed use.

## Prepared surface route

V0.10 retains the existing production surface stack:

- `EdgeMesher`;
- `FaceMesher` and scalable exact Delaunay2D;
- `ProductionSurfaceMesher`;
- `SurfaceMesherAdapter` for prepared loops/chains, holes and embedded curves;
- `SurfaceSizingPlanner`.

The prepared surface allocation probe completed with zero escaped failures:

```text
surface construction: 6236 injected allocation points, 0 escapes
surface publication:     15 injected allocation points, 0 escapes
```

## Geometry contract

`GEdge` now exposes default parameter bounds, periodicity and closest parameter.
`GFace` exposes default periodic/singularity queries and point projection.
These are minimal additions used by existing meshing algorithms; they are not
an attempted reproduction of the OCCT class hierarchy.

## Remaining boundary

The generic legacy `GFace::boundaryEdgeTags()` model still represents one
outer loop. Multiple loops and holes are supported through the richer prepared
`SurfaceMesherAdapter` request, not automatically inferred from the simple
built-in model interface. The OCCT adapter must translate its wires into that
prepared representation.
