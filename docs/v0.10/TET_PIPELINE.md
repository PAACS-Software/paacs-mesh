# V0.10 Tet production pipeline

## Normal execution

```text
PAACS Model + MeshSpec
        ↓
ModelBoundaryOwnershipProvider
        ↓
ProductionSurfaceMesher / SurfaceMesherAdapter
        ↓
BoundaryConformingMesher (T2-B)
        ↓
optional TetRefinementMesher (T3-A)
        ↓
optional TetQualityMesher (T4-A)
        ↓
optional Tet10CurvingMesher (T4-B)
        ↓
V0.10 final audit / quality / hash
        ↓
atomic publication to mesh::Mesh
```

## Delaunay foundation

The V0.9.3 `SerialExact` compact Delaunay path remains the default initial Tet
constructor. Public indices remain 64-bit. The private compact arena, exact
filtered predicates, canonical ordering, reusable workspace and safe legacy
fallback are unchanged.

The slower V0.9.3 deterministic-parallel policies remain explicit options and
are not selected by the V0.10 production orchestrator.

## Constraints and identity

The T2-B route retains:

- explicit artifact edge/facet identities;
- exact edge and interacting facet recovery;
- authoritative region ownership and interface sense;
- indexed local scheduling;
- topology-based represented-domain extraction;
- exact boundary bijection;
- atomic whole-model publication.

V0.10 does not substitute coordinate proximity for constraint identity.

## Refinement

T3-A remains a bounded, local, deterministic refinement stage over an immutable
recovered boundary. Its scheduler keeps a canonical violation queue and updates
only topology affected by committed deltas. The retained probe demonstrated:

```text
committed insertions: 1
blocked rows:         13
work units:           11,749
per-insertion scans:  0
context rebuilds:     0
```

The example intentionally ends as `Incomplete`, proving that the boundary and
private authoritative result survive when the requested refinement objective
cannot be completed under the budget.

## Quality

V0.10 exposes two levels:

1. report-only accepted T4-A plus deterministic interior smart-Laplacian
   smoothing and final exact audit;
2. opt-in bounded T4-A cavity improvement through
   `enableBoundedTetQuality()`.

The optimizer is deliberately not described as universal sliver removal.
Protected edges/facets and interfaces remain immutable.

## Representation boundary

The initial Delaunay constructor uses the compact private arena. Existing
qualified recovery/refinement stages continue to use their authoritative
transactional topology workspaces. V0.10 removes public-route fragmentation,
but it does not complete a zero-copy compact-arena-to-recovery unification.
That deeper internal consolidation can be revisited only if profiling after
real CAD integration justifies the qualification risk.
