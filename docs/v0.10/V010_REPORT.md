# V0.10 implementation report

## 1. Authoritative base

The source was derived from the complete V0.9.3 package. The V0.9.3 compact
Delaunay, exact predicates, deterministic parallel experiments, structured
engines and existing qualification tests are retained. The standalone patch
in this package applies directly to that source tree.

## 2. Integration strategy

The codebase already contained several increasingly capable production routes:

- T2-B boundary-conforming Tet4;
- T3-A local protected refinement;
- T4-A bounded Tet quality processing;
- T4-B Tet10 placement and whole-element Jacobian certification.

V0.10 adds `meshPreparedTetModel()`, which chooses the requested terminal stage
and runs the corresponding established route into a private `mesh::Mesh`.
After the route completes, V0.10 performs the final boundary/validity audit,
quality reporting, hashing and one non-throwing publication swap.

The original compatibility APIs remain present. V0.10 does not redirect old
legacy calls silently.

## 3. Public request levels

```cpp
PreOcctTetOutput::BoundaryConformingTet4
PreOcctTetOutput::RefinedTet4
PreOcctTetOutput::QualityTet4
PreOcctTetOutput::CertifiedTet10
```

`makePreOcctTetRequest()` constructs conservative resource budgets and wires
the nested T2-B/T3-A/T4-A/T4-B request hierarchy. `enableBoundedTetQuality()`
turns on the accepted local improvement path and computes the report capacity
required by its bounded outcome contract.

The normal quality route defaults to report-only T4-A plus bounded deterministic
interior smoothing. This allows a valid conforming mesh to publish without
pretending all remaining slivers have been eliminated. A caller can opt into
the stricter bounded-improvement outcome.

## 4. Geometry-independent adapters

`ModelBoundaryOwnershipProvider` derives face ownership from each solid's
explicit `boundaryFaceTags()`. Built-in/prepared models can use this directly;
an OCCT adapter can later supply a richer authoritative provider without
changing the Tet engine.

`ModelTet10PlacementProvider` handles:

- curve midsides through `GEdge::closestPoint()`;
- surface midsides through `GFace::paramAt()` and `pointAt()`.

`Kernel.hpp` now includes additive defaults for parameter bounds, periodicity,
closest parameter, surface periodicity/singularity and point projection.
Derived geometry backends are not forced to implement unused capabilities.

## 5. Canonical curves

`buildCanonicalCurves()` executes the existing production edge mesher once,
then records the natural node and Line2 order for every logical edge. A face can
request the same node sequence forward or reversed. This prepares the exact
shared-edge identity rule required by future OCCT faces.

## 6. Sizing

`MetricTensor3` provides a finite SPD representation for directional target
size. A fixed deterministic eigenvalue computation supports validation and a
`MetricSizeFieldAdapter` lowers the tensor conservatively to the minimum
principal size for existing scalar algorithms.

The existing compiled scalar spatial fields and mesh controls remain unchanged.
V0.10 does not claim general anisotropic 3D Delaunay refinement.

## 7. SDK and packaging

The API version is `0.10.0-final-pre-occt-candidate`.

- `MeshingSession::meshPreparedTet()` exposes the production route through the
  existing progress/cancellation/session model.
- The narrow C ABI adds built-in box and cylinder Tet4/Tet10 entrypoints while
  retaining the pre-release `paacs_v8_` session symbol prefix.
- The CLI can create, archive, reload and audit Tet4, Tet10 and Hex8 meshes.
- The installed package exports `PAACS::mesh` under `PAACSMeshV010` and now
  correctly imports its `Threads` dependency.

## 8. Scope of source changes

The patch contains approximately 2,300 added lines across 20 files. It adds
integration and geometry-provider layers rather than rewriting the accepted
recovery/refinement implementations.
