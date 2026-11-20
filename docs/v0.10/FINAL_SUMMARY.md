# PAACS Mesh V0.10 — Final Pre-OCCT Meshing Core

## Status

V0.10 is a complete-source **engineering candidate** based directly on the
attached V0.9.3 deterministic-parallel Delaunay package. It preserves V0.9.3's
compact exact `SerialExact` Delaunay implementation as the normal production
backend and does not alter its accepted canonical point-cloud results.

The principal V0.10 change is integration. The previously separate accepted
boundary-conforming Tet4, refinement, bounded quality, Tet10, surface, layer,
transition and hybrid components are now exposed through a coherent
geometry-independent workflow that an OCCT adapter can feed without exposing
OCCT types to the meshers.

## New production-facing pieces

- `PreOcctMeshingCore`: one atomic Tet request selects boundary-conforming
  Tet4, refined Tet4, quality Tet4 or Jacobian-certified Tet10.
- `ModelBoundaryOwnershipProvider`: translates the abstract PAACS Model's
  explicit solid/face topology into authoritative volume ownership records.
- `ModelTet10PlacementProvider`: supplies geometry-independent curve and
  surface placement through `GEdge` and `GFace` capabilities.
- `CanonicalCurveMesher`: meshes every topological edge once and provides
  forward/reversed reuse without coordinate-based resampling.
- `MetricTensor3` and `MetricField`: a deterministic SPD sizing foundation,
  with a conservative scalar adapter for existing meshers.
- V0.10 `MeshingSession`, C API and CLI routes for classified Tet4/Tet10
  generation from prepared models.
- Installable CMake package `PAACSMeshV010` with a verified external consumer.

## Retained engines

The existing generation engines remain distinct:

- exact constrained-Delaunay Tet;
- mapped Hex8;
- prepared mixed Tri3/Quad4 sweep to Prism6/Hex8;
- Pyramid5/Prism6/Hex8/Tet4 transition catalogue;
- structured polyhedral generation;
- identity-based hybrid assembly.

Tet output can be supplied to `HybridMeshingExecutor` through its existing
`PreparedFixedMesh` input. No generic monolithic element generator was added.

## Representative completed results

| Workflow | Result |
|---|---|
| Box quality Tet4 | 26 nodes, 43 volume Tet4, zero inverted, exact boundary audit, hash `16723556976053534829` |
| Box certified Tet10 | 118 nodes, 43 volume Tet10, zero invalid Tet10, hash `16508864188371921865` |
| Cylinder quality Tet4 | 182 nodes, 428 volume Tet4, zero inverted, exact boundary audit, hash `8085090672163299017` |
| Canonical curves | 12 box edges built once; forward/reverse reuse verified |
| Prepared layers | Tri3 source generated two positive Prism6 layers |
| Transition | Prism6 converted to positive Tet4 through retained catalogue |
| Bounded refinement | One local insertion committed before a structured incomplete outcome; no global scans/context rebuilds |
| Delaunay 200k | 1,346,539 Tet4, canonical hash `0xfcd63f974312b432`, 61.6 MB compact arena |

## Verification completed

- Clang 17 Release build with warnings as errors.
- Focused CTest: 4/4 passed.
- Full production library built under ASan+UBSan; focused V0.10 and Delaunay
  sanitizer tests: 2/2 passed with LeakSan enabled.
- Existing prepared-surface allocation probe passed.
- Existing volume-engine probe passed.
- Existing T3-A refinement probe demonstrated one committed insertion and
  structured rollback behavior.
- Delaunay goldens and canonical 50k/200k hashes preserved.
- Tet4 output identical for one and four requested workers on the tested box.
- CLI generation/readback/audit for Tet4, Tet10, cylinder Tet4 and Hex8.
- C API consumer generated a classified Tet4 box.
- Installed CMake package consumed successfully by a separate project.
- Clean V0.9.3 patch application and packaged-source rebuild are recorded in
  `PACKAGE_VERIFICATION.md`.

## Honest release boundary

This is not a claim that arbitrary CAD is already solved. The following remain
for the OCCT/hardening phase:

- OCCT/STEP/BREP import, healing and authoritative orientation extraction;
- arbitrary real-CAD surface corpus qualification;
- universal sliver removal;
- general anisotropic 3D Delaunay refinement;
- complete curved-layer corner/junction classification;
- full cross-platform, TSan and exhaustive allocation-failure qualification.

The important readiness result is architectural: the next phase can implement
OCCT topology and geometry providers around the existing PAACS algorithms
rather than replacing the meshing core.
