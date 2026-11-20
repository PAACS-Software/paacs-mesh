# OCCT readiness review

## 1. Exact Tet engine

**Ready as a geometry-independent foundation.** The compact exact Delaunay
constructor remains private, deterministic and public-64-bit compatible. The
constrained pipeline is exposed through a single V0.10 request.

## 2. Constraints

**Ready for the supported prepared contract.** Existing T2-B infrastructure
recovers explicit edges and interacting Tri3 facets, preserves classifications
and proves exact boundary bijection. Arbitrary dirty-CAD cases require an OCCT
corpus and may identify additional local recovery templates.

## 3. Refinement

**Architecturally ready and bounded.** T3-A accepts the normal size field,
updates locally and protects the immutable recovered boundary. Curvature,
proximity and thickness data can enter through scalar or compiled fields
without changing the refinement transaction model.

## 4. Surface meshing

**Ready through the prepared-surface interface.** The scalable 2D constructor,
constraint recovery and multi-loop/hole adapter already exist. OCCT must
translate wires, seams and classifications into the prepared request.

## 5. Canonical edges

**Ready.** V0.10 explicitly builds each logical edge once and supports reversed
reuse by adjacent faces.

## 6. Tet10

**Ready for provider integration.** T4-B consumes an abstract placement
provider and independently performs identity checks, deterministic elevation,
Jacobian certification and atomic publication.

## 7. Geometry providers

**Ready at the required minimal level.** Curve and surface evaluation,
projection, derivatives/normals, curvature, periodicity and parameter bounds
are represented without OCCT dependencies.

## 8. Topology

**Ready for translation.** Persistent entity tags, surface chain identities,
explicit solid boundary faces, region ownership and portal identity already
exist. OCCT must populate them from its `TShape + Location` topology.

## 9. Hybrid bodies

**Ready for prepared geometry.** Tet4/Tet10, mapped Hex8, prepared prism/Hex
sweeps, transitions and polyhedral bodies can be combined through the existing
plan/executor.

## 10. Expected next architecture

Adding OCCT should now principally mean:

```text
STEP/BREP import and healing
        ↓
persistent PAACS topology extraction
        ↓
GEdge/GFace evaluation and projection
        ↓
authoritative ownership/orientation providers
        ↓
existing V0.10 meshing core
```

The expected work is real-CAD adaptation and qualification, not replacement of
the core meshing algorithms.
