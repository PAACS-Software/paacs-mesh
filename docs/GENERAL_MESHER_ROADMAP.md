# Next milestone: dependable general tetrahedral meshing

PAACS now has useful specialized routes and reusable controls. These are a good
foundation for a geometry router: use a verified mapped or swept mesh where it
fits, and a constrained tetrahedral route elsewhere. Modern language/toolchain
choices support implementation, but geometry coverage, quality, resource use,
and dependable failure behavior determine whether this is a general mesher.

## What the new probes actually establish

The complex sweep corpus exercises an I-beam, a block with six polygonal bores,
a star duct, twelve cooling channels, a notched frame with two holes, a thin
bracket and a mapped block. It checks selected/opposite caps, internal bands,
both ends, rotations, scales, face tags, independent volume/topology and thread
identity. Polygonal profiles are faceted models; these are not curved STEP
qualification results.

The separate probes deliberately leave the recognized-sweep family:

| Geometry | Current observed result | Meaning |
|---|---|---|
| Tapered hollow duct | Accepted by the constrained Tet4 route; volume, classifications and boundary audited | One non-sweepable, multiply connected polyhedral solid is supported. |
| Pocketed block | Refused at boundary edge discovery: `NonStarEdgeCavity` | Valid input geometry reaches a recovery configuration outside the current accepted edge-cavity operation. |
| Branching junction | Same bounded `NonStarEdgeCavity` refusal | Junction geometry exposes the same important coverage gap. |

The pocket/junction fixtures use a voxel-union boundary with coplanar face
subdivisions, at target size 0.7 and meshing-work budget 500,000. They are
specific reproducible cases, not evidence that every pocket or junction fails;
a different boundary partition or resolution can change recovery behavior.

Both refused cases keep the destination unchanged. Their probe exit code zero
means the probe completed and checked rollback, **not** that a mesh was produced.
Their `published` value is false. The tapered-duct CTest gate explicitly requires
publication. Fixtures and numerical stage codes are retained in the source and
`evidence/rc10/frontier-final.json`.

## Priority 1 — robust constraint recovery

Make the pocket and junction fixtures accepted without dropping protected
edges/facets or filling voids. The next implementation should broaden valid local
recovery cavities, with explicit geometric witnesses; where a cavity cannot be
safely retriangulated, support a bounded refinement/recovery cycle. New boundary
vertices require coordinated curve/surface updates and a refreshed ownership
contract. A retry that simply discards a troublesome boundary is unacceptable.

Acceptance gate: published positive tets, exact classified boundary incidence,
independent volume/void checks, deterministic output, cancellation and finite
work/memory bounds across scales, translations and nearby perturbations. Retain
small failing fixtures whenever a new recovery configuration appears.

## Priority 2 — curved and imperfect CAD

The optional import adapter is not the same as broad meshing support. Add a
curved trimmed-surface corpus with seams, periodic charts, fillets, sliver faces,
small gaps and narrow passages. Edge discretization must be shared between faces;
curvature and proximity should drive refinement through one consistent sizing
policy. CAD projection, tessellation error, model units and tolerance handling
need explicit contracts. Offer repair/defeaturing choices with a change report;
never silently alter the represented solid.

Acceptance gate: quantified boundary error, preserved feature identities and
watertight meshes on a maintained public CAD corpus. Higher-order nodes require
curved Jacobian certification, not only midpoint insertion. Exercise the optional
OCCT path in CI and across supported platforms.

## Priority 3 — one spatial sizing policy across consumers

PAACS 1.0.0 supplies a shared bounded field interface, composable primitive-distance
controls, global bounds/factor and independent achieved-size diagnostics. Adaptive
boxes enforce whole-cell field bounds. Geometry Tet4 can apply boundary sizing
and opt into bounded interior refinement, but a tighter cube fixture stops at
`BoundaryRefinementRequired`. Coupling boundary subdivision and interior recovery
is now a concrete next gate. See [the field guide](SPATIAL_CONTROL_FIELDS.md).

Extend CAD face/edge selection and volume controls through arbitrary 3D domains,
with curvature/proximity options, bounded growth and a spatial query index.
Preserve the distinction between requested target, achieved spacing, and hard
acceptance limits. Add periodic/embedded constraints with authoritative ownership;
scalar fields alone cannot establish those topological relationships.

Acceptance gate: controls resolve the same geometric feature across supported
routes, independent of coordinate orientation, with bounded query cost and a
clear account of any target violation. Feature detection stays off by default
until enabled with angle/tolerance/minimum-size controls.

## Priority 4 — shared interfaces and fluid wall layers

Multi-material bodies need one authoritative interface mesh used by both sides,
consistent region labels and conforming junctions. Touching bodies currently
require an explicit hybrid plan; bounding-box overlap alone does not construct
such an interface.

For CFD, extend prepared straight prism layers to curved wall fronts, with
normal smoothing, concave-corner treatment, collision/proximity detection,
local layer termination and a conforming connection to a tetrahedral core.
Hex/prism/pyramid/tet transitions need shared topology. These are substantially
more demanding than extruding an already valid flat patch.

Acceptance gate: no layer intersections, positive Jacobians, exact core/interface
coverage, preserved wall tags and reported achieved first-cell height, growth,
layer count and total thickness. Solver wall-unit targets depend on flow inputs;
a geometric mesher cannot guarantee a target y+ on its own.

## Priority 5 — solver and release qualification

Keep quality improvement boundary-preserving, deterministic where promised,
and budgeted. Improve the low-quality tail as well as the average. Test structural
patch/bending problems and manufactured diffusion solutions before claiming
solver readiness; add controlled flow cases when the curved-layer/core pipeline
is available. Compare error versus degrees of freedom and assembly/solve cost,
not just element count or a visually attractive mesh.

Profile and reduce temporary allocations before advertising general-route speed.
The rc.11 small constant-field Tet4 fixture observes 191,938 allocations; the new
field compiler uses one, and queries use none. This points to inherited recovery
and source-meshing work as a profiling target, rather than field arithmetic alone.

Build a repeatable corpus dashboard: geometry coverage, refusal reason, boundary
error, minimum/percentile quality, meshing time, peak RSS and cancellation cost.
Use the same geometry, tolerances, target sizes, quality criteria, thread count
and machine across trials. Record failures as results. Add Windows/macOS
and supported compiler CI. Establish CPU scaling before moving kernels to GPU;
exact predicates, topology mutation and ownership need a deliberate accelerator
contract, not a blanket parallelization claim.

A useful next release gate is a published, reproducible corpus of mixed CAD parts
covering holes, pockets, junctions, thin walls, fillets and assemblies, with every
result either meeting explicit limits or returning a documented bounded refusal.
The first concrete blockers are the two recovery fixtures above. A general
all-hex guarantee is a separate research objective, not a prerequisite for a
useful general Tet4/Tet10 product.

See the [current meshing roadmap](ROADMAP.md) for product milestones and
the required qualification evidence.
