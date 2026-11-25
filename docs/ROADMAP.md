# PAACS Mesh roadmap

PAACS will qualify meshes for structural, thermal and fluid workflows through
reproducible CAD and solver studies. The priorities below extend the current
mapped, swept and tetrahedral routes. Release support and evidence are recorded
in the [release policy](RELEASE_POLICY.md) and [validation report](../VALIDATION.md).

## Current starting point

| Area | Implemented scope |
|---|---|
| Geometry | Explicit primitives, recognized affine blocks and polygonal sweeps, and bounded Tet4 for supported polyhedral models. The optional OCCT adapter provides native CAD queries and meshing routes with their own documented coverage. |
| Volume elements | Tet4, straight-sided Tet10, Prism6, mapped and swept Hex8, and Pyramid5 transitions on supported routes. |
| Wall layers | Graded stacks from prepared planar triangular, quadrilateral and mixed patches. A separate prepared curved triangular-front component supports normal smoothing, obstacle clipping and thickness backoff. |
| Higher order | Shared midside nodes for straight-sided Tet10. A bounded provider-driven curving API supports shared-node placement and whole-element Tet10 Jacobian certification. |
| Controls and output | Composable spatial sizing, profiles, layer schedules and achieved-size reports on supported routes. Native archive v1 and legacy VTK and MSH 4.1 ASCII writers are available. |

See [element choices](ELEMENT_CHOICES.md), [native CAD](NATIVE_CAD.md),
[prism layer controls](PRISM_LAYER_CONTROLS.md), [spatial fields](SPATIAL_CONTROL_FIELDS.md)
and [archive format](ARCHIVE_FORMAT.md).

## Diagnostic development milestone

The unreleased [inspection interface](MESH_INSPECTION.md) adds opt-in Tet4 mean
ratio, interior dihedral and signed scaled-Jacobian measurements; pagination of
every native cell; spatial distributions; exact trimmed-CAD centroid distances;
quality-field VTK; and deterministic whole-element cutaways and x-ray views.
Comparisons use a common world grid and preserve validation provenance, units,
source hashes and mesh-local cell IDs. These measurements distinguish shape,
geometry, achieved sizing and solver adequacy.

The [cylinder refinement study](REFINEMENT_FINDINGS.md) establishes the next
construction work: report unfinished surface-spacing work, allocate bounded
insertion allowances from the surface metric, and avoid starving bulk initial
construction when a fine wall consumes proposal budgets. Investigate the
near-boundary slivers separately. Preserve all native validity gates and compare
quality distributions, achieved spacing, CAD error, work and time on fixed and
held-out geometries before changing defaults.

Next diagnostics work includes complete local size/growth certificates,
source-and-final-attempt-paired workflow retry totals, additional element-family
metrics, and solver-provided error indicators. A poor shape metric alone is not
a reason to automatically refine. Extend bounded local editing only with
validated rollback and preserved material/CAD ownership.

## 1. Broaden CAD coverage and apply sizing consistently

Qualify the complete path from shared edge sampling through trimmed face
meshing, volume recovery and native CAD checks. Include spline and periodic
faces, holes, pockets, branching passages, thin walls, small gaps and touching
assemblies. Retain body and feature IDs. Build each conforming shared interface
once and use it on both sides with opposite winding. Preserve independently
meshed bodies and report material interference when the selected assembly
policy permits it. Report any geometry repair and its measured displacement.

Apply global size, selected CAD face and edge sizing, volume controls, curvature
and proximity through a common physical-size policy. Measure achieved spacing,
adjacent growth and geometric error. Keep targets, warnings and hard limits
distinct in the report.

**Qualification.** Publish fixed and held-out STEP cases with original-body
identity, boundary and void checks, CAD deviation, achieved size, runtime,
peak memory and bounded failure diagnostics.

### Explicit CAD preparation (planned)

Add repair and defeaturing as separate operations, both disabled by default.
Ordinary import, meshing and recovery must not silently enable either one.
Repair addresses diagnosed representation defects under explicit tolerance
limits. Defeaturing intentionally removes selected features, such as small
fillets, chamfers or unwanted holes. Neither operation guarantees better mesh
quality or simulation accuracy.

Preserve the imported baseline and perform preparation on an owned copy. Return
a validated prepared model with a change report and retained, modified, split,
merged, deleted and new entity history. Protect requested body identities,
placements, interfaces, named selections, sizing controls and analysis-critical
holes, gaps and thin walls. A deleted or ambiguous protected selection must be
reported rather than rebound to a nearby face.

Start with explicitly selected features and a before/after preview. Record
original/prepared geometry and option hashes, units, validity, topology counts,
volume changes, tolerance changes, measured geometric deviation, unavailable
checks, work and time. Validate meshes against their prepared geometry and
report deviation from the original separately. Modified variants do not count
as successful meshes of the unchanged source.

The current strict importer rejects invalid CAD. Repair therefore needs a
separate bounded preparation path for raw input; it must return a valid model
before that geometry can enter the mesher. Keep the current import and mesh
validity gates in place.

**Qualification.** Use owned fixtures with short edges, sliver faces, inconsistent
edge parameters, small gaps, fillets and holes. Check unchanged baselines,
protected features, assembly identity, truthful entity history, rollback and
resource limits. Compare prepared and unchanged variants at measured achieved
sizes, with quality distributions and geometry changes reported together.
Continue improving surface and volume construction on unchanged geometry.

## 2. Build Prism6 wall layers with conforming transitions

Connect the prepared wall-layer components to CAD-tagged surface patches and
a separately generated Tet4 core. Preserve first-layer height, growth and
patch labels through concave and convex corners, narrow gaps and opposing
walls. Resolve collisions with bounded clipping, termination or a reported
refusal. Construct conforming Prism6, Hex8, Pyramid5 and Tet4 transitions where
the interface topology calls for them.

**Qualification.** Check positive full-cell Jacobians, layer intersections,
matching core interfaces and achieved first-layer height, thickness and layer
count by patch. The [hybrid straight-duct example](../examples/flow_hybrid_duct/README.md)
checks four wall-normal stacks, a conforming tetrahedral core and laminar
pressure-drop convergence on a prepared duct. Extend validation to curved
walls, junctions and automatically generated core interfaces. Determine
first-layer requirements from the flow regime and wall treatment.

## 3. Extend Hex8 coverage and develop Hex27

Expand reliable Hex8 meshing through sweep recognition, compatible end-face
meshes and multi-block decomposition. Cover swept bends, holes, branched
parts and interfaces where a compatible block layout can be established.
Develop conforming hex-dominant transitions for the remaining volume.
Evaluate the existing Tet4-to-Hex8 conversion separately on its supported
single-region input. Automatic all-hex meshing of broader CAD remains research.

Introduce Hex27 on verified Hex8 layouts. Define the triquadratic node order,
share edge and face nodes across adjacent cells, and extend storage, archives
and writers. Start with affine blocks and progress to curved mappings.
Any additional higher-order wedge or pyramid families need the same topology,
Jacobian and interoperability checks.

**Qualification.** Check whole-cell Jacobians, shared high-order faces and
interfaces, CAD error, solver ordering and convergence on mapped and swept
cases. Require an explicit result for layouts the block strategy cannot cover.

## 4. Qualify CAD Tet10 and higher-order geometry

Integrate the Tet10 placement provider with native CAD ownership. Project shared
midside nodes through the correct curves and trimmed faces across seams and
interfaces. Try bounded high-order quality improvement when a projection
distorts a cell. Report straight-sided and curved output separately.

**Qualification.** Check positive Jacobians over the complete reference
element after curving, measure boundary deviation, and independently verify
shared-node topology and downstream high-order connectivity. Record whether a
refusal came from geometry, an invalid element or exhausted certification work.

## 5. Enable bounded local mesh editing

Extend local insertion, splitting, swapping, relocation and coarsening to
operate on selected regions while preserving CAD ownership and material
interfaces. Support solver-provided size or error indicators and retain a map
between original and edited entities for field transfer.

**Qualification.** Certify each affected patch before committing it. Verify
conservation for supported field transfers, reproducibility under fixed
settings, rollback on refusal, and agreement between local checks and full
mesh audits on regression cases.

## 6. Prove interoperability and simulation accuracy

Qualify a downstream structural and a flow-solver export workflow. Verify
units, element ordering and orientation, named boundary and material sets,
shared interfaces and high-order connectivity after import. Extend the native
archive and exchange formats as those workflows require.

Build from the [cantilever](../examples/elasticity_cantilever/README.md),
[plate-with-hole](../examples/elasticity_plate_hole/README.md) and
[hybrid straight-duct](../examples/flow_hybrid_duct/README.md) convergence studies.
Add curved Tet10, thick-cylinder stress and manufactured diffusion before
moving to more complex CFD cases. Across at least three mesh levels,
report solution error or quantity-of-interest convergence, degrees of freedom,
mesh generation and solve time, and memory. Include family-specific metrics
such as Jacobian and dihedral distributions for structural meshes and face
skewness, non-orthogonality and wall-normal spacing for finite-volume meshes.

## 7. Extend the mesh and model workflows

Qualify standalone curve and surface meshes, then add shell and beam solver
exports with explicit thickness, section and orientation metadata. Build on
CAD ownership and shared-interface checks to support periodic node mappings,
embedded features and independently meshed contact pairs. Add partitioning
with stable global IDs, ghost entities and named groups for parallel solves.

Expose CAD selection, size assignment, background meshing, section inspection
and local remeshing through the application layer. Save the CAD source reference,
settings, selections and mesh together so the workflow can be reopened and
reproduced. Extend adaptive remeshing and verify how field transfer affects
the solver's quantity of interest.

## 8. Distribute the engine through Python

The first [Python interface](PYTHON.md) wraps native STEP import, audited Tet4
meshing and tagged export, with combined face and box sizing, minimum size and
growth on the supported single-solid route. Published RC3 wheels support
CPython 3.12–3.15 on Windows and Linux x64/ARM64 and macOS ARM64, with macOS 15
as the deployment floor. See the [compatibility matrix](PYTHON_COMPATIBILITY.md)
for independent build and installation qualification. Optional NumPy
snapshots expose node coordinates,
typed connectivity and CAD classifications without modifying the native mesh.
Expose additional native controls
and cancellation without introducing a separate meshing implementation.

Qualify future interpreter and platform additions independently.
Keep CAD dependencies optional in the C++ build and bundle the required runtime
in STEP wheels. Publish wheels to PyPI after installation qualification.
CAD inspection and source-bound semantic selections expose named groups and
explicit geometric predicates, preserving original body identities and
reporting ambiguous selections. Bounded inspect–select–mesh–audit workflows
capture reproducible recipes, explicit retry controls and attempt reports.
Headless tagged previews support visual selection checks; project publication
preserves source, intent, recipe and tagged exports together. Extend these
functions with a thin automation adapter and richer interactive clients.

**Qualification.** Install each wheel in a fresh virtual environment without a
compiler or manually configured library paths. Compare its meshes, audits and
face/body tags with the native engine on the same inputs and settings.
