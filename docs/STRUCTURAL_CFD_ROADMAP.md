# Structural and fluid meshing: proposed development roadmap

Status: 1.0.0. These are development priorities, not claims of general
qualification. Balanced primitive strategies now cover cylinder and solid
circular sweeps. Opt-in surface-edge analysis is implemented, with angle and
length controls; it does not yet recover general CAD features or alter sizing.
Composable spatial fields and explicit distance-based region/face controls and conforming adaptive box Tet4
meshes are now available. Planar structured/unstructured prism/hex stacks
preserve prescribed wall layers and fill a graded interior. See
[local controls](LOCAL_SIZE_CONTROLS.md) and [layers](PRISM_LAYER_CONTROLS.md).

The next bounded implementation targets are:

1. Apply the new physical sizing/growth contract to general CAD-owned surfaces
   and volumes, with feature chains/patches and opt-in local variation.
2. Connect planar layers to a separate tetrahedral core, then qualify curved
   wall-front collisions, concave junctions and local layer termination.
3. Extend uniform layouts to annular/tapered solids; add independent
   tetrahedron dihedral/conditioning and finite-volume face metrics,
   followed by a first actual structural and flow solver convergence case.

The current release keeps feature analysis off by default and never disables
geometry preservation, Jacobian checks, or the boundary audit. See
[SURFACE_FEATURES.md](SURFACE_FEATURES.md) and
[UNIFORM_CIRCULAR_SWEEPS.md](UNIFORM_CIRCULAR_SWEEPS.md).

## What uniformity should mean

The new cylinder option targets similar cell volumes and balanced radial/axial
spacing. It retains curved-boundary accuracy controls, positive-cell checks,
full boundary incidence auditing, and Tet4/straight-sided Tet10 choices. It
does not make every edge equal. Curvature tolerances can require cells smaller
than the requested nominal size. Thin disks can still contain stretched cells.

Uniform isotropic tetrahedra are useful for some structural and bulk-flow
regions. They are not an appropriate universal objective for thin structures,
viscous wall layers, wakes, shocks, or strongly directional solution gradients.

## Priorities and measurable acceptance gates

| Order | Upgrade | Structural benefit | Fluid benefit | Evidence required before calling it ready |
|---|---|---|---|---|
| 1 | Qualify general curved CAD and feature recovery | Fillets, holes, narrow ligaments, seams, contact faces | Watertight fluid cavities, pipes and intersecting passages | Curved STEP/BREP corpus with cavities, multiple bodies, trim seams, tiny edges and bad inputs; independent CAD distance and inside/outside checks; preserved entity IDs; explicit refusals rather than silently discarded geometry |
| 2 | Unified size-field and growth controls | Refine stress concentrations, enforce elements through thickness | Refine walls, gaps, wakes and shear layers | One size contract across strategies; curvature/proximity limits; bounded adjacent-size transitions; requested versus achieved size report; convergence and memory-budget tests |
| 3 | Solver-specific quality reports | Jacobian bounds, dihedral angles, conditioning, slivers | Face non-orthogonality, skewness, cell-volume jumps and interpolation weights | Published metric definitions, independent oracles, per-region distributions and offending IDs; separate validity from suitability; configurable thresholds rather than one universal quality number |
| 4 | Robust prism/hex wall layers with junction handling | Efficient thin solids and layered materials | Viscous near-wall resolution with specified first-cell height and growth | Extend existing straight/curved layer components to concave corners, opposing walls, layer collisions and core transitions; report achieved thickness and layer count per patch; no inverted cells or hidden gaps |
| 5 | Curved higher-order elements | Better geometric fidelity for Tet10 and high-order analysis | Curved-wall high-order discretizations | Project shared edge/face nodes through CAD ownership; optimize or roll back; whole-reference Jacobian certificates after every accepted change; curved seams/contact tests and solver convergence |
| 6 | Local quality improvement and adaptive refinement/coarsening | Resolve local stress/strain error with fewer DOFs | Follow wakes, shocks and thermal gradients | Bounded split/collapse/swap/relocation work; preserve CAD features and material interfaces; monotone accepted objective changes; conservative transfer where needed; topology/Jacobian checks and rollback |
| 7 | Conformal multi-region and periodic meshes | Bonded interfaces, material regions, explicit contact pairs | Conjugate heat transfer, rotating/periodic sectors | Matching interfaces or explicitly represented nonmatching pairs; correct orientations and tags; periodic node/face maps; round-trip tests; no duplicate interface volume |
| 8 | Solver interoperability | Practical FEA model import with element ordering and sets | Boundary patches and volume regions accepted by CFD solvers | Select and document supported VTK/MSH/solver formats; ordering/orientation and high-order round trips; preserve units/tags; run an actual downstream solver, not just its file reader |
| 9 | Accuracy-based qualification suite | Patch tests, cantilever bending, thick-cylinder stress, plate with a hole | Poiseuille flow, manufactured solutions, laminar cylinder flow, heated pipe | At least three refinements; displacement/energy/stress or velocity/pressure/flux errors; mesh independence; compare cost at equal solution error, not just equal cell count |
| 10 | CPU scaling, then portable accelerator backends | Larger models within memory budgets | Large fluid domains and repeated adaptation | Profile generation, recovery, sorting, audit and high-order certification separately; compact storage, deterministic partitioning and bounded scratch space; 1/2/4/8-thread measurements; platform CI and race checks before GPU claims |

Cross-platform qualification is a release gate throughout: add Linux GCC/Clang,
Windows MSVC and macOS arm64 jobs with explicit supported compiler versions,
installed-consumer tests, strict floating-point behavior, cancellation and large
index checks. Keep C++26 as the requested reference and make any compatibility
mode explicit; language support alone does not prove platform portability.

Build the solver verification harness alongside priorities 1–3; do not wait
for all meshing features before testing solution accuracy. The ordering above
expresses capability priorities, not a requirement to finish each row serially.

The current code already contains size controls, layer components, quality
checks, Delaunay/recovery algorithms, and high-order machinery. The work is to
integrate and qualify them for these use cases, not to assume that having a
component proves an end-to-end production capability. The optional qualified
STEP/BREP import route is still limited to planar faces and straight edges.

## How the metrics differ

For structural tetrahedra, preserve positive full-cell Jacobians, inspect shape
and sliver measures, and test the intended element formulation. A straight
Tet10 increases interpolation order but does not itself improve the faceted
geometry. Stress near a mathematical singularity is not a reliable target for
pointwise mesh convergence; choose a defined quantity of interest away from
the singularity or an appropriate integral quantity.

For finite-volume CFD, a tetrahedral mean-ratio score does not describe face
non-orthogonality or skewness. OpenFOAM's documented controls explicitly include
these separate measures [1]. Its layer controls specify patch selection,
layer count, expansion and thickness [2]. First-cell height should follow the
flow regime and wall treatment; a generic mesher cannot promise an appropriate
y-plus without flow/viscosity/wall-shear information. Intentional anisotropy
must be judged in its requested metric, not rejected solely for a high Euclidean
aspect ratio.

## Where a defensible distinction could emerge

**Research hypothesis:** combine geometry-specific fast paths, a general
fallback, explicit memory/work budgets, and auditable validity guarantees with
solver-driven adaptation. Let the planner explain which strategy was chosen,
what geometric or quality constraint forced refinement, and the expected cost.

The useful result would be fewer failures and lower time/memory at a specified
geometry tolerance and solver error. MFEM documents target-matrix mesh
optimization and simulation-driven adaptation [3,4]. Novelty needs a literature
review and an experimentally supported contribution.

An incremental audit after a local edit is a promising performance direction,
but it needs a sound invalidation rule and retained global invariants. Keep the
full audit as an oracle and release gate until differential tests justify the
incremental path. Never trade a silently incomplete check for a headline speed.

## Comparison protocol

1. Publish the geometry corpus, units, element order, CAD tolerance, requested
   size/metric field, feature policy, quality thresholds and allowed budgets.
2. Compare pinned builds under the same CPU, thread and memory limits. Include
   difficult and failed cases. Match achieved size, boundary fidelity and
   solver error where possible.
3. Record end-to-end time, per-phase time, peak process memory, output size,
   cell/node counts, metric distributions and downstream error. Retain all raw
   repeats, seeds and failure reasons. Report trade-offs rather than one winner.
4. Separate geometry generation, meshing, validation and file I/O. Include
   equivalent validity checks in the total comparison. Keep acceleration tests
   separate from reference correctness tests, with deterministic fallback.

For accelerators, first identify large, independent kernels such as quality
evaluation and batched geometric queries. Retain CPU fallback for ambiguous
predicates and irregular recovery. A portable backend does not establish
performance or correctness on every GPU: qualify specific vendors, drivers,
precisions and transfer costs, including integrated GPUs and small workloads.

## Primary references

Accessed 2026-09-30. These sources describe existing capabilities; the roadmap
and acceptance gates above are project recommendations.

1. [OpenCFD: Mesh quality](https://doc.openfoam.com/2212/tools/pre-processing/mesh/generation/snappyhexmesh/meshquality/).
2. [OpenCFD: Layer addition](https://doc.openfoam.com/2212/tools/pre-processing/mesh/generation/snappyhexmesh/layers/).
3. [MFEM: Meshing miniapps](https://mfem.org/meshing-miniapps/).
4. [Dobrev et al., Simulation-Driven Optimization of High-Order Meshes in ALE Hydrodynamics](https://arxiv.org/abs/2001.11536).
