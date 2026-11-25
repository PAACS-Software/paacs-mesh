# Native STEP meshing

Build the optional `PAACS::occt` target with `PAACS_WITH_OCCT=ON`.
See the [automatic STEP build](BUILDING.md#automatic-step-build) and the
[installed CMake package](BUILDING.md). Native STEP meshing is a preview.

## Import and mesh

```sh
paacs-native-cad-preview part.step part-preview.json 0.05
paacs-cad-mesh part.step 2 part.paacsmesh --native-cad --report-json
```

The preview JSON contains display triangles and original CAD face tags.
The mesher uses the STEP geometry and topology directly. STEP coordinates,
sizes and geometric error limits are in millimetres.

Normal automatic Tet4 selection first considers verified whole-body primitives,
affine blocks and linear sweeps for compatible constant-size requests. Local
face/volume sizing and other eligible requests use native wall and gap point
placement. Bounded selection retains its historical dispatch. Hex requests and
conforming multi-body requests retain their existing routing. Normal Tet4
`--independent-parts` requests resolve AUTO separately inside each conforming
child, so verified primitives keep their established method and other eligible
parts can select the native front. Original body/face controls are mapped and
checked per part; they are not bypassed. Global work uses the smaller of the
meshing and validation ceilings before any child work is divided. Body identities,
interference reports and original-source/final Native audits remain mandatory.
Inspect `auto_dispatch_scope` and each body's `route_resolution`; wrapper sizing
metadata does not imply that every part used the same method. Explicit false,
explicit true and Bounded requests retain their existing eligibility and refusal
behavior. A child success does not qualify an incomplete assembly.
`--no-advancing-interior-front` selects the previous placement and recognized
primitive, planar or swept routes.
Check `archive_published`, `audit_valid` and each body's `route` in the JSON
receipt. A refusal leaves an existing output archive intact.

`--native-source-fallback` opts into an experimental provider switch after a
remaining named geometric chart/wire refusal. The existing live Tet4 construction and
all Native checks remain active; `--no-native-source-fallback` retains the default
disabled policy. It does not force the initial front or replace verified
primitive/mapped routing. Explicit local surface controls and construction caps
are not bypassed. Inspect each source attempt's provider, witness and configured
angle in `source_shell`; native tessellation is not an achieved-size certificate.
An opted-in live request that refuses its source stays refused rather than
silently switching to legacy star/partition construction.
The three-attempt source bound includes provider and chord retries together
(Bounded explicit opt-in uses at most two; default Bounded retains one).

`--boundary-first` separately opts into recovering the source boundary before
interior insertion in a compatible live spatial Tet4 route. It explicitly
overrides the order of an enabled initial front; default `--no-boundary-first`
keeps the existing order. It does not enable the source fallback, force the front
on, or change verified primitive/mapped choices. Inspect
`requested_boundary_first`, `initial_front` and `interior_front` for the request
and actual stages rather than inferring success from the flag.
`material_classification` reports returned-batch attempts even when a later
stage is interrupted. A Native query is counted when started; an interrupted
attempt need not yield a decision. `profiled_batches` and `ray_witness_attempts`
identify measured strict serial boundary-first work. The setup, Native-query
and ray-witness seconds are nested within material/volume time; zero profiled
batches means unmeasured, not zero classification cost.
Positive original trimmed-CAD ray witnesses are available
for strict serial boundary-first classification; ambiguous points still use the
ordinary native classifier. Inspect actual admission/query counters, and retain
all final checks. These nested phase times are part of volume/total time.
See [source preparation qualification](SOURCE_PREPARATION.md) for the corrected
sampling/correspondence rules and the current limits of the experimental routes.
Explicit interior site/candidate caps remain upper bounds for boundary-first
insertion; its current post-recovery builder also retains a 1,024-site ceiling.

## Native size controls

A single Tet4 solid can use selected CAD faces and volume boxes together.
For the checked-in cylinder, face 1 is the curved wall.

```sh
paacs-cad-mesh tests/fixtures/cad/cylinder.step .7 cylinder.paacsmesh \
  --native-cad --face-size 1 .38 --size-box -.3 -.3 .7 .3 .3 1.3 .31 \
  --min-size .34 --growth 1.55 --report-json
```

- `--native-sizing` selects the native size-field route for a face-only
  request. A box selects it automatically. Use `--face-size TAG SIZE` and
  `--size-box X0 Y0 Z0 X1 Y1 Z1 SIZE` to add controls.
- `--min-size` clamps the requested field. CAD curvature and boundary recovery
  can require smaller edges. The field uses the smallest active size.
- `--growth G` controls how quickly the requested size increases with distance
  from a selected face or box. Inspect achieved lengths and adjacent Tet4 size
  ratios when the transition matters.
- Recognized cylinder and tube cap profiles remain available through
  `--face-size` with `--no-advancing-interior-front` on the analytic route.
  Native size-field controls currently
  apply to one Tet4 solid.

Use [native_step_sizing.py](../scripts/native_step_sizing.py) to save a face
tag and box with the STEP input hash. Its `mesh` command reopens the project,
runs the native CLI and verifies face and body groups in a C3D4 solver include.
The `plot` command uses exact intersections of the stored Tet4 cells.

## Recovery, assemblies and limits

| Option | Effect |
|---|---|
| `--recovery-mode advanced` | Try bounded native edge and face recovery on routes that reach a certified shell. The default is `standard`. |
| `--local-partition-recovery`, `--adaptive-tet-recovery` | Enable optional bounded recovery routes. |
| `--curvature-points-per-turn N` | Enable curvature-based surface sizing on eligible native routes. |
| `--improve-tet-quality` | Try bounded quality improvement on native pieces and after spatial sizing, keeping the prescribed surface fixed. |
| `--initial-topology-quality`, `--no-initial-topology-quality` | Explicitly enable or disable severe numerical-topology cleanup before recovery; unset follows Normal enabled / Bounded disabled. This does not force the front or certify Delaunay topology. |
| `--advancing-interior-front` | Force wall and gap initial point construction for one supported native Tet4 solid. Unsupported combinations return a diagnostic. |
| `--no-advancing-interior-front` | Disable the new default and use the previous construction and route selection. |
| `--independent-parts` | Mesh supported bodies separately and report material interference. The result certifies each part individually. |
| `--max-work`, `--max-working-bytes` | Set work and temporary storage limits. |
| `--warning-quality`, `--min-quality` | Report a weak mesh or require a strict quality floor. |
| `--trace-native-stages` | Record active stages and diagnostics on standard error. |

Normal disables the practical logical work cap; Bounded retains 50 million
units for construction and validation. `--max-work` sets both ceilings; an explicit limit is retained
across automatic fallback attempts. Work units measure bounded operations,
not elapsed seconds.

With `--no-advancing-interior-front`, for spatial Tet4 sizing,
`--improve-tet-quality` also supplies interior sites to the first constrained
Delaunay construction. Enclosed boxes contribute points
before the fill, and oriented surface triangles propose inward points using the
local size field. Spacing and native-solid checks admit these proposals. This
avoids relying on small later insertions to replace long connections inherited
from a boundary-only fill; it does not guarantee a requested quality floor.

The subsequent bounded quality pass first tries to repair critical cells with
face and edge swaps. If direct swaps cannot replace a nonconvex edge star, it
can reconstruct a bounded connected cavity from an existing nearby vertex.
The replacement must preserve the oriented shell, every vertex, protected
constraints and body ownership, and pass quality and native-solid checks.
Partial interior-node movements and ordinary edge swaps then improve spacing.
These stages preserve CAD boundary coordinates and tags, and retain
the native geometry and validity checks. Placement is limited by node, work,
query, temporary-storage and output admissions. If the extra-site seed is refused,
the ordinary boundary seed is retried within the remaining allowance; both
attempts count toward work and time. Trace stages include `interior_sites`,
`interior_seed_fallback` when needed, and `volume_quality`. Inspect the achieved
size and quality distributions, since bounded runs can leave weak cells unresolved.

Automatic method selection records its capability proof and charged attempts in
`route_resolution`. It never selects a structured method that ignores face,
box or spatial controls. The original Native gates still certify the output;
failed methods and any fallback share the caller's remaining budgets.
`--advancing-interior-front` forces the native size-field route;
`--no-advancing-interior-front` disables it. The last explicit flag wins.
In C++, `NativeCadMeshOptions::useAdvancingInteriorFront` is `std::nullopt` for
automatic selection, `true` to force it and `false` for the previous route.
Automatic selection respects disabled native fallback, requested geometry
strategies, body controls and other unsupported construction contracts.

For method comparisons, set `--max-chord-error` and
`--max-relative-volume-error` explicitly (`--max-volume-error` is retained as
an alias). Automatic chord defaults differ by construction method; inspect the
body `fidelity` receipt and measured CAD errors. A source retry may use a tighter
chord but never a looser requested volume bound. A structured fallback to the
general front retains the primitive chord limit rather than silently increasing it.
If the automatic route refuses a surface before initial point construction has
started, it can try the previous route within the remaining work allowance.
The receipt records this in `initial_front_fallback`, including the first
failure and its work and time. Both attempts contribute to total cost.
Cancellation, exhausted budgets and strict quality refusals are terminal.
An explicit `--advancing-interior-front` request does not use this fallback.
It constructs initial interior points **before** constrained volume
meshing, including for a uniform global-size request with no face or box controls.
A first sweep constructs wall sites from the fixed surface triangles. Step
lengths follow actual neighbouring edges and grow toward the requested size.
A persistent Delaunay queue then proposes sites in the largest remaining
interior gaps, relative to the requested size, before resuming the exposed-face
front. Existing points can be reused where
advancing regions meet. These stages construct the initial point cloud; they
do not add a refinement region or change the prescribed size field.

Initial sites use their requested size-field values for separation from other
interior sites. A small launch step near the surface does not permanently give
an interior point a smaller exclusion radius. Distances to fixed boundary nodes
retain the local surface scale, so the original fine boundary can be respected.
For a uniform request, initial interior pairs remain at least `0.65 * size`
apart. Constraint recovery can subsequently add geometry-required points; this
initial-site rule is not a minimum-edge guarantee for the final mesh.

Fresh positions are ranked using neighbour separation and the actual connected
conflict cavity in a persistent auxiliary Delaunay triangulation. Cavity discovery
can visit the complete admitted complex; a small fixed cell count must not bias
which point positions are considered. Work and storage still bound the search.
The auxiliary constructor and subsequent constraint-seed constructor check work,
temporary payload, cell limits and cancellation during construction. Their
internal retries share each constructor's remaining allowance.
The selected patch must pass topology and embedding checks before any insertion.
For the owned auxiliary complex, an exact local shell proof certifies a complete
positive cone fan and its contacts with the surviving mesh. Unsupported patches
retain the full pair and survivor checks. This avoids quadratic comparisons
between cones that share the inserted point. Exact source-surface checks, native
point admission and original body identity are retained. Cached material states
can be inherited only when the replacement preserves the same union and no
internal original source facet is erased.

Gap candidates use the current tetrahedral circumcentres, source-surface
clearance and all-known-point separation. Exact ray parity against the original
surface and strict native admission reject cavities and exterior points.
The queue updates only cells changed by an insertion, discards retired entries
lazily and rebuilds after compaction. Its storage is reserved inside the existing
working allowance; it does not trigger a full triangulation after each point.

Retired auxiliary tetrahedra are reclaimed between candidate batches. Compaction
copies the live connectivity and cached material classification without moving
points or retriangulating. It keeps the admitted storage capacity and stages the
replacement within the same working allowance. Old candidate plans are invalid
after the cell indices change; a failed staged rebuild leaves the old mesh intact.

Only the accepted point coordinates enter ordinary constrained construction.
Temporary front cells and the auxiliary triangulation are not published material
meshes or coverage certificates. Rebuilding and recovering constraints can change
Delaunay connections, so initial proposal quality is not a final quality guarantee.
The normal native validity and CAD-error checks still control publication.

The body receipt's `initial_front` object records initial sites, reused apices,
propagated faces, complete-cavity sizes, certification work, queries and elapsed
time. A zero `cavity_limit` means automatic discovery of the complete admitted
cavity. Propagation depth describes ancestry, not distance from the CAD surface.
`wall_nodes` and `void_nodes` distinguish the first sweep from gap placement;
closure may supply additional sites. `source_sweep_complete` records whether
every original source face was visited. `void_queue_exhausted` is not a coverage
certificate: candidates can be refused by geometric or quality guards.
`maximum_unresolved_void_ratio` also includes historical refused candidates
that later insertions may have repaired.
Storage rejections, compaction counts and reclaimed-cell counts distinguish
memory pressure from geometric refusals. `compaction_peak_modeled_bytes` estimates
the auxiliary payload; it is not a process-memory measurement.
`budget_limited` or `stalled` can retain admitted sites with unresolved construction.
Inspect the achieved sizes and quality; there is no certificate for every edge
length, growth ratio or global minimum quality.

A uniform request does not run the later adaptive interior stage. Explicit face
or volume sizing controls can additionally activate it; `interior_front` then
reports its separate work and unresolved targets. Volume controls supply a size
field rather than a Cartesian point grid. The optional `--improve-tet-quality`
pass remains separate and is not required to run initial point construction.

The default material checks reject exterior tetrahedron centroids. On eligible
Tet4 recovery routes, `--allow-centroids-within-chord-tolerance` applies a
bounded CAD distance check and additional interior samples. Check the receipt
for chord and volume error, body ownership, and whether all native samples
were checked. Sampled CAD error describes the checked locations.

The [release validation](RELEASE_READINESS.md) gives the matched corpus results.
The [native CAD development notes](NATIVE_CAD_HISTORY.md) describe the recovery
experiments, evidence and detailed option behavior.
