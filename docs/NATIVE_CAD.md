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

The CLI chooses a mesh route based on the solid. It can use a recognized
primitive, a supported planar or swept strategy, or native Tet4 recovery.
Check `archive_published`, `audit_valid` and each body's `route` in the JSON
receipt. A refusal leaves an existing output archive intact.

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
  `--face-size` on the analytic route. Native size-field controls currently
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
| `--advancing-interior-front` | Select experimental size- and shape-guided interior construction for one native Tet4 solid. |
| `--independent-parts` | Mesh supported bodies separately and report material interference. The result certifies each part individually. |
| `--max-work`, `--max-working-bytes` | Set work and temporary storage limits. |
| `--warning-quality`, `--min-quality` | Report a weak mesh or require a strict quality floor. |
| `--trace-native-stages` | Record active stages and diagnostics on standard error. |

For spatial Tet4 sizing without the advancing interior front,
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

`--advancing-interior-front` selects the native size-field route for a single
Tet4 solid. It constructs initial interior points **before** constrained volume
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
