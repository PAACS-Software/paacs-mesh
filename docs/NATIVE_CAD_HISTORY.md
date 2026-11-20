# Optional native CAD adapter

This local development extension uses `PAACS::occt`, built with
`PAACS_WITH_OCCT=ON`. The portable `PAACS::mesh` target remains independent of
OCCT. It adds new APIs without reintroducing retired topology routes or the
previous V1 model format. It is separate from the original 1.1.0 qualification.

See the dated development evidence below for the latest checks. The earlier
[CAD coverage iteration and images](LOCAL_CAD_COVERAGE_2026-10-03.md) and
[native validation](LOCAL_NATIVE_CAD_VALIDATION_2026-10-03.md) remain historical evidence.

Paths under `evidence/` and `out/` refer to private development records
retained outside the source package.

## Native Tet4 recovery mode

`--recovery-mode standard` is the default. Use `--recovery-mode advanced` to
try bounded source edge and triangle recovery in the full PAACS Delaunay
complex when the native Tet4 fallback reaches a certified surface. The
existing planar, primitive, and other native routes can still finish first.
In C++, set `NativeCadMeshOptions::recoveryMode` to `NativeCadRecoveryMode::Advanced`.

```powershell
paacs-cad-mesh part.step 2 part.paacsmesh --native-cad --recovery-mode advanced --native-validation-threads 4 --report-json
```

The advanced step keeps the surface coordinates fixed. It carves material
using the orientation of recovered source triangles and withholds cells that
conflict with missing triangles. The resulting seed goes through native
material classification, the existing closure and contact checks, and the
final native CAD and volume gates. The report's `constraint_recovery` section
records the stage result, missing constraints, accepted patches, seed size,
work, and elapsed time. A seed report describes an intermediate result.

The live recovery mesh retains reciprocal tet neighbours after each accepted
patch. Missing edges use an exact transverse walk, with indexed discovery for
ambiguous contacts. Face recovery follows neighbours from the target's edge
stars. The edge-cavity solver receives just the affected vertices and can
stage one Steiner point per replacement, within the shared work budget and
a 64-node ceiling. Original source nodes and triangles remain fixed during
this stage. Carving exports only the staged nodes used by material tetrahedra.
Those nodes require strict native material admission before publication.
The JSON report separates walks, indexed fallbacks, candidate attempts,
inserted Steiner nodes and Steiner nodes retained after carving.

When the coarse partition of a single solid cannot find a star kernel,
advanced mode can try direct boundary recovery before the remaining partition
strategies. This trial uses at most one quarter of the remaining modeled work
and 70 million work units overall. Its constraint stage has a 45 million unit
cap. The partition strategies can still run after the trial refuses a mesh.
A strict minimum element quality keeps the established route order.

The final recovery route shares the meshing work and memory budgets and can
spend up to a third of the remaining modeled work on constraints, capped at 450 million work units. Raise
`--max-work` when the advanced stage reports that its allowance was exhausted.
If it returns no seed,
the established front can run with the remaining budget. Assembly policy,
strict quality controls, and the explicit centroid chord allowance remain
independent options. `--no-native-fallback` disables the native fallback
including this stage.

`--reuse-native-queries` also reuses the OCCT distance solver's target-face
setup on consecutive queries to the same imported face. Each point still
runs the trimmed-face distance solve and the same acceptance checks.

`--local-partition-recovery` enables a bounded Tet4 trial when a coarse split
leaves a piece with a handle. It constructs one partition, recovers the
tetrahedra inside its pieces, and refills any cavities exposed by the original
CAD material checks. The shared interface triangulations stay fixed. The
completed mesh must pass the original-solid and interface checks before
publication. `--native-validation-threads` can accelerate larger centroid
batches. The C++ option is `NativeCadMeshOptions::useLocalPartitionRecovery`.

`--adaptive-tet-recovery` enables experimental route ordering and local
refill. After a coarse star-kernel refusal on a spline-rich solid, bounded
Cartesian partitions run before expensive speculative whole-solid fill. A
whole-solid star refusal can also start an early five-division trial at the
configured chord limit. The published mesh must pass the same measured
original-face chord and material checks. A local refill can continue from
retained tetrahedra and search a wider cavity
within the original work cap. After particular partition budget refusals the
existing concave-feature split can run earlier. Original-CAD material, volume,
boundary and body-ownership checks still gate publication.

When two exterior triangles on the same original CAD face bound a
material-rejected tetrahedron, the adaptive partition route can use the
opposite diagonal exposed by the retained tetrahedra. It checks the new
facets against the CAD and certifies the revised shell before refilling it.
Original CAD edges and shared partition cuts stay fixed. On contact
assemblies, the adaptive route can attempt a bounded recovery on an
individual non-star piece before splitting the whole component again.

With advanced recovery enabled, a single-body surface that fails shared
curve conformity can also retry with curvature-based edge sampling.
Every affected face uses the same edge samples. Its trimmed chart is
triangulated using a metric derived from the physical CAD surface, then
checked for chord error and geometric embedding.

The adaptive mode can reuse exact prior strict serial centroid checks when
the original BREP, body, tolerance and point bits match. With
`--reuse-native-queries`, it also memoizes successful trimmed-face distance
answers under an exact native face and point key. The C++ option is
`NativeCadMeshOptions::useAdaptiveTetRecovery`. Assembly policy remains
explicit through `--independent-parts` and its interference report.

`--improve-tet-quality` tries an optional quality pass on small native Tet4
partitions. It tries local edge and face swaps and bounded circumcenter
insertion where the incident cells are poorly shaped. Each change must
improve the worst local scaled Jacobian and pass the original CAD material
checks. The pass preserves original and shared-cut triangles, has separate
work and native query limits, and is off by default.

`--improve-surface-spacing` (`NativeCadMeshOptions::improveNativeSurfaceSpacing`)
tries a size-graded boundary on single solids with up to 32 face uses
after the direct native star fill declines. It samples each CAD edge once and
uses the same samples on all incident faces. Face points are placed in a
local physical metric, with chord-driven refinement near curved boundaries.
An accepted mesh reports `native_spaced_surface` and passes the regular
embedding, CAD material, volume and boundary checks using the constraint
seed route. If the bounded trial
fails, meshing continues through the established routes. The option is off
by default and can be set independently of `--improve-tet-quality`. Use
`--improve-surface-spacing --no-improve-regular-quality` to start with
this spacing route.

`--regular-surface-front` (`NativeCadMeshOptions::useRegularSurfaceFront`)
tries a graded advancing front on the same single-solid Tet4 route. It
uses the size at each edge's endpoints to predict a near-equilateral apex
in the physical face metric. A sample at the predicted apex reduces that
step where curvature or a short trim segment calls for smaller triangles.
The front leaves shared curves fixed, then repairs chord error and tries
one quality-improving sweep of interior face points. Accepted meshes
report `native_frontal_surface`. The earlier spacing mode is tried next
when the front's mesh cannot pass the existing checks. Use
`--regular-surface-front --no-improve-regular-quality` to start with this
version of the regular front.

Advanced recovery enables the bounded quality pass by default on the
eligible single-solid Tet4 fallback. The pass swaps
internal face diagonals and moves face-interior points in the physical CAD
metric. It improves the worst local triangle shape while retaining the
local count of triangles with angles between 45 and 75 degrees. CAD curve
samples remain fixed. After volume recovery, a bounded Tet4 swap pass
targets low scaled Jacobians without adding nodes or changing boundary
facets. Accepted meshes report `native_frontal_quality` and pass the same
CAD and volume checks. On refusal the regular and spaced front routes
remain available. Set `NativeCadMeshOptions::improveRegularMeshQuality` to
false or pass `--no-improve-regular-quality` to disable it. Standard
recovery uses its existing route order by default. Set the option to true
or pass `--improve-regular-quality` to request the quality route explicitly.

Adaptive local partition recovery can validate a facet sample with a point
on the original trimmed face that lies within the chord limit. If the
projection cannot establish that bound, PAACS uses the existing
trimmed-face distance query. Same-face midpoint checks can reuse a verified
witness for the same joined edge. The audit still accounts for all four
samples per triangle and reports the largest witness distance as an upper
bound on the measured sample error.

`--live-material-recovery` enables a bounded single-solid Tet4 trial in
Advanced recovery. It closes the original triangular boundary first, then
tries to repair CAD-rejected tetrahedra using local edge-ring swaps, interior
point insertion and closed interior vertex-star relocation. Each accepted
patch keeps its oriented boundary and checks replacement centroids against
the original solid. Standard CAD, boundary and volume checks still gate
publication. The prior recovery routes remain available on a refused trial.
The `live_material_*` receipt counters describe the live repair attempts.
The final publication and validation fields describe the mesh that was
actually saved.

`--local-delaunay-refinement` enables an experimental interior refinement
trial for strict Tet4 meshing. In Advanced recovery it also enables
the live material route. It can use standard, regular or
curvature-sized surfaces when the corresponding route recovers
a Tet4 fill. It uses the circumradius divided by the local size to
prioritize insertion. Initial vertex sizes come from incident edges,
and each inserted vertex interpolates the size in its containing
tetrahedron. The trial grows a cavity from the selected tetrahedron,
checks its visibility and protected boundary, then improves the
worst radius ratio through bounded edge swaps and interior relocation.
On the curvature route, it first tries placement and optimization
on the recovered discrete mesh, then repairs any original-CAD material
rejections. If the candidate cannot pass those checks, the route can
try local support insertions and the existing global refill. Original-CAD
material validation still gates publication. A positive minimum quality
requirement or an outside-centroid tolerance leaves the trial inactive.
The API option is
`NativeCadMeshOptions::useLocalDelaunayRefinement` and defaults to false.
The `local_delaunay_refinement` receipt field records the requested
option. Per-body `local_delaunay_support_nodes` and
`local_delaunay_refinement_nodes` count local insertions,
`local_delaunay_quality_edits` counts accepted swaps and relocations, and
`local_delaunay_used_global_fallback` records a return to global refill.

`--curvature-points-per-turn N` enables a bounded single-solid Tet4 trial
in Advanced recovery. The trial limits edge and surface spacing by local
curvature as well as the requested size and chord error. Adjacent faces
reuse the same sampled edge chain. If the imported face trim disagrees
with an interior edge sample, the trial can fit its parameter on the
same trim within the original edge-use tolerance. The API field is
`NativeCadMeshOptions::curvaturePointsPerTurn`. The range is 0 through 64,
with 0 disabling the trial. Add `--live-material-recovery` to try local
material repair on the curvature-sized boundary. When the initial fill
uses only boundary nodes and rejects more than 16 centroids, the combined
trial can place up to 96 original-CAD-admitted inward points beside the
rejected cells and refill the same shell once. Each point and its
corresponding face-to-point tetrahedron centroid must pass original-solid
classification.
An initial distance-to-shell query prioritizes possible repair sites.
Every centroid of the retained volume still needs an original-solid
admission. The serial retry can reuse positive decisions for bit-identical
coordinates when the original solid and tolerance match.
For a small number of remaining cells, the combined mode can also place
an interior point along the inward direction shared by two protected
facets. It keeps those facets fixed and checks every replacement cell
before committing the cavity.

`--trimmed-ray-witnesses` is an experimental option for serial, strict
validation in the live Tet4 route. It accepts a positive decision from
two opposing ray pairs against the original trimmed CAD faces.
Conservative face bounds limit the intersections. Unsupported solids,
ambiguous rays and points without a positive witness use the existing
solid classifier. The receipt records ray admissions separately. The
ordinary classifier remains the default.

`--trace-native-stages` writes flushed stage events to stderr. The trace
records started and finished scopes, including the active stage when an
external deadline terminates the process. The JSON mesh receipt is unchanged
in format apart from additive fields.

Before the live adjacency and Steiner change, a volume-only CAD23 trial at a
450-million-unit stage allowance reduced the missing edge count from 1,202 to
223 and the missing triangle count from 1,662 to 255 before carving. Its
22,794 retained tetrahedra cover 95.0% of the enclosed volume. The independent
discrete audit found zero nonpositive cells and zero illegal tetrahedron
pairs. There are still 291 source triangles absent from the retained
material boundary and 453 additional front triangles. A 150-second refill
trial reached its deadline. A smaller work cap was exhausted during
preparation before either cavity search began. The original STEP still
refuses at the stock surface edge check. See the [shared section through the
actual cells](../evidence/local-2026-11-08/recovery-mode-checks/CAD23-constraint-recovery-section.png).

The live recovery replay on that same surface reached 96.3% enclosed volume
in 66 seconds. It staged 54 Steiner points and retained 40 in the carved
seed. The uncarved complex had 155 missing edges and 177 missing triangles.
After carving, 209 source triangles were still absent from the material
boundary.

An eight-file ABC STEP holdout excluded families in the original 51 and
previously manifested inputs. With partitioning disabled to exercise native
volume recovery, the final build published two advanced meshes and zero
standard meshes under the same 90-second deadline. CAD237 recovered its
entire constraint set. CAD393 used the partial seed and the existing
refill to complete. Both archives passed native validation and independent
exact orientation and contact audits, with quality warnings. Five files
stopped at the assembly interface check in both modes, and CAD205 reached
the deadline in both. See the frozen sample (`evidence/local-2026-11-08/live-recovery-holdout/manifest.json`)
and matched results (`evidence/local-2026-11-08/live-recovery-holdout-final/results.json`).

On the original CAD06 STEP with partitioning disabled, the advanced seed
and native refill produced 546 tetrahedra and a 476-triangle boundary.
The result passed the native checks and an independent exact orientation
and contact audit. Its CAD volume error was 0.120% and its minimum
mean-ratio quality was 0.013. CAD06 already passed the qualified corpus
through a different route. See the native receipt (`evidence/local-2026-11-08/recovery-mode-checks/cad06-advanced-direct-volume/receipt.json`).

On CAD27's archived qualification profile, both modes reached a 260-second
deadline on the current executable. The isolated advanced stage produced
an 87.4%-filled private seed in 11.5 seconds. The current end-to-end
timings are unresolved for this profile.

### Local facet and cavity recovery (9 November 2026)

Advanced recovery can walk a missing facet from its recovered edge stars and
replace the complete star of a tetrahedral edge piercing the facet. The
replacement must recover the facet or reduce its conflicts while preserving
the protected features and the oriented cavity boundary.

After carving, bounded local rules close small cavities before invoking
broader region search. One Tet4 can consume three or four open triangles.
Two Tet4 can advance together when the second reduces the front. The solver
validates both cells and their contact before insertion.
Spatial indices restrict exact contact tests to nearby candidates. The
native mode retains the existing publication checks and the same CLI toggle.
Core callers can opt in through `TetPlcCavityOptions::maximumLocalFrontTets`.

On the fixed CAD23 surface, piercing-edge recovery applied 19 patches and
reduced the missing facets before carving from 177 to 174. Local rules
reduced the carved front from 514 to 396 triangles in 0.34 seconds after
preparation. The combined refill reached 340 open triangles in 111 seconds.
Its 23,213 retained Tet4 passed independent exact orientation and pair
checks. The checkpoint still lacked 132 source triangles and no complete
mesh was published. The original STEP rerun refused at `nonconforming_edge`
before volume recovery. See the residual incidence (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-hybrid-rules/partial-incidence.json`).

The previous and updated advanced builds both published five of eight
meshes on a fresh, size-stratified ABC STEP sample prefiltered for single
solids. The other three reached the 90-second limit in both builds. Native
partitioning was disabled to exercise volume recovery. The local rules
contributed cells to three published meshes. Each passed independent exact
boundary incidence, orientation and pair checks. On CAD643, one matched run
reduced refill preparation from 27.94 to 1.95 seconds and total time from
37.45 to 10.55 seconds. See the frozen sample (`evidence/local-2026-11-08/hybrid-recovery-holdout/manifest.json`),
previous results (`evidence/local-2026-11-08/hybrid-recovery-baseline/results.json`)
and updated results (`evidence/local-2026-11-08/hybrid-recovery-final/results.json`).
CAD237 and CAD393 also passed again. CAD237 reproduced its archive byte for
byte. CAD393 used 41 local-rule tetrahedra and its updated mesh passed the
independent exact checks. See the replay results (`evidence/local-2026-11-08/hybrid-recovery-prior-passes/results.json`).

### Repairing native face charts from shared curves (9 November 2026)

The advanced native fallback can rebuild a damaged OCCT face triangulation
from the shared CAD curve samples. It detects conflicting vertex bindings,
invalid interior parameters, and missing wire segments. The rebuild uses
each wire occurrence's own parameters and constrains those segments in the
face triangulation. Distinct CAD vertices keep their original identities.
A small UV detour can collapse only at the same degenerate CAD vertex and
within its existing tolerances.

The advanced fallback keeps the physical chord, closed-shell, embedding,
and source immutability checks. It can use a directly evaluated point on the
trimmed face as an upper bound for the distance from a facet sample. When
that witness fails, the existing native projection still runs. The reported
maximum facet deviation on this route is an upper bound at the checked
samples. An accepted midpoint witness can be reused for a neighboring
triangle on the same CAD face. Rebuilds are capped at 32 charts and refuse
faces shared between different solids.

From the original CAD23 STEP, the production adapter rebuilt nine charts
and produced a closed, embedded surface with 7,877 nodes, 15,826 triangles,
all 454 CAD face tags, and 0.0551% surface-volume error. Its archive matched
the prior experimental surface byte for byte. One staged comparison reduced
surface time from 96.92 to 69.38 seconds and global projections from 6,417
to 379, excluding import and exact embedding. The whole-file advanced runs
reached their 300-second deadlines without publishing a volume mesh. See
the production surface figure (`evidence/local-2026-11-09/CAD23-production-surface.png`),
staged surface receipt (`evidence/local-2026-10-18/recovery-followup/cad23-production-wire-chart-midpoints/stdout.log`)
and whole-file receipt (`evidence/local-2026-11-09/cad23-wire-chart-cached-original/receipt.json`).

### Faster exact checks during advanced volume recovery (9 November 2026)

The advanced route ranks large edge-star replacements before testing the
selected replacement against the exact cavity shell. If the selected patch
fails, the existing search checks each candidate against the shell. A strict
segment crossing through a triangular base certifies the two equivalent
triangulations of a convex bipyramid for a 2-to-3 or 3-to-2 flip. This allows
the solver to reuse the old patch's exterior contact result. General patches
still check the new cells against each other and the retained mesh. Steiner
trials use bounding boxes and exact supporting planes to avoid unnecessary
triangle contact tests.

The live spatial index adds new cell IDs in batches and consolidates when
enough batches accrue. An interrupted build leaves the previous batch index
usable or invalidates an incomplete full rebuild before a retry. Material
labeling walks the live tet neighbours across unprotected faces.

On the frozen CAD23 surface, the matched control spent 13.99 seconds in large
stars, 14.15 seconds checking patch contacts, and 16.88 seconds in Steiner
trials. Both a post-Steiner sweep and a bounded local flip-sequence trial
ended with the same 155 missing edges and 174 missing facets as the control.
See the control profile (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-profile-control/stdout.log`)
and local flip-sequence trial (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-local-lookahead/stdout.log`).

At the same 450-million work cap, the final build reduced the CAD23
constraint stage from 57.67 to 21.53 seconds and from 238.56 to 85.82 million
counted work units. The 23,107-tet checkpoint is byte-identical to the
control. Its independent audit found no nonpositive tetrahedra or illegal
contacts. The checkpoint covers 96.32% of the volume enclosed by the source
triangles. Before carving, 174 constraints remain missing.
See the optimized profile (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-optimized-450000000/stdout.log`)
and exact checkpoint audit (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-spatial-batches/exact-embedding.log`).

At a 100-million work cap, the previous build exhausted its budget before
exporting a checkpoint. The updated build retained 22,765 tetrahedra covering
91.07% of the discrete volume in 17.48 seconds. These profiles exclude STEP
import, CAD surface generation and native volume acceptance. See the
smaller-cap control (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-profile-control-100m/stdout.log`)
smaller-cap result (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-optimized-100000000/stdout.log`),
and exact audit (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-optimized-100000000/exact-embedding.log`).

The matched native holdout reproduced all five published archives byte for
byte. The other three again reached the 90-second limit. CAD525 went from
3.85 to 1.44 seconds and CAD643 from 10.72 to 7.11 seconds in these runs. See
the native results (`evidence/local-2026-11-09/recovery-performance-holdout/results.json`).
Eight focused tests passed. The Steiner test exercised its recovery witnesses
and the first 32 allocation-failure positions.

### Fixed-surface Tet4 follow-up

The live recovery solver now checks the vertex table at construction and
reuses that invariant during facet-contact queries. Each query still checks
the local geometry. For a positive tet, the separating-plane test uses the
known inward sign. A 2-to-3 flip with an exact convex bipyramid certificate
computes its obstruction change from the four affected faces directly.

On the fixed CAD23 surface at the 450-million work cap, the constraint stage
went from 21.06 to 14.44–14.99 seconds. The 23,107-Tet4 checkpoint has the same SHA-256
as the independently audited prior result. It still has 174 missing source
facets before carving and covers 96.32% of the discrete volume. These timings
exclude STEP tessellation and the separate embedding audit. See the
PAACS Release-linked profile (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-local-exact-linked/stdout.log`).
At the 100-million cap, the byte-identical 22,765-tet checkpoint went from
15.80 to 10.73 seconds. See the
smaller-cap profile (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-local-exact-100m/stdout.log`).

Two narrow strips account for most of the residual gaps. Their perimeters are
already recovered and their missing internal diagonals reach 138 model units.
The gap plot (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-local-exact-final/cad23-boundary-gaps.png`)
shows the 203 missing source facets after carving, including 127 in those strips.
The directed star trial (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-directed-star/stdout.log`)
staged flips around the rim of a crossed face before committing a patch. It
left 173 facets missing and reduced enclosed coverage to 96.14%. Recovering
long segments before protecting any facets left 202 missing. An implicit
exterior trial also failed to close the boundary. The experimental schedules
and implicit exterior were kept out of the live path. See the
segment-first trial (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-directed-segments-first/stdout.log`)
and implicit exterior trial (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-implicit-exterior/stdout.log`).
The older two-apex replacement could use none of the 790 walked cavities
tested. Its exact shell precondition requires a single rim that connects both
endpoints of the missing edge. See the
two-apex trial (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-two-apex-rim/rim.jsonl`).

The fresh native build passed all six focused CTest targets. On the same
eight STEP files, all five previously accepted archives reproduced byte for
byte. The remaining three again reached the 90-second deadline. See the
native holdout (`evidence/local-2026-11-09/local-exact-native/holdout/results.json`),
native tests (`evidence/local-2026-11-09/local-exact-native/native-focused-ctest.log`),
and core tests (`evidence/local-2026-11-09/local-exact-native/core-focused-ctest.log`).

### Constraint recovery after a stalled flip pass

The advanced scheduler runs the long-link and existing-node cavity pass
when the short-link pass stalls. A six-tet regression requires this
escalation and verifies that the replacement preserves the oriented cavity
shell. The fresh native holdout reproduced all five accepted archives byte
for byte. Three files again reached the 90-second deadline. The current
CAD23 fixed-surface run still leaves 174 constraints missing before carving.
See the regression (`evidence/local-2026-11-10/constraint-escalation-native/front-test-ctest.log`),
holdout (`evidence/local-2026-11-10/constraint-escalation-native/holdout-comparison.json`),
and CAD23 profile (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-escalation-production/stdout.log`).

The CAD23 cavity probe found a six-triangle strip that splits an 82-tet
cavity into two closed sides. Delaunay recovery and the local cavity audit
completed its material side in 13 Tet4. Merging that side into the carved
checkpoint reduced missing source facets from 203 to 194 and extra front
faces from 311 to 300. Independent checks found no inverted cells, illegal
tetrahedron contacts, patch intersections with the original surface, or
reversed source facets. See the merge inventory (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-strip-2-material-merge/merge.json`),
embedding audit (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-strip-2-material-merge/exact-embedding.log`),
and patch-to-source audit (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-strip-2-source-audit/stdout.log`).

The two 62-triangle strips need more space for recovery. The bounded
expansion produced one manifold material shell. The other still has a
disconnected vertex link. Grouping several facets around common support
vertices reduced the first shell's open front from 158 to 64 faces before
the refill stalled. The exact shell audit found zero illegal triangle pairs.
See the cavity inventory (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-connected-strip-cavity-alternate/strip-cavities.jsonl`)
and grouped refill (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-long-strip-0-guided-refill/stdout.log`).

### Complete CAD23 volume on the frozen triangular surface

In the long-edge probe, a walk with only convex 2-to-3 moves completed
none of nine tested missing source edges. Allowing 3-to-2 moves along the
whole path completed one after 28 moves.

Advanced PAACS recovery now tests adjacent cells when the crossed cavity
has no legal replacement. An exact check for a possible cyclic edge link
rejects impossible shells before searching tetrahedra. This reduced missing
CAD23 facets before carving from 174 to 167. The native holdout retained all
five passes. Four archives were identical, and the changed archive passed
the independent embedding audit. See the
fixed-surface comparison (`evidence/local-2026-11-10/cad23-hull-recovery/comparison.json`)
and native holdout (`evidence/local-2026-11-10/residual-cavity-native/holdout-comparison.json`).

The staged material repair completed CAD23 with **24,218 Tet4** and all
**15,826 original oriented source triangles**. It filled both long strips in
groups, then reused the updated open front to identify independently
fillable chambers. Reopening eight adjacent retained tets and prioritizing
the remaining source triangles for interior point placement closed the
last surface gap. A final 201-tet replacement closed the internal void.
All tetrahedra came from PAACS. The exact audit found zero nonpositive cells
and zero illegal contacts. The independent boundary and volume certificate
found zero missing or extra facets and relative volume error of 1.08e-15.
See the complete archive and certificate (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-complete-volume-certificate`),
intersection audit (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-material-strip-complete/exact-embedding.log`),
and replay driver (`evidence/local-2026-10-18/recovery-followup/replay_cad23_material_strip.py`).

The existing hybrid volume refill independently closed the 50-face interior
front after the staged source repair. Its 24,033-tet result also passed exact
embedding checks. See the hybrid closure (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-source-closed-hybrid/stdout.log`)
and independent audit (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-source-closed-hybrid/exact-embedding.log`).

This is a frozen-surface result from the staged experiment. The native STEP
spot check classified 14 of 32 sampled interior nodes and two of 32 sampled
tetrahedron centroids outside the original solid. The larger measured
centroid distance was 0.01683 mm. The minimum mean-ratio quality was 1.16e-11
and the median was 0.187. See the
STEP sample (`evidence/local-2026-10-18/recovery-followup/constraint-cad23-complete-native-sample/stdout.log`).

## Historical CAD qualification and follow-up checks

The historical qualification passed **40/51**. Follow-up checks add
whole-file evidence for CAD01, CAD10, CAD14 and CAD28. **44 distinct files**
have a recorded pass across builds and profiles. That tally includes
independently validated assemblies that preserve each part and report
interference.

An early fixed-surface CAD23 candidate filled 81.4% of the enclosed
discrete volume. Its remaining front had 467 non-source triangles in
four connected groups. The checkpoint lacked 173 required edges and
195 required source facets. Later staged results are recorded above.

The next volume prototype should keep the full Delaunay candidate while
recovering missing source edges and facets. Use indexed connected-cell patches
and transactional replacements. Retain cells around unresolved constraints
until those constraints are recovered. After the recovered facets form barriers,
label material, exterior and cavity regions through cell adjacency. Then batch
original-CAD checks with isolated classifiers while retaining strict interior
node inclusion and the explicit centroid policy. A passing frozen CAD23 test
must retain every source facet, close the volume, and pass exact embedding and
original-CAD checks before the experimental route can be promoted.

## Recovery dispatch and surface reconstruction (20–21 October 2026)

The main mesher now tries the established partition policies before the two
newer adaptive surface-recovery policies. A traced CAD9999 run showed the newer
ordering spending its 600-second allowance inside surface recovery. Restoring
the established priority completed the same command in 468 seconds and produced
the identical historical archive. Geometry, topology and resource gates are
unchanged. These observations were not isolated performance benchmarks.

The rebuilt executable passed **all 89 tests**. Five original files were also
replayed and independently audited, covering eight original bodies:

| File | Original bodies | Stored Tet4 cells | Centroid policy |
| --- | ---: | ---: | --- |
| CAD9999 | 1/1 | 34,880 | Strict native IN/ON |
| CAD8181 | 1/1 | 51,022 | Strict native IN/ON |
| CAD24 | 1/1 | 20,907 | Strict native IN/ON |
| CAD01 | 2/2 | 15,667 | Explicit recorded CAD chord allowance |
| CAD28 | 3/3 | 14,716 | Explicit recorded CAD chord allowance |

The first three archives are byte-identical to their earlier accepted outputs
under the same flags. The last two use independent-parts assembly semantics:
each original part is retained and audited separately, with interference
reported. They do not establish a globally nonoverlapping assembly domain.
All five have independent archive, boundary, exact cell-contact and original-CAD
checks. The serial CAD24 native audit reached its 900-second deadline; its
complete checks then passed in eight disjoint index partitions, using isolated
native processes. That auditor reproduced every aggregate metric from the
serial CAD9999 control. No nodes or cells were sampled away.

Element quality remains a substantial limitation: CAD24's minimum reported
mean-ratio quality is approximately `5.4e-21`, despite passing exact orientation
and contact checks. These are geometry/validity preservation examples; they
are not simulation-ready quality claims or a fresh 51-file qualification.

CAD23 now has a **certified experimental surface**, outside the production
adapter. Native face 401 reused a tessellation node for two distinct original
CAD vertices separated by `0.493110 mm`. The prototype reconstructs complete
face wires with constrained triangulation and preserves separate seam/pole
parameter occurrences. Tiny parameter detours are contracted only at the same
degenerate native vertex, with bounded parameter extent and source-tolerance
physical witnesses. No distinct 3D vertices are merged.

Automatic failed-face selection reconstructs nine faces and reproduces the
same surface archive without a face-ID list. It retains all **7,877 original
3D node coordinates/classifications byte for byte** and all **454 CAD face
tags**, producing **15,826 triangles**. Closed orientation and embedding,
original native facet checks and source immutability pass. Maximum sampled
facet deviation is `0.216350 mm` against a `0.583153 mm` allowance; enclosed
surface-volume error is `0.055135%`. The same bounded pole treatment still
refuses CAD13's crossing trim curves. CAD23 volume recovery remains experimental
and has not earned an original-file Tet4 success.

The execution ledger (`evidence/local-2026-10-18/recovery-followup/surface-dispatch-status-2026-10-20.json`)
pins commands, builds, inputs, policies, outcomes and remaining work. See the
five-file cutaways (`evidence/local-2026-10-18/recovery-followup/dispatch-preservation-gallery-2026-10-20.png`)
and CAD23 surface repair (`evidence/local-2026-10-18/recovery-followup/CAD23-certified-surface-repair-2026-10-21.png`).
The incomplete CAD23 volume checkpoint (`evidence/local-2026-10-18/recovery-followup/CAD23-volume-recovery-checkpoint-2026-10-21.png`)
shows the remaining material and source-facet coverage gaps after 118 local
cavity repairs; it contains only `27.883%` of the enclosed surface volume and
must not be counted as a completed Tet4 mesh.

Subsequent private bulk-construction trials exposed another concrete gap:
source-contact checks alone can retain cells in empty concavities. The raw
bulk checkpoint contains `197.3%` of the source volume. A construction filter
using exact segment/triangle predicates against the certified closed source
rejects outside proposals before costly native queries. It passes analytical
checks for a cavity, disconnected bodies, thin material, boundary ambiguity
and a concave L-shaped solid. This is separate from the production test suite.
The actual before/after cells (`evidence/local-2026-10-18/recovery-followup/CAD23-bulk-containment-comparison-2026-10-21.png`)
show the empty holes preserved by that filter.

The filtered bulk checkpoint contains `54.896%` of the enclosed surface volume;
after 264 connected-region repairs its retained candidate reaches `69.634%`.
Including the admitted interior seeds in the bulk Delaunay triangulation also
reduces the open front after growth from 3,170 to 1,926 triangles. Its subsequent
175 local repairs retain 24,001 candidate tetrahedra and `81.377%` of the
enclosed surface volume. The two bounded trials reached their respective
1,200-second and 900-second deadlines during original-CAD material filtering.
Neither published an accepted volume. A sparse per-query profile completed
64 original-native checks in 7.81 seconds, with one classification taking
3.02 seconds; checked batching is a concrete next performance target.
These remain incomplete diagnostics. Every final original-CAD and volume gate
remains required. The discrete centroid diagnostic separately records three
classification changes caused by summation-order displacements no larger than
`1.44e-14 mm`; these must not be confused with the raw bulk's large unwanted
volume. Final acceptance uses the recorded CAD tolerance policy.

The successful native build used the C++23 compatibility configuration;
C++26 qualification remains pending with this toolchain.

With the recorded OCCT runtime paths configured, the suite can be replayed as:

```powershell
& 'C:/Program Files/CMake/bin/ctest.exe' --test-dir build/windows-msvc-release `
  -C Release --output-on-failure --timeout 1800 -j 3 --no-tests=error
```

`importOcctCad` in `geometry/OcctModel.hpp` retains valid bounded solid BReps,
including curved surfaces and disconnected located bodies. It exposes native
geometry evaluation, projection, inside/outside classification, entity
tolerances, oriented solid-face uses, wire-edge uses and use-specific pcurves.
A cylinder seam has two distinct pcurve uses of one shared edge. STEP
coordinates are transferred explicitly to millimetres; BREP units are
unspecified. Whole inverted solids representing unbounded material are refused.

Tags and `identity()` are stable for the lifetime of one imported snapshot.
They do not establish persistent correspondence across edited STEP files.
Located bodies preserve geometry and topology; full XDE product names,
occurrence provenance and persistent selector archives remain future work.

## Explicit CAD tolerance for PLC volumes (19 October 2026)

Strict native centroid inclusion remains the default. A linear triangulation
of a concave curved face can place otherwise conforming tetrahedra slightly
outside the exact solid. `--allow-centroids-within-chord-tolerance` enables an
explicit approximation policy for the `native_plc_cavity` route:

- Interior nodes must still classify IN/ON at the original native tolerance.
- An OUT centroid must have a measured distance to the original solid no
  greater than that body's configured chord allowance. UNKNOWN refuses.
- Four additional interior samples per tetrahedron must satisfy the same bound.
- Exact positive orientation, source-boundary conformity, embedding, contact,
  original entity identity, sampled facet error and volume checks remain gates.

The final independent-parts assembly rechecks centroids and the four interior
samples against each original solid. It enables the band only for a child PLC
result carrying the complete surface and sample certificates. Other routes
retain their existing containment checks. No parts are fused or coordinate-welded.

```sh
paacs-cad-mesh input.step 1 output.paacsmesh --native-cad --size-factor .15 \
  --independent-parts --allow-centroids-within-chord-tolerance --report-json
```

`--max-chord-error` sets the geometric allowance explicitly in model units;
STEP import uses millimetres. The JSON report records the enabled policy,
per-body limit, OUT centroids accepted within the band, actual maximum distance,
additional sample counts and distances, and original-source guards. These are
sampled CAD fidelity checks, rather than a continuous Hausdorff certificate.
Results using this option must be reported separately from strict-inclusion
corpus results. Element quality remains a separate warning or requested gate.

The fresh original CAD01 run now completes both parts: **12,523 Tet4 cells and
3,858 nodes**. Part 1 uses the PLC route and part 2 the established partition
route. All eight independent part checks pass, including exact cell contacts
and original-CAD fidelity. There are no inverted cells, improper contacts,
missing boundaries, outside interior nodes or samples beyond the chord bound.
Part 1 has 125 native-OUT centroids; their maximum distance is `0.0393174 mm`
against its `0.356678 mm` allowance. Part 2 has zero OUT centroids. Relative
volume errors are `0.4434%` and `0.2096%`. Quality warnings remain substantial.

The 88-test suite passes, followed by four freshly rebuilt native/assembly
checks after final policy propagation. CAD10, CAD27 and CAD28 still complete
under strict centroid inclusion; CAD28 is byte-identical, and CAD10 has an
exact coordinate/classification and cell-geometry bijection with its previously
independently audited output. These preservation checks do not replace a new
full 51-file qualification with the optional distance policy.

See the CAD01 execution receipt (`evidence/local-2026-10-18/recovery-followup/CAD01-next-file-receipt.json`)
and actual mesh cutaways (`evidence/local-2026-10-18/recovery-followup/CAD01-complete-CAD-tolerance.png`).

An experimental saved-surface replay also completes original CAD14 part 7
with 7,441 Tet4 cells. Its independent exact contacts and original-CAD audit
pass with zero OUT centroids and no interior nodes outside native material.
It combines PAACS Delaunay candidates, CAD-audited curved boundary diagonal
changes and coordinated native-wire refinement. This is a part-level trial;
it is not integrated as a production route or counted as another complete file.
The part 7 cutaway and transverse section (`evidence/local-2026-10-18/recovery-followup/CAD14-part7-completed-recovery-trial.png`)
show the preserved hollow region in the exported tetrahedra.

## Native meshing

`meshNativeCad` in `geometry/OcctMesher.hpp` first recognizes supported bodies
against the authoritative BRep. Whole analytic parts must match a newly
constructed OCCT body by volume and Boolean common volume before the existing
primitive engine is selected. Primitive classifications are then resolved to
the actual imported CAD entities. Planar bodies use the current geometry/core
pipeline. Separated bodies are assembled without merging coincident nodes.

| Input | Native route |
|---|---|
| Supported planar blocks, sweeps and polyhedra | Existing mapped/swept/constrained core strategies |
| Complete cylinder or coaxial tube | Existing cylinder/tube engine |
| Complete cone or frustum | Existing graded cone/frustum engine |
| Full solid torus or planar-capped solid torus segment | Existing uniform torus/elbow engine |
| Rotated/translated supported part | Same engine in the imported coordinate frame |
| Geometrically separated supported bodies | Independent body meshes with imported body identity |
| Sphere, partial cylinder and other bodies with a certified star-shaped native boundary | Native shared-edge surface triangulation and certified Tet4 fill |
| Hollow/concave bodies whose bounded native partition pieces are all certified | Native copied-body partition, shared-interface assembly and Tet4 fill |
| Unresolved shared interfaces, overlapping bodies, or unsuccessful bounded recovery | Explicit refusal with stage diagnostics |

Tet4 is the native CLI default. `--elements native` selects the existing
linear mixed/structured routes. Imported curved Tet10 is not offered: existing
primitive Tet10 has straight midsides, which would not be an exact CAD boundary.
Uniform sizing remains the default on the supported uniform primitive routes.
Cones, tubes and intentionally varying size controls have different measured
cell distributions; the visual report records their actual quality and volumes.

The native fallback copies and tessellates CAD geometry. Boundary nodes join by
native vertex/edge identity and sample ordinal; nearby coordinates alone never
establish a join. Face trims, periodic seams, poles and native orientations are
retained. Exact triangle contact predicates certify an embedded closed boundary.
A central fill requires a strictly certified kernel inside every triangle's
outward halfspace; a guessed center is insufficient.

The partition fallback splits a copied native body with bounded cutting planes.
It conserves native volume and maps derived entities back through OCCT history.
All pieces use one surface session. Shared cut triangles must have identical
native node keys, bitwise matching coordinates/classes and opposite winding
before they cancel. Original exterior faces and body identities remain in the
assembled mesh. Final checks cover exact exterior embedding, complete volume
incidence, positive cells, original CAD node projections, tetrahedron centroid
inclusion, sampled boundary error and volume agreement. They do not establish
continuous CAD containment or a Hausdorff certificate.

These routes prioritize valid coverage. They can produce very slender cells;
global size controls surface tessellation and does not certify volume edge
spacing. `--warning-quality .25` reports low quality without rejecting an
otherwise accepted mesh. `--min-quality .25` remains a strict publication gate.
Selected-face/spatial-field controls and mixed/high-order elements on these
fallbacks are explicitly unsupported. `--no-native-fallback` restores the
earlier planar/analytic routing; `--no-native-partition` disables partitioning
while keeping the direct native surface route. `--max-native-pieces` bounds
partition complexity (default: 512 pieces). Per-attempt partition diagnostics are kept.

`--separate-contact-bodies` explicitly retains separate nodes at point/edge
contacts after native material/face overlap checks and exact cross-body surface
contact checks. Shared face contacts require an explicitly supported interface policy.

`--separate-matched-box-faces` adds a narrow full-face contact route for native
axis-aligned boxes at eight-corner Tet4 resolution. The adapter proves each
actual solid against its native expected box, evaluates endpoints from their
original CAD vertex identities, and requires unique, exactly coincident,
oppositely wound triangles covering the complete native contact area. Body
nodes and IDs remain separate. Partial faces, material overlap, mismatched
triangulations and finer unsupported box resolutions refuse transactionally.
This is a separate contact boundary model; it does not construct a bonded
multimaterial interface. JSON receipts record matched face pairs/triangles and
the maximum coordinate correction from native vertex evaluation.

The star route also has bounded numerical halfspace proposals after its normal
kernel search fails. Each proposal must pass every existing exact facet-side
test before construction; numerical feasibility alone never admits a mesh.
All proposal attempts consume the shared work allowance.
Near-parallel constraints have an additional bounded numerical proposal retry;
it passes the same exact tests. No proposal changes boundary coordinates or
relaxes the acceptance predicates. Feature planes come from actual native
cylinder/cone and circle geometry, and each copied partition must conserve the
source volume before meshing. Unsupported partition attempts can proceed to
the next bounded policy while retaining diagnostics and cumulative work.

Each mesh passes final incidence and quality audits before publication. Reports
include BRep/mesh volume difference and boundary-node/face-centroid projection
errors. These samples and the primitive chord bound do not establish a general
Hausdorff certificate. Failure, cancellation and exhausted budgets preserve the
caller's previous mesh; the archive writer preserves an existing output file.
The native CLI rejects input/output aliases, including existing hard links,
before import so a mesh export cannot replace its source CAD file.

```powershell
build/windows-msvc-release/Release/paacs-cad-mesh.exe `
  tests/fixtures/cad/cylinder.step .3 out/cylinder.paacsmesh `
  --native-cad --report-json
```

Read `bodies[].start_face_tag` and `end_face_tag` from that receipt before
selecting an imported cylinder/tube cap. For the generated unit-cylinder
fixture, start face 3 has been tested with:

```powershell
build/windows-msvc-release/Release/paacs-cad-mesh.exe `
  out/visual-validation-2026-10-03/native/cylinder.step .3 out/fine-cap.paacsmesh `
  --native-cad --face-size 3 .15 --min-size .12 --growth 1.3 --report-json
```

Cylinder cap sizing varies the interior section and axial spacing. Tube cap
sizing controls axial spacing with the existing swept section.

### Native size controls

For native sizing on one imported Tet4 solid, select `--native-sizing` and
pass `--face-size TAG SIZE` for a tagged CAD face. Add one or more
`--size-box X0 Y0 Z0 X1 Y1 Z1 SIZE` controls in model coordinates. A box
also selects the native sizing route. Use `paacs-native-cad-preview` to obtain
face tags from the same STEP file. STEP import coordinates are in millimetres.

The field combines the global size, the distance to each selected trimmed
face and the distance to each box. The smallest active request controls shared
CAD curves, trimmed-face charts and interior tetrahedron insertion.
`--min-size` clamps requested field values before the growth ramp. CAD
geometry and boundary recovery can still require shorter edges.
`--growth G` controls the field slope, `2*(G-1)/(G+1)`. Achieved lengths and
neighbor ratios depend on the boundary, recovery and refinement work limits.
The native topology, orientation and CAD acceptance checks also apply to this
route. Assembly sizing and size-controlled native hex cells are currently
unsupported.

`scripts/native_step_sizing.py` saves a portable recipe with its STEP input
and input hash. Use its `mesh` command to reopen the recipe, apply its controls
in `paacs-cad-mesh`, and write the native archive and a tagged C3D4 solver
include. The script checks all exported body and face groups against the
archive and reports selected-face edge lengths, box-interior Tet4 sizes and
adjacent Tet4 size ratios. See `save --help`, `mesh --help` and `plot --help`
for the arguments.

`--optimize-profile` runs four bounded iterations of the existing interior-node
quality optimizer on a profiled-volume cylinder. Boundary coordinates and
connectivity remain unchanged. CAD containment/projection, volume, incidence
and quality are checked afterward. `--warning-quality .25` reports a warning
for a low final floor; `--min-quality .25` refuses publication if that strict
floor cannot be achieved. The measured example improves from 0.115 to 0.167,
so a requested 0.25 floor is still refused. Nominal profile-section/axial counts
describe construction; interior positions can change during optimization.

## Recovery and bounds

`--planar-recovery --general-recovery` explicitly enable the bounded planar
surface route and existing-node general edge-cavity extension. Defaults
preserve previous core routing. The cavity solver has limits on cells,
vertices, enumerated tetrahedra, search states and exact work, plus a cumulative
scheduler budget. Bounded adjacent expansion preserves exact cavity boundaries
and protected identities; it does not add Steiner nodes or relax constraints.
This extension does not promise that every polyhedron can be recovered.

The native CAD core can also prepare a bounded replacement cavity when its
boundary has an edge pinch or a disconnected vertex link. It adds incident
cells from the same temporary region, without crossing a protected face. A
closed boundary is only a proposal: the exact cavity-ball, replacement,
protected-entity and transaction audits still have to pass. This preparation
does not carve material or publish the initial unconstrained tetrahedral fill.
Only constraints incident to source cells being replaced are required in the
local refill; constraints owned entirely by untouched cells retain their
existing identities. CAD search selects a tetrahedron containing the required
edge first, with unchanged final acceptance checks.

The captured CAD 27 cavity now changes from 22 to 32 source tetrahedra with
the same 21 vertices and passes the exact ball audit. That patch still has no
admissible tetrahedron containing the missing edge, so full CAD 27 recovery
remains refused. See the actual cavity comparison (`out/visual-validation-2026-10-11/recovery-fix/cavity-repair-and-gap.png`)
and completed analytic empty-gap regression (`out/visual-validation-2026-10-11/recovery-fix/analytic-empty-gap.png`).
The latter is a separate U-shaped fixture, not an additional corpus pass.

Isolated recovery-point diagnostics now get past CAD 27's original failed edge.
A constructive insertion subdivides crossed tetrahedra/faces/edge stars and
uses exact-checked flips or small complete-star refills. Twenty added points
form 21 short boundary edges and 42 adjoining surface triangles; local patch
boundaries and signed volumes are preserved, and the whole temporary candidate
validates. A separate native C++ midpoint trial passes the original obstruction
with one added point. New points are classified on the native wire or face,
projected within the existing chord allowance, and checked at the original
entity tolerance. Both adjacent faces share the subdivisions. Inverted or
collapsed projected connectivity is refused; a fresh construction can be tried
under a bounded allowance.

These experiments are not integrated into the normal CLI and have not produced
a complete CAD 27 material mesh. See the
first edge insertion image (`out/visual-validation-2026-10-11/recovery-fix/segment-chain-progress.png`)
and 78 checked shared-point insertions (`out/visual-validation-2026-10-11/recovery-fix/adaptive-recovery-progress.png`)
and iteration report (`evidence/local-2026-10-11/recovery-fix/iteration-report.html`).
The bounded C++ midpoint run exhausts its phase work allowance after 78 checked
insertions. Batch subdivision increases missing edge counts; exact planar
diagonal flips reduce the initial count but do not prevent this growth. Native
interior seed tests and a whole-stage edge-before-facet trial do not complete
recovery. These outcomes receive no additional corpus credit.
The completed production replay preserves all 39/51 successes and all 77 tests
pass. Thirty-eight accepted archives are byte-identical; CAD 24 changes by two
tetrahedra with unchanged volume and reported geometry error, and independently
passes CLI validation and auditing. The core-only build without OCCT also
passes its box validation/audit. This is C++23 compatibility with MSVC
`/std:c++latest`, not C++26 qualification.

`--native-kernel-admission` also enables an appended recovery route for copied
native partition pieces that have no admissible star kernel. It constructs a
Tet4 mesh against the exact piece boundary. If a curved-face diagonal causes a
recovered tetrahedron to leave the native material, bounded attempts may change
that artificial diagonal. CAD curve segments, body identity and shared cut
triangles remain fixed. The complete original-body containment, projection,
facet-sample, volume and mesh audits still precede publication. Successful
results using this construction report `native_partition_core`.

Recovery reserves the general-cavity solve allowance and collects retained core
construction/import counters after each phase. Those retained APIs check
cancellation between phases; these counters do not bound all CPU work or imply
in-phase cancellation. Failed and discarded attempts retain their work charges.

The native CLI accepts file, node, cell, logical output, audit-storage,
validation/work and primitive auxiliary-storage bounds. `--max-working-bytes`
sets the primitive/native auxiliary payload admission bound. A very fine cap can exceed its
conservative slab bound and return `budget_exceeded`; the limit is retained.
Analytic construction uses node/cell admission limits and a construction floor,
with separately bounded profile work. OCCT parsing/Boolean/projection work and
allocations are outside these counters. No process-RSS cap is implied.

The corpus runner launches separate bounded-time processes and measures peak
working set. Import, preparation, meshing, CAD validation and final audits have
separate timing fields. The meshing field measures the whole body-engine call,
including its internal checks; native validation/audit fields measure additional
wrapper checks. On planar routes this also includes surface/volume preparation;
it is not a separately instrumented pair or a pure volume-generation benchmark.

## Bounded PLC cavity recovery and native integration (18 October 2026)

The portable core now exposes `bnd::recoverTetPlcCavities` in
`paacs/bnd/TetPlcCavityRecovery.hpp`. It completes cavities in a positive Tet4
partial mesh against an immutable, closed Tri3 PLC. It is independent of OCCT.
The source node prefix and its classifications must match exactly; additional
input points and volume cells must belong to the supplied volume region.
Source face/node classifications survive publication. Unused additional points
are compacted, without coordinate-only welding.

For a stuck cavity, the method proposes an interior point and enlarges the
cavity across faces that prevent that point from seeing the entire boundary.
Original PLC triangles are hard barriers. The enlarged cavity is replaced by a
star patch only after exact positive orientation, original-boundary contacts,
all proposed/retained cell contacts, and the complete oriented front ledger
pass. This is a local, verified reconstruction, rather than a fan of the entire
CAD body. Work, storage admission, patch size, candidate count, cancellation,
and allocation failures are explicit; refusal preserves both inputs and the
destination. Low quality produces warnings, not rejection. Recovery points can
also come from adjacent accepted cells; visibility expansion must remove the
intervening cells before committing the patch. An optional `interiorAdmission`
predicate adds native material checks for introduced interior nodes and Tet4
centroids. It cannot bypass any exact PLC check. Unavailable recovery reports
the source-triangle indices still exposed to unfilled components, without
publishing a partial mesh.

```cpp
paacs::bnd::TetPlcCavityOptions limits;
limits.region = {3, solidTag};
const auto recovered = paacs::bnd::recoverTetPlcCavities(
    originalBoundary, partialTetMesh, limits, completedMesh);
```

`partialTetMesh` contains Tet4 cells only, with the original boundary's nodes
and classifications as its exact prefix. The completed mesh includes the
original Tri3 cells followed by the volume cells. An empty partial mesh can
handle suitable star-shaped PLCs. This API does not replace initial volume
construction, CAD projection, or CAD acceptance policy.

The main native mesher now appends `native_plc_cavity` after established routes
have refused, sharing the same cumulative work and output limits. PAACS's own
Delaunay builder proposes native-inside seeds; a bounded advancing front and
the public cavity API construct a complete PLC mesh. Tetrahedra whose centroids
are outside native material are removed and reconstructed with native admission
by default. The explicit approximation option below has a separately recorded
distance policy; strict native inclusion remains the default.

If remaining source triangles block recovery, the adapter can refine their
actual CAD wire chains, recheck every incident pcurve and the complete surface,
discard cells crossed by the refined boundary, and rebuild that local front.
Original boundary-node coordinates/classes and CAD face identities remain
fixed. The adapter explicitly refines the mesh triangulation; each subsequent
core call still treats its supplied PLC triangles as immutable. Trimmed arcs of
periodic curves with distinct endpoints are supported. Explicit closed-wire
recovery also supports a complete, uniquely ordered native parameter cycle,
including its seam. Ambiguous or incomplete cycles return diagnostics. All
incident triangles and actual shared-face coowners use the same new chain.

CAD27 now passes the main CLI. A retained run contains 5,678 Tet4 cells, 1,790
nodes and 3,290 boundary triangles. The original 1,631 boundary nodes remain
unchanged; eight new boundary nodes lie on evaluated native curves. Independent
archive validation, boundary auditing and complete cell-pair checks pass, with
no inverted cells, illegal contacts or remaining cavities. Native checks find
zero outside centroids, zero outside interior nodes and zero invalid boundary
node classifications. The maximum measured facet-sample error is `0.054583 mm`
against the unchanged `0.139937 mm` allowance. Selected volumetric samples are
not a complete geometric-error certificate. The physical gap remains empty.

The final retained end-to-end CLI run took about 232 seconds, including the
earlier routes and copied-geometry/source-immutability validation. A prior
integrated trial took about 195 seconds. These runs overlapped broader tests.
Construction uses bounded
time slices and logical work limits; exact element counts can depend on the
budget checkpoint. The original retained-front closure-only timing of 1.1
seconds and public API timing of 15.9 seconds excluded initial construction.
Very thin cells remain flagged (minimum mean ratio approximately `6.40e-12`);
the caller's strict quality requirement, when supplied, is still a final gate.

New regression checks cover external admission, transactional refusals,
trimmed periodic arcs, native source identity, cancellation and resource limits.
The complete original 51-file replay accepts 40 files, including CAD27. All 39
previously accepted archives are byte-identical, with no regressions. The other
outcomes are two meshing failures, five unsupported cases, three resource-limit
refusals and one invalid input. All inputs remain unchanged and no failed case
publishes an unexpected archive. The 80 registered tests and an installed
core-only API consumer pass.

Direct fallback trials on the other 11 files complete 11 of 34 individual
bodies (two in CAD14 and nine in CAD3636). These do not count as additional
whole-file successes; shared-interface assembly was not exercised. The linked
report separates work/storage refusals from source-surface, initial-front and
native-material gaps. In particular, CAD01 exposes introduced interior points
outside native material; production refuses that candidate rather than
publishing it. The one-page gallery uses only accepted whole-file archives from
the integrated replay.

See the main-mesher comparison (`out/visual-validation-2026-10-18/native-recovery/cad27-main-mesher-recovery.png`),
40-file gallery (`out/visual-validation-2026-10-18/native-recovery/coverage-40-cad-gallery.png`),
integration report (`evidence/local-2026-10-17/main-recovery/main-recovery-report.html`),
and [API example consumer](../standalone/plc-cavity-consumer/main.cpp).

## Preserved assemblies and diversified recovery

`NativeCadAssemblyPolicy::IndependentParts` is an explicit policy for retaining
the original parts of an overlapping assembly. The default remains
`ConformingDomain`. The CLI enables the part policy with `--independent-parts`
and always emits its JSON receipt in that mode:

```powershell
paacs-cad-mesh assembly.step 2 assembly.paacsmesh --native-cad --independent-parts --report-json
```

Each part is meshed independently on copied geometry. Mesh classifications are
mapped back to the original CAD identities, and interior nodes and tetrahedron
centroids are checked against the original body. Parts retain distinct node
identities; coordinate coincidence does not weld them. All parts must pass
before the destination or archive is published. A failure can report how many
parts were validated, but never publishes an incomplete assembly.
This explicit part policy currently supports Tet4; hex recovery remains separate
future work.

Native intersections are reported as `material_overlap`, `no_material_overlap`
or `unresolved`. An inconclusive Boolean is not evidence of disjointness.
Overlapping solids are not fused. `mesh_scope: independent_parts` and
`interbody_nonoverlap_certified: false` distinguish this result from a globally
non-overlapping domain. Its volume is the sum of independently audited parts,
including overlap. The archive stores geometry and CAD classifications, not a
validation certificate: retain the JSON receipt with it. A successful boundary
incidence audit alone does not certify global assembly non-overlap.

The reusable `TetPlcCavityOptions` API now supports fixed-seed interior proposals
through `maximumRandomCandidatePoints` and `randomSeed`. The default random
candidate count is zero. A seed changes proposed interior points without moving
the source boundary. `maximumMergedComponents` optionally lets one visibility
patch absorb adjacent unfilled components by removing intervening accepted
cells. Original PLC triangles remain hard barriers. Positive orientation,
complete boundary incidence, source identity and exact contact checks remain
mandatory. `patchesCommitted` counts patches; `cavitiesClosed` also counts the
additional components absorbed by those patches.

With `--native-kernel-admission` enabled, the native adapter uses at most two
diversified retries by default, configurable with `--recovery-retries 0..4` and
`--recovery-seed N`. Work, storage and output
allowances still apply. Source-barrier diagnostics guide native wire refinement
or bounded projection of selected facet centroids onto their actual trimmed CAD
faces. Every reconstructed surface and volume must pass the existing gates.
Introduced points outside native material are removed with all incident
tetrahedra before repair; deleting only their centroids is insufficient.

The iteration's retained evidence, commands, exact coverage comparison and
qualification results are in the
recovery/assembly report (`evidence/local-2026-10-18/recovery-next/recovery-assembly-report.html`).
The recovery cutaways (`out/visual-validation-2026-10-18/recovery-next/bounded-cavity-recovery.png`)
use a synthetic regression fixture and do not add a CAD-file success. Separate
periodic-wire experiments are retained in the evidence directory; they are not
part of the production candidate or a claim of completed CAD28 recovery.

The subsequent focused CAD28 iteration now completes all three original solids
through `--independent-parts --native-kernel-admission`. Its published archive
contains 11,576 Tet4 cells and 4,011 nodes. Independent per-part contact and
original-CAD audits pass; native interference checks report no material overlap.
This remains an independent-parts result: it does not establish a conforming
interface or a globally non-overlapping discretized material domain. The complete
51-file corpus has not been replayed for this iteration.

The new reusable surface repair detects opposite copies of the same triangle
on two CAD charts. It changes a bounded number of interior diagonals on original
exterior faces, preserving all native curve segments, node coordinates and shared
partition triangles. Every changed pair must meet the native face chord bound;
the repaired boundary must still pass closed-manifold and exact embedding checks.
Only the final original-CAD, positive-cell, volume and boundary audits permit
publication. CAD28's curved clamp needs one such diagonal change. The helper has
no CAD filename, body number or source-coordinate special cases.

Two bounded feature arrangements are tried before the older fine circular
fallbacks: four principal radial planes with native planar supports, then eight
concave-axis radial planes with up to sixteen failed-region cuts. They use the
existing native history and exact interface assembly. A private
`meshConstructionCandidate` admission allows at most `1e-4` relative native
split-volume drift, further capped at one percent of the requested mesh-volume
error allowance. This is a disposable construction aid, not a tolerance-certified
CAD subdivision. The default `1e-6` native partition admission is unchanged.
Original-CAD node/facet/material checks and the requested final mesh-volume
allowance remain mandatory. JSON reports expose `partition_volume_admission`,
`partition_volume_change` and `surface_diagonal_repairs` explicitly.

The focused report (`evidence/local-2026-10-18/cad28-focused/cad28-report.html`)
retains commands, source hashes, failure witnesses and independent audit receipts.
Its published-mesh cutaways (`out/visual-validation-2026-10-18/cad28-focused/cad28-complete.png`)
are exact sections of the actual exported tetrahedra. Element quality remains
poor in places; this establishes geometry coverage within the recorded validation
depth and tolerances, not readiness for simulation.

CAD3636's two tiny failing parts expose a different issue. Native nearest-point
queries find discrepancies of about `0.00351 mm` and `0.000308 mm` between
incident boundary representations. These exceed their respective volume/area
scales of `0.000426 mm` and `0.0000269 mm` (a scale estimate, not wall thickness).
Finer source tessellation plateaus without repairing the discrepancy. Tightening
tolerance metadata on a disposable copy also failed to produce a valid embedded
surface. Production does not silently change tolerances or omit these parts.
See the measured native-curve comparison (`out/visual-validation-2026-10-18/recovery-next/thin-parts-native-edge-disagreement.png`).

Guo and Zhang's
[Recovery techniques for finite element methods](https://arxiv.org/html/2412.03787v1)
concerns recovery of solution gradients and Hessians for error estimation and
adaptive computation. It is relevant to later simulation-driven adaptation,
rather than a geometric PLC boundary-recovery procedure.

## Retained project recovery reuse (2026-10-18)

The older project contains a bounded one-Steiner edge-cavity solver and a
staged-point transaction that were absent from this condensed package. These
are adapted into the current indexed recovery route. The existing-node solver
runs first; after a certified no-solution result the optional fallback tries
deterministic interior-point proposals. Both searches share the same work cap.
Every accepted replacement must retain the complete directed cavity boundary,
protected features, positive cells, cyclic target-edge star and legal contacts.

Core callers opt in with `TetEdgeGeneralCavityLimits::maximumSteinerCandidates`.
Native CAD callers use `--recovery-steiner-candidates 0..64` (default `0`),
alongside `--native-kernel-admission`. This extends the fixed-PLC recovery
route, rather than the separate advancing-front cavity filler. The native
adapter retains generated node identity during material extraction and checks
introduced nodes and cell centroids against both original and partition CAD.
No geometry healing, part fusion or quality-threshold relaxation is involved.

Retained local regression cavities establish the added recovery capability;
they do not count as complete CAD-file successes. The new core regression
also exercises protected-feature refusals, exhausted bounds, deterministic
replay and allocation failures. Existing routes remain selected when this
option is zero.

The local candidate replay checks do not establish whole-native-pipeline
bitwise reproducibility. CAD27 produced different interior triangulations on
separate runs of the previous frozen executable as well as the updated one;
a paired old/new replay produced the same archive on both. Its unchanged
source boundary and independently audited volume remain the acceptance gates.
This existing native reproducibility limitation is retained in the review
evidence rather than hidden by archive-hash comparisons.

## Connected region recovery

The optional core `TetPlcCavityOptions::patchSearch` refills a closed region
using its existing vertices and constrained outer triangles. It builds
tetrahedra from front faces, backtracks when a choice blocks completion, and
caches immutable candidate and pair-contact decisions. A pinched edge can
reopen its surrounding live tetrahedra as a group. Otherwise the region can
grow through bounded neighboring cell layers. Source triangles remain hard
barriers; other unfilled components require explicit merge allowance.

`maximumStates=0` disables this route in the core API. Candidate tetrahedra,
vertices, layers, replacement cells, work and working payload are bounded.
The implementation permits up to 256 local vertices and 512 replacement
tetrahedra; exhausting those bounds is a refusal, not proof that a CAD body
cannot be meshed. `useVisibilityExpansion=false` allows a caller to skip an
already exhausted single-apex trial. The search changes a local face ledger
at each step; it prepares the surrounding value object only for a completed
proposal. The outer caller publishes only a completely validated result.

The native final PLC fallback tries this route only after the earlier material
repairs refuse. `--recovery-patch-states 0..65536` controls the per-region,
per-layer state cap (default 2048). Its current native limits are 128 local
vertices, 16384 candidate tetrahedra and three neighboring layers. A zero
state cap disables the final route. Source geometry, surface coordinates,
native curve segments and shared interface triangles stay fixed. No invalid
part is skipped and overlapping bodies retain the selected assembly policy.

JSON `region_recovery` counters distinguish local patches from a completely
accepted body. Its timings include repeated PLC preparation and validation
calls, not just the final region search. The core partial-mesh audit uses an
x-axis broad phase for source-node containment and tetrahedron pairs while
retaining the same exact geometric predicates. Phase timings are elapsed wall
time; work counters and conservative payload admission are not CPU or RSS caps.

The core `centroidAdmission` callback can specify a separate cell policy;
leaving it empty retains the existing `interiorAdmission` rule for centroids.
Retained cells, visibility refills, kernels and region searches use the same
cell policy. Introduced interior nodes always retain `interiorAdmission`.
The native adapter forwards its explicit CAD chord policy through these repairs
and still checks the completed mesh's additional interior samples.

The local U-shaped and blocked-chamber tests establish concave refill,
connected-cell replacement, backtracking, preserved physical gaps and
transactional refusal. They do not establish a new whole-file CAD success.

## Longitudinal partition recovery

`--longitudinal-recovery` enables a bounded native partition/front fallback and
also enables native kernel admission. The equivalent API setting is
`useLongitudinalRecovery=true` together with `useNativeKernelAdmission=true`.
It runs after direct, coarse-partition and sweep attempts. The default remains
disabled, and `--no-native-partition` disables these additional cuts too.

The route combines native feature planes with 8, 16 or 32 lengthwise sections.
For a body without native cylindrical faces, it uses three Cartesian divisions
per axis plus the lengthwise sections, with the strict split-volume conservation
check and no sampled trim admission. Its report route is `native_cartesian_front`.
Each attempt starts from an owned copy of the original CAD and spends the same
job budget. Small pieces use a plain Tet4 advancing front and existing connected
region recovery; oversized or exhausted pieces cause another bounded attempt.
The Cartesian route allows up to fifteen seconds of local front growth before
region search, while the feature-axis route retains its three-second limit.
Both retain the same hard work, memory, output and cancellation checks.
Shared cut surfaces are tessellated once. Joining uses native entity identities
and exactly opposite triangles, never coordinate welding. Source curves and
piece boundary triangles remain fixed during volume construction.

The existing sampled native trim witness can admit a split-volume discrepancy
only after checking original supports, coverage, original tolerances and a
physical drift envelope. JSON records its sample count, physical drift limit
(`maximum_partition_trim_drift`) and measured original-tolerance ratio. This
witness is not a continuous CAD certificate: the complete
joined mesh must still pass original-CAD node, strict centroid, boundary,
embedding, ownership and volume gates. The new route does not enable the optional
centroid distance band. Quality warnings remain separate from mesh validity.

`--independent-parts` preserves each original assembly part and reports material
interference. The artificial subdivisions of one part join conformingly; an
interfering original assembly is not labelled a globally nonoverlapping domain.
The corpus runner forwards these controls and records working/region budgets.

## Actual geometry and mesh images

`paacs-native-cad-preview INPUT.step OUTPUT.json LINEAR_DEFLECTION` creates a
display tessellation from a copied BRep, applying face locations/orientations
and preserving face tags. It checks that the original model's triangulations
are unchanged. Point, triangle and JSON-output bounds apply to publication;
they do not bound the OCCT tessellator's internal allocation. Display geometry
is labelled separately from a supported PAACS volume mesh.

See [mesh_visual_report.md](../scripts/mesh_visual_report.md) for exact
cell-plane sections, cutaways, size/quality distributions, pinned plotting
dependencies and input/script/figure hashes. The scientific renderer reads
retained archives directly; it does not substitute illustrations for results.
