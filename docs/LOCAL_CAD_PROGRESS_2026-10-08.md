# Local CAD development progress — 9 October 2026

The following describes the frozen 39/51 checkpoint. Later
[wire recovery trials](CAD_RECOVERY_TRIALS_2026-10-10.md) and
[facet/kernel recovery trials](CAD_FACET_KERNEL_TRIALS_2026-10-10.md) add opt-in
construction methods and independently audited comparisons. They do not
increase whole-source coverage. The baseline CLI/core remain unchanged;
the current adapter has five selected fresh regression passes. The earlier
full-suite qualification below belongs to the frozen checkpoint.

The fixed original test set now has **39 whole-source mesh successes out of
51 files in one complete frozen-build, common-profile screen**. The preceding
native-wire-authority build reproduced 38/51. Its one regression was source
16: a failed component strategy suppressed an established whole-assembly
strategy. The current production dispatcher restores that strategy under the
same cumulative limits. The complete screen reproduced all preceding 38
accepted archives and the restored source-16 archive byte for byte.
No isolated split piece, surface-only diagnostic, repaired replacement or
private trial is counted as a whole-source success. The 40/51 target remains
outstanding. The historical [25/51 screen](LOCAL_CAD_COVERAGE_2026-10-03.md)
remains available for baseline comparison.

All 25 original baseline archives were reproduced byte for byte in the
current complete screen, preserving the earlier successful meshes.
The current eight-file correction passed 75/75 native tests (1,326.53 seconds
CTest elapsed time) and 15/15 freshly installed examples (2.83 seconds).
The unchanged core library retains its 4/4 selected core-only checks.
The native boundary test now contains 106 checks. Later volume
experiments are isolated in `build/` and local `evidence/`; they are not covered
by that production qualification.

The dispatcher tries the established whole-assembly construction once after
a typed component-construction refusal. It verifies unchanged destination and
source bytes before the retry, preserves all earlier work and phase timings,
and shares the remaining global admission. Cancellation, allocation failure,
overlap, invalid requests, unsupported controls, strict quality refusal and
internal/source failures remain terminal. The retry changes candidate kernel
selection; original-CAD ownership, material, sampled facets, interface,
embedding and volume checks remain mandatory. A portable thin annulus,
touching cap and isolated cube reproduce the issue: the component strategy
exceeds 64 pieces, while the whole strategy succeeds with 52 pieces. Tests
also cover one-less cumulative work, both strategies failing, strict quality,
cancellation and overlapping original bodies.

The retained additional successes are source folders 0, 2, 3, 15, 16, 21, 24,
909, 2727, 4545, 6363, 8181, 9090 and 9999. Complete original solids and body
identity are retained. Acceptance checks include finite positive cells,
connectivity and oriented boundary conformity, matched interfaces, exact
surface embedding, native material tests and a 1% whole-body CAD volume gate.
Native geometric checks use nodes, tetrahedron centroids and facet samples;
they are not a continuous CAD or exhaustive global tetrahedron-intersection
certificate. The preferred 0.25 quality value is a warning threshold.

## Current failures and development focus

The complete screen reports 39 successes, two failures, six resource-limit
refusals, three unsupported assemblies and one invalid source. The following
diagnoses combine its recorded failures with separately identified private
experiments; they do not imply every failure has a completed root-cause proof.

| Remaining source folders | Observed gap and next requirement |
| --- | --- |
| 1, 13 | Native edge/support or trimming discrepancies. A mesh must respect the existing owning-entity tolerance, or an explicitly reported geometry-repair variant is needed. |
| 10, 14, 3636 | Positive-volume original-body overlap. An explicit material/Boolean ownership model is needed before declaring the assembly a valid conforming mesh. |
| 25 | Native invalid BRep, including an unorientable shape. Preserve the source and distinguish a repaired variant from original-source coverage. |
| 11, 12 | General trimmed curved-surface topology, including a collapsed-pole/native-wire witness. Surface reconstruction and volume recovery remain unqualified. |
| 23, 26, 27 | Boundary recovery, native material conformity or partition preparation. A valid boundary alone does not establish a valid volume fill. |
| 28 | Shared curved-edge sampling, consistent face preparation and recovery of the complete contact assembly. Local repairs have passed; complete original-body publication is still required. |

The six top-level resource-limit reports must not be interpreted as exhaustion
of the common 300-million work cap. Their recorded work ranges from 123,719
to 7,685,096 units. The dispatcher retains a local partition/search admission
refusal while later bounded strategies continue; those later strategies also
report volume mismatch, invalid split shapes, lost history or native-kernel
failure. Several feature strategies refuse an excessive plane arrangement
before producing any pieces. Raising the global work cap would not resolve
these demonstrated geometry failures.

A separate bounded cut-selection experiment ranks authoritative native feature
and planar cuts by source face area, retaining Cartesian/principal cuts and
at most 24 or 40 total planes. It links the unchanged qualified libraries,
preserves original source bytes and keeps all existing partition and volume
gates. All four originals still refused: 11 reported native-kernel/invalid
split-shape failures, 23 invalid split shapes, 26 volume mismatch/invalid
split shapes, and 27 volume mismatch. The smaller arrangement avoided
plane-cap rejection but did not clear these geometry gates. No trial exhausted
its 40-million work cap or 600-second timeout. It published no mesh and was
not promoted into the dispatcher.

Case 27 has an independently replayed private surface with 1,962 prescribed
nodes and 3,936 directed triangles, Euler characteristic −6 and genus 4.
The earlier 1,809-node
surface had 208 explicitly typed NativeLine segments. Some additional straight
CAD wires are stored as two-pole splines. An exact collinear-pole check and
authoritative curve parameters extended the independently replayed diametral
test to 357 straight native segments. Existing coordinates and classifications
remained fixed. Subsequent bounded curved-face point insertion recovered more
internal diagonals, but the 32-attempt experiment still refused its volume fill.
Nineteen proposed interior points were requalified against the current closed
PLC, two exact parity rays and every authoritative straight-wire diametral
constraint. Earlier seed and planar conditioning trials remain recorded
refusals; none is accepted coverage.

A separate read-only probe found that a face-10 wire endpoint projects within
its raw owning-edge tolerance but outside the smaller raw face tolerance.
All ten existing face-10 surface points were reported OUT by the solid
classifier although whole-solid distance was zero on that exact owning face.
The private boundary-evidence experiment admits a new boundary point only through raw
owning-face projection, strict UV trim and whole-solid nearest support on the
actual sole owning face. Old nodes retain their own entity tolerances and
coordinates. Interior-point and tetrahedron-centroid material checks remain
strict. Its explicit cumulative ceiling is 1.5 billion work units, including
918,005,397 units retained before the new solve; earlier closed research scopes
are reported separately. It refused after 17 Core attempts and 16 further
CAD-evaluated boundary points; no archive was published. The final retained
surface has 1,978 nodes and 3,968 triangles. The remaining wire is a spline
whose evaluated midpoint differs from the physical chord by only
1.78e-15 mm, but its control poles are not exactly collinear. A new private
test permits authoritative midpoint evaluation without claiming exact
straightness, retaining the raw owning-edge chord bound, both actual pcurve
trims, protected interfaces and all volume checks. The cumulative work
ceiling remains unchanged. The evaluated-wire test recovered that wire and
advanced through 33 Core inputs, adding 30 further boundary points and two
validated quad flips. It retained a 2,008-node, 4,028-triangle surface but
still refused a complete volume; source bytes were unchanged and no archive
was published.

The next read-only probe identified two additional assumptions at wire 39:
faces 4 and 47 legitimately share three native edges, and a fixed source
vertex is offset from its straight spline wire within the vertex's own
tolerance. A nearest native curve midpoint differs from that mesh chord by
5.26424e-5 mm, has strict support on both pcurves and classifies ON. The private
four-attempt follow-up selects the actual edge through endpoint membership
on both uses and uses existing endpoint owning-entity tolerances for that
chord displacement. New points still require raw owning-edge support and
unchanged strict trims; all old vertices remain fixed. Its cumulative ceiling
remains 1.5 billion units, with 1,438,970,861 already charged before the solve.
These are private recovery experiments, not new public meshing routes.

The four-attempt endpoint follow-up recovered wire 39 using its existing
vertex tolerance (0.0001052849 mm); new curve points remained on their raw
1e-7 mm owning edge. It added three points and made one quad flip before its
attempt bound. Source bytes were unchanged and no archive was published.
The cumulative retained scope ended at 1,479,813,538 work units.

The normal CLI recheck also refused case 27. It stopped in partitioning, with
volume mismatch and lost-history diagnostics before the new recovery could
fill a volume. A private build tested the existing two-accuracy span-aware
native integration retry without enabling sampled trim recovery or changing
the original 1e-6 conservation limit. It also refused: the additional 276
recorded work units did not clear conservation or CAD-history checks. It
published no archive and preserved the source bytes. This numerical experiment
has no production or whole-corpus credit and was not promoted.

A separate unchanged-production probe extracted the default original-27
boundary (1,631 nodes, 3,274 triangles), bypassing Boolean partitioning without
changing the CAD. Its logical 128 MiB admission was too small for the existing
core estimate; an explicitly recorded 256 MiB probe entered recovery and
still refused. Diagnostic recovery identified edge [0,1], connecting native
vertices 1 and 2 across actual faces 1 and 13, with InvalidCandidate. That
original straight wire is about 19.011 mm long and has no intermediate default
surface nodes. A retained private prefix comparison also refused. The four
bounded probes consumed 18,016,704 of their separate 20,000,000 work units;
none published an archive. This separates boundary sampling and recovery from
the earlier partition failures; raising memory admission alone did not solve it.

A read-only inventory of eight unresolved originals found hidden planar
spline supports only at faces 74 and 75 in sources 11 and 12. Both also contain
207 spline surfaces and other curved categories, so two planar faces do not
prove a sweepable whole solid. No recognition change was promoted from this
inventory. General recovery must allow necessary additional points while
preserving the shared CAD boundary; this is distinct from optional element
quality improvement (see the primary [TetGen boundary-recovery explanation](https://www.wias-berlin.de/software/tetgen/1.5/doc/manual/manual002.html)).

The current private recovery test corrects a separate adapter issue: endpoint
membership in the same CAD curve does not make an internal face diagonal a CAD
wire. Actual cross-face incidence, or two explicit same-face seam uses, supplies
wire authority. Complete volume validation is still required after this change.

The main code now uses one private actual-use authority helper in both
artifact construction and recovery. Cross-face incidence establishes an
ordinary native wire; repeated nondegenerate uses preserve periodic seams.
An internal face diagonal whose endpoints belong to the same native curve
does not receive wire authority. Cylinder regressions exercise actual cap
diagonals, true rims and both seam branches. This correction does not add
points, alter tolerance policy or enable the later private recovery algorithms.

Case 28 uses a private whole-assembly experiment with consistent shared native
wire samples, one simultaneous incident-face preparation patch and a bounded
automatic repair of copied curved-face triangulation chords. A further local
cylinder construction is permitted only when OCCT's retained intersection
branch supplies no nearby stationary projection. Every proposed movement is
limited by the unchanged Original owning-edge and derived-entity tolerances,
all incident cut planes, face trims and native material checks. Original CAD
vertices are fixed. Complete body, interface, geometry, volume and source
guards remain mandatory; a local success is not coverage credit.

The latest complete case-28 trial refused at piece 139 after passing the earlier
stopped pieces. Investigation found that a nearest support projection can reject
a point whose exact owning-wire pcurve is still within the unchanged raw edge
tolerance. The private construction now checks the same authoritative Original
curve parameter, its strict owning-face UV trim and its physical image on both
faces; it does not widen face or UV tolerances. Sampling passed for the isolated
piece. Refining its shared cylindrical face reduced maximum sampled CAD error
from 0.00249991 mm to 0.00000896153 mm. The volume still failed: eleven
independently classified interior points did not resolve the shared curved
boundary failure. Recovery needs consistent boundary changes and volume
remeshing on both adjacent owners before another complete whole-assembly
result can receive coverage credit.

The production private session API now accepts an explicit bounded selection
of native edges to refine during preparation. It prepares every incident face
together on the same owned CAD copy, preserving matching shared keys and
opposite interface winding. Empty selection retains the previous output and
work exactly. Existing CLI calls keep that default behaviour; an automatic
edge-refinement retry is not yet connected. The two-body curved-interface
regression measured 9.76 times lower sampled surface error, with 64 to 200
matching interface triangles on each owner. Invalid selections, incident-face
and work limits, cancellation, source preservation and repeatability are
covered by the production tests.

The mesh figure renderer now computes the immutable model scale once per
drawing pass. Previously it scanned every node for every boundary face,
making dense cutaway rendering quadratic. Nine retained panels remained
byte-identical after the correction. The current common-profile gallery
independently measures 4,390,219 volume cells across all 39 meshes. It uses
verified byte-identical cached panels for the preceding 38 successes and a
fresh actual-cell section for source 16;
it does not count diagnostic surface meshes or add any new coverage success.

## Fixed set and wider corpus

The `abc_0000_step_v00` download contains **10,000
STEP files in 10,000 numbered folders**, one file per folder. Folder
`00000003` is already in the fixed 51-file set. The other 9,949 files form the
wider corpus; this inventory does not execute or qualify them.

Continue improving the fixed set first. Then use the same versioned mesher
and explicit resource profile on an unseen wider screen, preserving refusals,
timeouts, source hashes and any reported repair variants. A larger numerical
success count is useful evidence, not a claim that arbitrary CAD is supported.

The local toolchain uses MSVC's C++23 compatibility mode with `/std:c++latest`
and OCCT 8.0.1. It remains separate from the package's GCC/C++26 release
qualification. Logical working-memory admission does not cap OCCT's allocator
or process RSS.

## Retained artifacts

- `evidence/local-2026-10-03/coverage-40/contact-preservation-build-qualification.json`: current tested build, source and binary pins, and completed common-profile qualification.
- `evidence/local-2026-10-03/coverage-40/contact-preservation-complete-screen-qualification.json`: all 51 original-source receipts; 39 accepted meshes, preserved preceding 38 and restored source 16.
- `evidence/local-2026-10-03/coverage-40/contact-preservation-current-failures.json`: per-body and per-attempt reasons, work counters and source-preservation evidence for all 12 remaining originals.
- `evidence/local-2026-10-03/coverage-40/bounded-native-cut-originals/receipt.json`: four negative bounded cut-selection experiments against the unchanged qualified production libraries; no whole-source credit.
- `evidence/local-2026-10-03/coverage-40/contact-preservation-all-tests.xml`: current 75/75 native test results.
- `evidence/local-2026-10-03/coverage-40/installed-contact-preservation-validation.json`: current installed-consumer receipt.
- `evidence/local-2026-10-03/coverage-40/contact-preservation-baseline25-preservation.json`: byte-for-byte preservation of all 25 baseline archives under the current complete screen.
- `evidence/local-2026-10-03/coverage-40/contact-preservation-source-footprint.json`: comparable C/C++ source inventories; generated experiments are excluded.
- `out/visual-validation-2026-10-09/complete-common-profile39/index.html`: measured one-page gallery of the current complete 39/51 common-profile run.
- `evidence/local-2026-10-03/coverage-40/renderer-scale-cache-preservation.json`: nine scientific panels remained byte-identical after the drawing performance fix.
- `out/visual-validation-2026-10-04/coverage-40-retained-progress39/gallery.png`: actual retained whole-mesh gallery, with recorded build/profile provenance.
- `out/visual-validation-2026-10-08/source27-interior-point-diagnostic.png`: private surface and interior-point image; no volume-success claim.
- `out/visual-validation-2026-10-08/source28-shared-interface-diagnostic.png`: actual surface progress and rejected-cell diagnostic.
- `out/visual-validation-2026-10-08/production-shared-surface-preparation.png`: actual production two-body surface regression; no corpus-volume credit.
- `evidence/local-2026-10-03/coverage-40/corpus-inventory-2026-10-08.json`: full file inventory and unchanged fixed-source hashes.

No 40/51 or general-CAD completion is claimed in this progress report.

These local generated receipts, mesh archives and images remain in the
workspace. The development source ZIP excludes build products, OCCT/third-party
SDKs, raw CAD downloads and local experiment evidence.
