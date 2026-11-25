# Source preparation and boundary-first construction

This development round follows the [initial construction qualification](GENERAL_CONSTRUCTION_QUALIFICATION.md).
It corrects curve sampling and parameter correspondence, and adds two independent,
experimental construction switches. Neither switch is enabled by default.

## Curve sampling

The adaptive curve integrator combines CAD continuity breaks with a small
uniform seed grid. On four edges of one pinned STEP input, an artificial seed
lay just one or two representable parameter values from a native continuity
break. The resulting interval could not supply the five distinct parameters
required by quadrature. The curve was smooth there; this was an artificial
partition problem, rather than a missing geometric feature or exhausted budget.

Native continuity breaks and endpoints now take priority. An artificial seed is
inserted only if both resulting intervals can supply the exact integration
stencil. Every native break remains present. Search and insertion work is charged,
and existing interval, memory, metric-error and geometric checks remain active.
A genuinely unresolved native interval still produces a numerical refusal.

Samples at a shared continuity break also belong to their native spans. The
left and right derivative limits can differ even though both positions agree.
The integrator retains both endpoint samples through the final arclength and
size-metric calculation; their same-parameter connector has zero length. Native
interior evaluation remains in use. Explicit local spline derivatives apply
only at internal span endpoints, including the second derivatives used for edge
curvature. Placements, numerical error thresholds and resource checks remain
active. Historical Bounded sampling retains its existing implementation.

Owned regressions use straight edges represented as cubic splines with Greville
poles. One has a continuity knot next to an artificial seed; another has two
mandatory continuity knots separated by three floating-point steps. This
separates a removable integration cut from a native interval the algorithm
cannot resolve. No imported filename participates in construction decisions.

## Pcurve correspondence

An imported edge's direct parameter map can miss its stored face pcurve within
the original source tolerance. Regular chart construction now makes the existing
bounded interior fitting method available independently of curvature sizing.
It searches only between the neighboring canonical sample brackets and must
meet the original source tolerance. Endpoints, degenerate uses and seam branches
retain their existing restrictions. Source geometry and canonical 3D nodes do
not move, and all subsequent chart, shell and CAD checks remain mandatory.

## Experimental routes

`native_source_fallback=True` permits the audited native tessellation provider
after a named geometric curve or chart refusal. It retains the live volume
constructor. Unknown failures, cancellation, invalid topology and resource stops
cannot select this fallback. Explicit face/box/body controls and surface point
caps are not silently bypassed. Source providers and chord retries share the
existing fixed attempt allowance.

A native source is not a global-size certificate. Its configured angle,
geometric checks and lack of verified source spacing are reported separately.
The source volume, oriented embedded shell and all final mesh/CAD checks must
still pass. A refused source does not authorize publication of a partial mesh.

`boundary_first=True` separately recovers the source boundary before inserting
interior points on the positive material complex. It does not force the live
constructor on or replace a verified structured route. Explicit interior site
and candidate ceilings remain upper limits; the existing post-recovery stage
also retains its 1,024-site algorithmic ceiling.

Strict serial boundary-first classification can use positive certificates from
the original trimmed CAD. Ambiguous queries retain the ordinary classifier.
Neither the geometry tolerance nor the final validity gates change. These
options are available in Python, saved recipes and the CLI; see
[Python usage](PYTHON.md) and [native options](NATIVE_CAD.md).

## Evidence and interpretation

Source attempts retain their provider and precise refusal witness. Material
classification reports cumulative requests and attempted point visits, plus
the last batch's completion and source-identity status. An interrupted native
query may be counted without yielding a decision. Optional serial timing divides
classifier setup, ray witnesses and ordinary queries; unprofiled stages are not
reported as zero-cost measurements.

Interior insertion retains private staged counts, work and elapsed time on
interruption. Those counters do not make a partial complex a published mesh.
The material and insertion phase times are contained in volume/total time and
must not be added to those parent timings again.

Local pilot requests use the same pinned source, global size, chord/volume
criteria, one meshing worker, one validation worker and resource limits. Every
attempt, including an unsuccessful retry, belongs in the cost. Short diagnostic
pilots do not update the original full-corpus coverage score. Downloaded CAD
sources and their figures remain in ignored local evidence; owned fixtures can
be redistributed with the tests.

## Paired original 51-file development qualification

The frozen Windows/Python 3.12 candidate `fb6997d0` was compared with
`be303775` on all 51 pinned STEP inputs, once per side in two separate profiles.
Both used size `0.15 × body-bounds diagonal`, growth 1.3, one meshing worker
and one validation worker, a 120-second cooperative allowance and a 150-second
process deadline. Native validity and CAD tolerances were unchanged. These are
development-build results, not installed-wheel or cross-platform qualification.

| Profile | Native and independently inspected, baseline → candidate | Native gains | Native regressions |
| --- | ---: | --- | --- |
| Default | 30/51 → 30/51 | None | None |
| Both experimental switches | 31/51 → 33/51 | CAD27, CAD7272 | None |

All 30 default passes were preserved. Fifteen Native-accepted attempts initially
hit the 60-second file-inspection allowance. Extended archive-only inspection
qualified all 15: two actual audits took 92.17 and 88.64 seconds, with 13 exact
proof reuses. Reuse requires identical archive bytes, numerical implementation,
options and limits; each attempt retains its own Native/CAD/source and tag
binding. No geometry gate or non-time resource limit changed, and no remeshing
was performed. CAD4149 already passed Native gates on both builds; the original
inspection timeout was not a meshing failure or new coverage gain.

The 90-second counts depend on the clock. Every counted mesh has completed
independent inspection. Native and pipeline clocks count their own work across
every generation attempt. The controller clock additionally charges the complete
referenced inspection proof, including earlier failed inspections:

| Profile | Native API call | Import + mesh + export | Full controller process |
| --- | ---: | ---: | ---: |
| Default, baseline → candidate | 28 → 28 | 28 → 28 | 23 → 24 |
| Experimental, baseline → candidate | 29 → 31 | 29 → 31 | 25 → 27 |

The target of 34/51 within 90 seconds has not been met. The two profiles have a
descriptive union of 34 accepted files; that is not an implemented dispatcher,
a single coverage score or a combined runtime result. The original 204 screen
attempts consumed 12,192.72 controller seconds. Including 12 held-out requests
and 40 fixed-slow repeats gives 256 generation attempts and 16,372.86 controller
seconds; the completed inspection follow-up adds 182.58 seconds of actual work,
counting shared audits once and retaining verification overhead. Per-case
qualification still charges the full referenced audit, so proof reuse does not
turn a slow qualification into a sub-90-second result. Native and pipeline clocks
remain unknown for 26 generation attempts; known subtotals are not full totals.

In the experimental CAD03 pair, Native time fell from 113.52 to 84.34 seconds
with the same recorded archive hash, quality and CAD errors. Import/mesh/export
took 85.58 seconds; controller-observed process time was 90.23 seconds. The
worker's shorter reported time is a different clock. These are single paired
observations, not timing medians. CAD27 passes the Native and independent checks
with 12,210 Tet4 cells at 93.71 Native seconds, but its minimum mean ratio is
about `4.32e-11`. Positive orientation and geometric validity do not make these
very poor cells suitable for simulation.

File 7272 also changed from refused to accepted, producing 7,176 Tet4 cells in
47.83 Native seconds with three source-facet additions. Its minimum mean ratio
is about 0.0226; it is another validity result, not a solver-quality claim.

Three repeats were completed for the five preselected slow files. Only three
files per profile qualified on both builds in every repeat, with identical
recorded archive hashes, achieved sizing, CAD errors and tags:

| Profile / input | Native median, baseline → candidate (s) | Controller median, baseline → candidate (s) |
| --- | ---: | ---: |
| Default / CAD20 | 79.90 → 81.55 | 87.97 → 89.56 |
| Default / CAD04 | 106.69 → 107.16 | 111.58 → 111.90 |
| Default / CAD2727 | 30.22 → 30.76 | 39.26 → 39.74 |
| Experimental / CAD24 | 44.88 → 41.15 | 50.60 → 46.87 |
| Experimental / CAD20 | 79.56 → 80.54 | 87.57 → 88.74 |
| Experimental / CAD2727 | 23.83 → 23.83 | 32.59 → 32.59 |

Default CAD24/CAD8181 and experimental CAD04/CAD8181 did not qualify across all
repeats. Their costs remain in the totals. The complete five-case median cannot
be reported, and halving it has not been demonstrated. The qualified pairs show
an 8.31% Native median reduction for experimental CAD24; the other matched
Native medians were approximately flat or 0.44–2.07% slower. These local medians
do not establish scaling or broad speed superiority.

A separate four-run CAD16 CLI trial compared stock queries with existing face-tool
reuse, leaving adaptive recovery and memoization disabled. Median CLI time fell
from 22.74 to 21.41 seconds (5.82%). All runs retained identical archive bytes
and CAD-error signatures; the archive also matches the Python matrix result.
All 54,088 projection queries remained. Facet-witness time fell from 6.73 to
6.05 seconds, while centroid checking showed no benefit. This is a modest
initialization saving under a separate CLI timing policy, not a new default or
an extrapolated corpus speedup.

Each build/profile accepted one of three preselected held-out files. That
success, CAD8645, had already qualified previously. CAD3393 and CAD5235 were
first executed this iteration and failed. This is useful failure evidence,
not general unseen coverage.

The implementation changes address named stages. Source facet projection and
swept-patch remapping repair original-CAD material disagreement while preserving
unaffected cells and all subsequent audits. Canonical chart repair retains
shared curve nodes; inconsistent original trims still refuse, including CAD13.
Bounded source-coverage diagnostics and exact-bit trim memoization reduce
repeated core work. The earlier source-first probe remains Advanced-only
because its pilot did not improve total time; these changes are not a quality
optimization or a default-mode promotion.

The [numeric record](data/source-preparation.json) retains historical pilots,
all three clocks, unknown measurements, exact proof identities and the completed
repeat results. The later diagnostics-only commit `c829ba34` is separate from
the frozen matrix binary. Its Python STEP wheel and mesh qualification workflows
both completed successfully. No merge, tag or package release follows from this
development qualification.

## Historical development qualification

The base build at `be303775` passed 266 Python tests, with one expected POSIX-only skip,
all 104 STEP CLI cases, 94 native meshing assertions and the expanded native
boundary, curve-resource, interior-front, source-fallback and classification
tests. Python subprocess tests explicitly loaded the same source and native
extension. This is development-build evidence, not installed-wheel or
cross-platform qualification.

The subsequent cavity/chart candidate passed all 32 selected C++ regression
suites, including 2,156 boundary checks, 449 facet-refinement checks and 32
source-coverage checks. Its matching Python run passed 266 tests with one
expected platform skip and 82 redirected subprocesses; all 104 STEP CLI checks
passed. These regression results do not replace the paired 51-file qualification.

The production sampler completed all 500 curves in three pinned inputs. Staged
diagnostics recovered four CAD03 curves at the artificial-cut fix and nineteen
CAD13/27 curves at the derivative-side fix. These are source-stage results,
not a paired full-corpus baseline. Full volume results remain separate:

| Input | Result with both experimental switches | Tet4 | Native seconds | Remaining gap |
| --- | --- | ---: | ---: | --- |
| CAD03 | Accepted; archive orientation, ownership and tags also checked | 8,931 | 117.44 | Exceeds the historical 90-second target; poor cells remain |
| CAD13 | Refused without publishing | — | 24.50 | Crossing chart constraints; alternate source is nonmanifold |
| CAD27 | Refused without publishing | — | 83.96 | Remaining volume cavity after material classification |

Each request used size `0.15 × body-bounds diagonal`, a 120-second cooperative
allowance, a 150-second process deadline, one meshing worker, one validation
worker and the same geometry/resource settings. There were no external retries
within these requests. Earlier failed development attempts are also retained in
the [numeric record](data/source-preparation.json); single samples are not
timing medians. These historical pilot results do not replace the later paired qualification above.

CAD03's regular chart source passed on its first attempt; it did not need native
source substitution. Boundary recovery resolved the remaining five edges and
ten faces in 0.45 seconds. Surface construction took 46.03 seconds and material
classification took 59.35 seconds, including 59.06 seconds of ray witnesses.
Those are measured remaining costs, rather than evidence that cavity recovery
dominates this successful run. The final boundary sample error was 0.08753 mm
against a 0.10262 mm limit; relative volume error was 0.3516% against 1%.
Minimum and median mean-ratio quality were approximately 0.00166 and 0.44059.
Acceptance establishes the recorded validity/geometry checks, not solver adequacy.

An owned in-memory prism with a genuinely unresolved three-step native knot
interval separately exercises the source-provider fallback. It produces 495
Tet4 cells, preserves tags/source bytes, and survives archive reopening and a
boundary audit. Minimum mean ratio is 0.0227; source spacing remains unverified.
It is a recovery regression, not a typical-quality benchmark. Its STEP
roundtrip is unqualified and is not counted as STEP coverage.

## Diagnosed source recovery

A closed polygonal boundary can pass its chord and volume criteria while a
small volume next to a concave curved face cannot satisfy original-CAD material
admission. Increasing the cavity search budget cannot repair that discrepancy.
The experimental boundary-first route can refine diagnosed source facets by
projecting their centres onto their original trimmed CAD faces. This is geometry
recovery; it does not introduce a user sizing region or relax CAD tolerances.

Each proposed point retains the original trim, evaluation, displacement and
triangle-foot requirements. Unavailable proposals and planar no-ops can be
omitted from this private candidate set. The accepted subset is committed only
after complete surface closure, embedding and shared-owner checks. Cancellation,
allocation failure and work, node or memory exhaustion still refuse the batch.
The route has at most four rounds and 32 added source nodes per round. Explicit
per-face sampling caps remain hard limits.

After a committed `ABC → ABP / BCP / CAP` edit, the changed domain is contained
in `conv(A, B, C, P)`. The local remapper tests retained tetrahedra against these
sweeps with the existing identity-aware geometry predicates. It preserves legal
contacts, drops affected cells and keeps other coordinates and body ownership.
This produces a partial candidate, not a completed mesh. The following recovery
still performs its full source, contact, centroid and cell-pair checks, followed
by the original final CAD, volume and mesh audits.

Standard recovery diagnoses these source edits after its bounded connected-region
attempt. Advanced recovery can probe them earlier, before the larger region
search. The early probe is experimental: one local pilot completed correctly but
used more source points without reducing total time. It is not the Standard
default. Both routes retain the same publication checks and explicit limits.

The source-shell receipt describes the initially admitted provider. Subsequent
committed source changes are reported separately: `surface_triangles` describes
the latest retained source, and `recovery_source_facet_nodes` counts cumulative
committed private additions across attempts. These additions are not necessarily
all present in the final mesh. Source geometry refinement, partial remapping and
source validation times are nested in native volume time. Boundary-error maxima
can include earlier admitted source configurations and are conservative rather
than final-only measurements.

Shared-edge chart repair separately bisects only canonical curve intervals whose
polygon chords cross. All incident faces reuse the inserted curve nodes. The
existing native coordinates, face/body identities and construction limits stay
in force. Inconsistent stored 2D trim curves are a distinct issue: a successful
projection onto an underlying surface does not by itself establish membership
in that original trimmed face. The original trim gates remain enabled.

Original-CAD ray admission also retains its existing two direction-pair
certificates. Analytic faces can share one line solve between eligible opposite
halves. Nonanalytic faces retain their original directed solves, and a refused
positive half avoids unnecessary deferred negative queries. Margins, trim and
transition tests, ambiguity handling and ordinary-classifier fallback remain
unchanged. This reduces repeated work; it does not certify element quality or
solver accuracy.

Surface construction can also reuse identical original trim queries within one
face. Keys retain the exact floating-point parameter bits, and misses call the
original classifier. The optional cache holds at most 4,096 entries per active
face. Its storage is included in the working-memory receipt and released before
otherwise admissible geometric growth needs the space. It does not reduce the
existing geometric planning allowance or change trim tolerances.
