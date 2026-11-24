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

## Windows development qualification

The matched build passed 266 Python tests, with one expected POSIX-only skip,
all 104 STEP CLI cases, 94 native meshing assertions and the expanded native
boundary, curve-resource, interior-front, source-fallback and classification
tests. Python subprocess tests explicitly loaded the same source and native
extension. This is development-build evidence, not installed-wheel or
cross-platform qualification.

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
timing medians. The original 51-file score has not been rerun.

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
