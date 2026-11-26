# Assembly diagnostics

Body geometry, separate part meshes and shared mesh interfaces need distinct
results. `assembly="preserve_parts"` retains original body identities and reports
material interference. Separately valid parts can overlap; their acceptance does
not establish a globally non-overlapping volume domain.

`CadModel.body_relations(pairs, ...)` inspects explicitly selected body pairs.
It distinguishes separation, face/edge/point contact, material overlap and
unresolved queries. The API accepts at most 64 pairs per request and retains
original model identities and placements. A geometric contact classification
does not certify matching surface triangles, a bonded interface or contact
physics. Incomplete query reports retain returned evidence and leave unqueried
pairs unknown. The time allowance is cooperative, and tracked memory excludes
some native-kernel allocations. See [Python usage](PYTHON.md) for the public API.

## Completed assembly preservation cohort

The final local build at `bea63fa` completed all 21 pinned profile/input
requests. The selection uses three pilot receipts, six qualified prefix receipts
and twelve tail receipts. The Native binary, runtime and original requests
match across these runs; reviewed controller versions differ in bounded I/O
handling. The failed Default CAD29 export is excluded from qualification and
its full process cost remains charged alongside the successful retry.

| Profile | Previous acceptances preserved | New Native acceptances | Native regressions |
| --- | ---: | ---: | ---: |
| Default | 8/8 | 0 | 0 |
| Experimental | 8/8 | 0 | 0 |

The experimental profile also retried five previously failing assemblies; all
five still cancelled. These results preserve 16 existing Native acceptances
across the two profiles. They do not update the original 51-file coverage score
or establish independent archive qualification for every accepted attempt.
Selected process time totals 1,013.052 seconds; including the failed export
gives 1,033.196 seconds. Internal method retries are included in each call.

| Input | Remaining observed gap |
| --- | --- |
| CAD00 | Two of six parts validated. First-child AUTO refusal took 50.844 s before a 10.072 s Legacy success. A later child returned successfully, but its outer original-body validation did not complete. |
| CAD01 | Source generation took 70.611 s and volume work 43.956 s; one boundary face remained unrecovered at cancellation. |
| CAD10 | The first child succeeded in 55.426 s with 20,296 child centroid admissions. Its outer original-body checks remained unfinished. Reported part interference is descriptive under `preserve_parts`, rather than the refusal reason. |
| CAD14 | Preparation/interference queries consumed the 144.530 s process before child meshing. Seventeen of 21 pair rows returned; four unqueried pair outcomes remain unknown. The unfinished preparation timer is unavailable, rather than measured zero. |
| CAD28 | Source generation took 9.730 s and volume work 109.852 s, including 94.837 s in regional search. Twenty-four edges and thirty faces remained missing at cancellation. Material classification had completed. |

These are deadline-limited observations. Child success does not finish the
wrapper's original-body checks, and nested phase times must not be summed with
their parent call. The remaining inputs are not classified as permanently
unsupported.

## Local validation and mesh quality

The coherent Windows build passed five focused Native suites in 22.38 seconds
and 88 Python tests in 45.038 seconds. This is local development validation;
it does not establish installed-wheel, cross-platform or completed CI results.

CAD05's final recorded archive is byte-identical to its independently inspected
historical and v5 archives: 3,200 Tet4 cells and 1,190 nodes. That exact-payload
orientation, incidence, tag, quality and achieved-size evidence can be reused;
each build retains its own Native original-CAD acceptance.

CAD16's experimental mesh passed independent exact orientation, internal-face
winding, exterior and body/face tag checks, alongside its historical comparison:

| Measurement | Historical | Final |
| --- | ---: | ---: |
| Tet4 cells | 9,016 | 9,208 |
| Median mean ratio | 0.1835 | 0.2728 |
| Cells with mean ratio below 0.1 | 2,628 | 1,639 |
| Minimum dihedral angle (degrees) | 0.02813 | 0.001392 |
| Public mesh call (seconds) | 21.177 | 29.656 |

The median and low-quality population improved while the worst angle worsened.
These single attempts demonstrate mixed quality and longer runtime, with changed
achieved geometry. The file checks inherit original-CAD validity from the Native
receipts; they do not prove global/interbody nonintersection or solver adequacy.

CAD20's shorter calls also changed the approximation: Default Tet4 count fell
from 28,930 to 13,860, and the experimental count to 13,897. Measured CAD errors
increased within the existing limits. Equal requested size therefore supports
no equal-achieved-size speedup or fair-halving claim. The numeric record retains
the counts, body errors and clocks separately.

## Isolated-body evidence

A private diagnostic harness selected body 2 from the six-body CAD00 input. It
copied the body with complete topology-history mappings and checked bitwise
world placement, source/child BREP guards, original-body boundary geometry and
original-body centroid admission. Export remapping changed classifications
while preserving coordinates and connectivity. The whole-file global size
remained 91.5599 mm.

The original automatic chord policy produced a Native-accepted body mesh in
18.632 process seconds, including 7.157 seconds in the child mesher. Its 4,346
Tet4 cells have Native-reported minimum mean ratio `4.21e-7` and median 0.126.
This is severe poor quality despite passing the geometric and validity gates;
it does not establish simulation adequacy or acceptance of the six-body input.
A separate explicit global chord request failed after 5.346 process seconds.
That request changed the original automatic chord semantics and remains a
separate failed variant. The initial loader failure also remains in the costs;
it occurred before meshing.

A controlled pair used the same executable and produced identical archive bytes:

| Original-body validation mode | Centroid admission (s) | Complete process (s) |
| --- | ---: | ---: |
| Ordinary classifier | 10.975 | 19.585 |
| Opt-in trimmed-ray witnesses | 0.488 | 8.967 |

The ray route retains original-CAD evidence and the ordinary fallback for
ambiguous queries. This single pair isolates validation cost; it does not
qualify a default change, timing median or assembly-wide speedup. Replaying the
old whole-assembly dispatch still cancelled at 120.050 Native seconds and
120.871 process seconds.

## Dispatch preservation

The whole-assembly AUTO setting had been resolved before individual-body
dispatch, preventing a child from using the same route as its standalone
request. Simply deferring AUTO to every child also changed established source
selection. In the first preservation pilot, CAD05 changed from a Native
acceptance in 3.238 seconds to cancellation after 120.021 seconds. Its source
provider changed from native tessellation to regular charts; the requested
automatic chord had not tightened.

The cohort stopped after that regression. An already-started CAD15 request was
aborted with unknown Native outcome. Its corrected elapsed wall time is 37.293
seconds; the original 38.100-second receipt remains preserved with an erratum.
A monotonic duration for that aborted request is unavailable. The completed
CAD05 process and aborted CAD15 request cost 158.180 seconds together.
A subsequent v5 pilot completed two attempts, followed by ten of nineteen
planned remainder attempts. All eight Default preservation cases and Live
CAD05/CAD15 retained Native acceptance. Live CAD16 then refused with a geometry
error after 22.839 Native seconds and 23.763 process seconds, so the remainder
paused. The pilot's Live CAD00 request still cancelled near 120 seconds. These
12 completed development trials cost 468.442 process seconds; the nine remaining
requests are unexecuted and unqualified. Their receipts remain separate from a
complete selection. These historical trials establish no additional whole-file
coverage by themselves.

A v6 pilot subsequently retained Live CAD16 acceptance in 29.656 Native
seconds. The three-case pilot completed, and the remainder reached seven
requests before an evidence I/O failure on Default CAD29. Its saved Native
report records successful geometry, audit and quality gates for both bodies,
but an atomic worker-receipt replacement raised `PermissionError [WinError 5]`.
The official worker receipt remained `running`, no archive was exported and
final worker guards were incomplete. The derived erratum retains every raw
receipt and temporary-file hash. This is an unqualified infrastructure attempt;
a Native algorithm regression is not established. Its 20.144 process seconds
remain charged. The ten v6 development attempts cost 357.272 process seconds.
That failed export remains separate from the completed selection above.

## Contact-map progress

Initial selected-pair queries completed for CAD00, CAD01, CAD10 and CAD28.
CAD10 reported one face contact, one material overlap and one unresolved pair.
The initial CAD14 query hit a 90.016-second process deadline, leaving all 21
pair outcomes unknown.

A later CAD14 progress-preserving query returned 15 rows: 11 separated pairs,
three face contacts and one unresolved pair `(3, 7)` with
`material_query_failed`. One material Boolean had started for that pair. Six
unqueried pairs remain unknown. The final source guard passed, and cancellation
returned without a hard kill. However, the requested 30-second cooperative
allowance took 81.753 Native query seconds and 84.256 process seconds. Prompt
cancellation remains unresolved. This is partial diagnostic progress, without a
complete map, conforming interface or mesh-coverage gain.

## Bounded intersection probe

A private query isolated CAD14 pair `(3, 7)` and compared oriented bounding-box
broad-phase pruning off and on. Both cancelled without establishing a material
relation. Intersection time was 71.733 versus 71.508 seconds, with roughly
71-second progress gaps inside edge-face intersection. This corrects the earlier
unconfirmed face-face hypothesis for this pair. The probe did not demonstrate a
useful improvement or identify the specific inner geometric operation.

Both source and copied-body guards passed, and both processes shut down without
a hard kill. All three process costs, including a 1.071-second loader failure
before the two queries, total 149.400 seconds. The
[separate probe record](data/cad-intersection-probe.json) retains the source,
runtime and receipt hashes, cancellation observations and unknown outcomes.
No production Boolean option, default, geometry tolerance or validity gate
changed because of this experiment.

The [numeric record](data/body-interface-diagnostics.json) retains each attempt,
exact artifact hashes, verified worker-memory observations, unknown measurements
and separate nested clocks. The completed diagnostic and meshing attempts,
including development cohorts, final selected attempts and the intersection probe,
consumed 2,146.099 seconds; adding the corrected aborted wall time gives
2,183.392 seconds. Receipt identities prevent repeated selected attempts and the
excluded I/O failure from being counted twice. Historical references and unknown
preparation/compile costs remain separate. Raw CAD, BREP snapshots, meshes and
figures remain in ignored evidence.
