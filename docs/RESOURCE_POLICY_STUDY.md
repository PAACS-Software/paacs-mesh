# Resource-admission study

Progressive admission removes a reproduced obstacle to initial interior-point
construction. Focused cylinder runs show better shape distributions under the
same explicit caps, and finer achieved wall spacing when those caps are unset.
They also expose remaining memory stops, poor individual elements and long
construction times. The default remains `resource_policy="bounded"` pending
qualification. Experimental `auto` and `thorough` scheduling remain opt-in.

These are unreleased development results. A package version of `1.2.0rc3` in a
local report does not identify the tested native binary. This study distinguishes
the frozen parent build from the later P6 native snapshot; source, binary,
archive and saved-report fingerprints are retained with the local evidence.
It does not describe a new published wheel.

The [implementation and qualification plan](RESOURCE_POLICY_PLAN.md) defines
the remaining checks. [Construction controls](CONSTRUCTION_BUDGETS.md) and
[diagnostic contracts](DIAGNOSTIC_SCHEMAS.md) explain the request and receipt
fields. The [earlier construction study](CONSTRUCTION_BUDGET_STUDY.md) concerns
different snapshots and must not be combined with these rows as one timing
experiment.

Local measurements use Windows 11 build 26300, an AMD Ryzen 9 8945HS with
16 logical processors and 31.3 GiB of visible physical memory, CPython 3.12.6,
OCCT 8.0.1 and MSVC 19.44 with strict floating-point compilation. Each attempt
uses one meshing worker and one validation worker; timed mesh cases run
serially. These are observations on one host. A timing median requires at
least three complete repetitions with matching source, binary and mesh hashes.

## Defect and corrections

The old initial-point auxiliary admission coupled its permitted future cell
count to bootstrap storage and import-work reservations. A larger ceiling could
therefore prevent the same small initial seed from being admitted. Raising
global work alone did not address that dependency. The parent fixed-cap cylinder
runs below admitted no initial-front points and reported `auxiliary_bootstrap`.

Normal admission now constructs the initial numerical seed once within the
available work and temporary-memory allowance, then admits the returned cells
and their actual import work. It separates the future ceiling from currently
reserved storage. Capacity grows between candidate batches; compaction and
growth retain allowances for old and replacement arrays while both coexist.
Capacity-only growth preserves IDs, while compaction invalidates insertion
plans and rebuilds affected queue state. Resource retries revisit the proposal
instead of interpreting a storage refusal as a failed geometric candidate.

The bulk queue also grows by demand. The enclosing front, queue and auxiliary
complex share a modeled allowance. Initial placement no longer reserves a
future volume front or apex cache that is created after its temporary data is
released. Later recovery retains a separate input/workspace envelope. Normal
curve sampling uses adaptive numerical integration and model/resource-derived
segment and face limits instead of the historical fixed curve and face-count
ceilings. Shared-edge identity, trims, orientation and subsequent CAD checks
remain required.

Completed work receipts are collected before propagating an observed
cancellation. Ordinary incremental callbacks still check cancellation before
another operation. Explicit work, memory, node, cell and output limits remain
active; Normal with omitted work keeps accounting and an overflow guard without
a practical logical-work ceiling. Its omitted modeled memory remains 512 MiB.
Logical work is not CPU instructions, and a stage-memory estimate is not process
resident memory or a hard operating-system limit.

Focused native tests cover bootstrap independence from an unused larger
ceiling, transactional growth/refusal, cancellation, ID/plan handling, adaptive
curves exceeding 256 segments and a valid 112-face boundary. These checks do not
establish general CAD coverage or simulation adequacy.

## Curve metric correction after coarse-case diagnostics

The new adaptive curve route also introduced a numerical compatibility
regression on some coarse inputs. In the P6 capped comparison, CAD06 and CAD08
abandoned their regular surface trials during shared-edge sampling with
`nonconforming_edge`. The parent accepted those trials. Normal instead selected
partition-star or general cavity fallbacks, changing both surface and volume
meshes. Initial interior-point placement had not been attempted, so these
particular first-stage failures were not evidence against its live numerical
complex or volume constructor.

An isolated, source-pinned curve probe identified `curve_parameter_resolution`
on 24 of CAD06's 33 shared edges and six of CAD08's 41. Every refused curve had
a valid historical bounded schedule under generous unused resource allowances.
The failure was not a work/memory stop or a nonfinite-data refusal.

Two sizing effects caused relative quadrature error to persist through
subdivision. Vertex feature sizes applied only at exact endpoints created a
pointwise jump to the immediately interior size. On remaining CAD06 curves,
the final rejected samples exposed curvature-derived size jumps a few ULPs
before nominal C2 breaks. Derivative-speed, length and geometric-mapping checks
passed there; the sizing-metric error did not. Sampling only one ULP to either
side of the rounded C2 break had missed those transitions.

The correction prepares a continuous endpoint feature envelope and lowers each
five-sample stencil under the existing growth slope before estimating metric
error. Raw cached samples remain unchanged, and the final global arclength
growth pass still runs. Physical chord distance is a conservative bound on
arclength for these local envelopes. This deliberately estimates a prepared
metric instead of trying to certify integration of a discontinuous raw sizing
field. It adds charged bounded traversal; it does not change source geometry,
weaken geometric/validity checks or introduce a permissive fallback.

After this correction the isolated probe completed all 33 CAD06 and all 41
CAD08 shared curves. An owned wide-parameter line joining unequal-radius
features passed an independent analytic metric/spacing reference. An owned C1
spline with a curvature jump passed a dense analytic-arclength reference, with
metric agreement within 0.1% and per-segment metric-spacing agreement within
1%. Fine circle and Bezier cases retained their 629 and 664 segment counts;
the historical Bounded schedule checks also passed.

These are curve-level and focused regression results. P7 full-mesh results,
quality/timing comparisons and the subsequent corpus rerun remain pending.
Successful curve sampling alone does not establish that either changed volume
mesh has recovered the parent's quality or runtime.

## Completed parent corpus baseline

The frozen parent completed all 51 selected input attempts with one meshing
worker, one validation worker and one case at a time. The profile uses Tet4,
advanced recovery, size factor 0.15, an explicit 300 million logical-work cap,
200,000 nodes, one million cells, 128 MiB output/audit allowances and a
90-second process deadline. Inputs are pinned by path, size and SHA-256.

| Parent outcome | Files |
| --- | ---: |
| Reported successful publication | 32 |
| Process deadline | 15 |
| Unsupported | 2 |
| Invalid geometry | 1 |
| Other failure | 1 |
| Total | 51 |

The separate archive screen passed all 32 published meshes, examining
**3,919,286 positive Tet4 elements**. Its scope is archive structure, exact
orientation/shape measurements and recorded body-identity checks against the
saved evidence. It did not re-import CAD, repeat native geometric validation,
perform a global intersection audit or demonstrate solver convergence. The
remaining 19 attempts published no mesh to screen.

The P7 paired candidate corpus and its independent archive screen are pending.
The parent count alone does not establish preserved coverage, gains or
regressions for Normal. Nested stage timers also cannot be added together as
independent costs; exclusive timings and every retry must be retained in the
completed comparison.

## Fixed-cap cylinder comparison

The owned fixture is a radius-1 mm, height-2 mm cylinder. Global size is 0.5 mm,
the cylindrical wall request is 0.1 mm and growth is 1.3. No volume sizing region
is added. Both runs request `construction_effort="legacy"`, surface/chord caps
of 1,024/384 per face, 65,536 interior candidates, 2,048 initial sites,
100 million work units and 512 MiB modeled memory. Both use one meshing worker
and one validation worker.

| Build/policy | Tet4 | Mean ratio below 0.1 | Minimum / median mean ratio | Wall median / p95 / maximum edge (mm) | Native seconds |
| --- | ---: | ---: | --- | --- | ---: |
| Parent, bounded | 4,581 | 851 | 0.0000762 / 0.214 | 0.102748 / 0.172800 / 0.294896 | 24.950 |
| P6, normal | 8,326 | 10 | 0.018282 / 0.702 | 0.103403 / 0.172576 / 0.318633 | 20.704 |

Native seconds include meshing and its validation, excluding import, export and
inspection. These are single attempts, not repeated medians. Requested caps
are comparable, but surface triangulations and volume meshes differ. The time
change is not a comparison of equal element counts or identical achieved
volume meshes.

The P6 initial front admitted 630 points, then stopped at its allocated work
limit after 9.384 seconds. Its source sweep and bulk queue did not complete.
The cylindrical face still reached its 1,024-point cap with 765 pending entries.
The wall median changed little, and its maximum edge became longer. This result
supports the admission correction and an improved measured shape distribution;
it does not establish complete requested sizing or improvement in every metric.

## Inherited Normal cylinder probe

A separate P6 attempt requests a 0.08 mm wall, the same 0.5 mm global size and
growth 1.3, with stage caps unset and Normal's inherited work/memory policy.
The effective modeled memory is 512 MiB. This probe explicitly allows
180 native seconds; its **112.093-second** result is not a 90-second pass.

| Measurement | Observed result |
| --- | ---: |
| Tet4 / nodes | 29,895 / 6,279 |
| Cells with mean ratio below 0.1 | 1 |
| Mean-ratio minimum / p01 / median | 0.010309 / 0.403815 / 0.850844 |
| Minimum dihedral angle | 0.057035 degrees |
| Minimum scaled Jacobian | 0.00100357 |
| Unique wall-edge median / p95 / maximum | 0.079988 / 0.083605 / 0.100709 mm |
| Initial-front time / logical work | 81.099 s / 371,955,983 |
| Largest recorded combined placement estimate | 536,869,376 bytes |

All three surface construction queues exhausted with no pending entries. The
initial front completed its source sweep, but its bulk queue remained
unfinished and its stop reason was `combined_working_memory`. Queue state and
the historical unresolved-void ratio are not final sizing certificates. The
wall maximum still exceeds the request, and no complete growth or solver-error
certificate is inferred from the edge population.

The published result passed its recorded native validity checks. Its maximum
sampled boundary error was 0.001189 mm and relative volume error was about
0.0794%. These geometric measurements are separate from the very small minimum
dihedral angle and scaled Jacobian. One element below a threshold does not
describe the severity of that element or establish a simulation-ready mesh.

## Quality warnings and validity

The cylinder benchmark omits both `warning_quality` and `minimum_quality`.
Its saved native reports explicitly record a zero warning threshold,
`warning_enabled=false`, zero required minimum and null requested thresholds.
Consequently, `quality_warning=false` does not mean that poor-shaped elements
were absent. Positive orientation and a completed native audit are distinct
from a requested shape-quality target.

An earlier exploratory snapshot produced an exact-positive, almost-coplanar
cell with independently measured mean ratio around 1.4 × 10⁻¹⁰. Its ordinary
edge lengths concealed its almost-zero volume; validity alone did not establish
solver adequacy.
That observation remains a quality defect to investigate, not a reason to
discard orientation or CAD-validity checks. Native aggregate and independent
per-element metrics are retained separately; near cancellation can affect
their floating-point magnitudes without changing the intended mean-ratio
definition.

## Pending qualification

| Evidence | Status |
| --- | --- |
| Frozen parent 51-case run and independent archive screen | Complete: 32 publications and 32 archive-screen passes. |
| Candidate 51-case run, gains/regressions and independent screen | Pending. |
| Repeated work/memory grid with every refusal and retry retained | Pending. |
| Box, sphere, through-hole and oblique scaling cases | Pending broader qualification. |
| Fixed-camera cutaways with verified archive cell IDs | Pending for this study. |
| Clean installed-wheel and platform qualification | Pending. |

The [plan](RESOURCE_POLICY_PLAN.md#reproducible-experiment-plan) supplies the
repeatable grid and inherited-policy commands. The
[trial runner](../benchmarks/construction_budget_trials.py),
[corpus runner](../scripts/run_cad_corpus.py) and
[archive screen](../benchmarks/audit_resource_corpus.py) retain different scopes
of evidence. Sanitized results and provenance for the completed paired study
remain to be added. No default promotion, release-tag change or package
publication follows from these focused observations.
