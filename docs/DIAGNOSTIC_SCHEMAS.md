# Diagnostic contracts and issue codes

This reference includes unreleased development interfaces. Discover an installed
build with `paacs-mesh capabilities --json`; package version alone does not
identify a locally built development artifact. Preserve its source commit too.
Readers should accept additive keys, check `schema`, and refuse unknown schema
versions when they cannot interpret the evidence safely.

## Capability discovery

`pm.capabilities()` and the CLI return `paacs-mesh-capabilities-v1`. The object
contains `package_version`, `engine_version`, `units`, `operations`,
`cli_commands`, formats, selection predicates, controls, recovery and assembly
restrictions, diagnostic limits and explicit limitations. Development builds
also describe `mesh_inspection`, `cad_proximity`, `mesh_visualization` and
`mesh_comparison`. The development `controls` object additionally lists
`construction_efforts`, `construction_default`, `construction_stage_caps` and
`construction_caps_override_overall_resources=False`. For example:

```python
import paacs_mesh as pm

caps = pm.capabilities()
assert caps["schema"] == "paacs-mesh-capabilities-v1"
if caps.get("mesh_inspection", {}).get("all_cells_available_by_pagination"):
    print(caps["mesh_inspection"]["metrics"])
```

Capabilities describe implemented operations and restrictions; they do not
promise that arbitrary CAD will mesh. `diagnostics.face_achieved_size=False`
means native reports do not certify complete local sizing. The opt-in inspector
can still measure boundary-edge populations, with their scope stated explicitly.

## Report-only diagnostics

`paacs-mesh-diagnostics-v1` contains `status`, `source_sha256`, `units`,
`validation`, `shape_review`, `summary`, `bodies`, `timings`, `issues`, `issue_count`,
`issues_truncated` and `unavailable`. Status is `valid`, `valid_with_warnings`,
`failed` or `incomplete`. `summary.native_hint` copies an available native hint.
`unavailable` describes this report-only operation, even when a separate opt-in
inspector could supply additional evidence.

Each issue has `id`, `code`, `severity`, `scope`, `body_tags`, `cad_face_tags`,
`cell_ids`, `evidence` and `candidate_actions`. Its ID derives from source hash,
code, scope and body tags; it is not a cell ID or proof of correspondence between
different geometries. Actions are suggestions, not automatic operations.
Severity `info` includes review advisories which do not change native validity
or the diagnostic validity status. `warning` and `error` retain their existing
preferred-quality, resource and failure meanings.

| Code | Meaning |
| --- | --- |
| `selection_failed` | Selection or CAD inspection failed before meshing; candidate evidence is retained. |
| `control_request_failed` | A structured sizing request failed. |
| `operation_failed` | Another structured request failed before a native mesh report was available. |
| `mesh_refused` | Native output was refused or no mesh was published. |
| `validation_incomplete` | Published output lacks required completion evidence. |
| `native_validity_failed` | Reported audit/quality validity or invalid-cell counts fail. |
| `preferred_quality_not_met` | Preferred shape quality was missed; warning policy applies. |
| `strict_quality_not_met` | An explicit strict quality threshold was missed. |
| `near_flat_shape_review` | A severe reported Native shape minimum warrants inspection, independently of configured quality thresholds. |
| `low_shape_distribution_review` | The reported Native shape distribution warrants inspection. Intentional anisotropy and solver requirements need separate assessment. |
| `source_shell_not_accepted` | Recorded source-shell geometry attempts did not complete acceptance. Their measured volume errors and limits remain available. |
| `source_shell_validation_interrupted` | Source volume criterion passed, but cancellation or resource refusal interrupted acceptance; this does not establish a volume defect. |
| `source_fallback_spacing_unverified` | A native fallback source was privately published without a global source-sizing certificate. Inspect achieved sizes; Native validity and failure status are unchanged. |
| `initial_construction_fallback` | The original automatic construction failed and a fallback was attempted. Surface and volume diagnostics are retained separately. |
| `volume_sizing_not_applied` | The published construction explicitly reports that it did not apply a volume size field. Geometry validity is separate. |
| `construction_budget_limited` | A native construction stage reports a bounded stop. This is a warning, not a sizing certificate or a replacement for native validity checks. |
| `body_meshing_failed` | A body reports a meshing failure. |
| `assembly_interference` | Preserved parts have material overlap. |
| `assembly_interference_unresolved` | Part interference could not be established. |
| `native_report_inconsistent` | Report status, publication or assembly claims contradict other evidence. |
| `workflow_failed` | A supplied workflow did not succeed. |
| `workflow_incomplete` | Workflow diagnostics are explicitly incomplete. |

Native `reason`, `surface_reason`, `volume_reason` and selection reasons are
separate fields, not aliases for these summary issue codes. Inspect their
recorded evidence rather than matching human message text.

### Automatic shape review

`shape_review` interprets complete Native `subtet_mean_ratio` summaries (ideal
shape = 1). It recommends review when the minimum is at most `0.001`, or when
the median is at most `0.1`, or when both the median is at most `0.2` and the
fifth percentile is at most `0.05`. The exact policy is available through
`capabilities()["diagnostics"]["automatic_shape_review"]`.

These are automatic review heuristics, separate from `warning_quality` and
`minimum_quality`. Their issues have severity `info`; disabled preferred
warnings remain disabled, configured strict thresholds retain their meaning,
and native validity is unchanged. A thin or intentionally anisotropic mesh can
require review without being invalid. Suggestions inspect cell shapes, regions
and solver/element requirements; they do not automatically prescribe refinement.

The review status is `recommended`, `no_automatic_trigger` or `unavailable`.
Missing, incomplete or unknown metric conventions do not receive invented
measurements. Neither the absence of a trigger nor a Native pass is shape or
simulation approval. `solver_adequacy` remains `unassessed`, and report-only
advisories identify no individual cells.

### Source-shell and volume-sizing receipts

Native body reports add `source_shell` for the latest source strategy. Its
`available` flag distinguishes recorded attempts from missing evidence;
`accepted` concerns source geometry, not completed element sizing. An attempt
records its chord limit, work/time, node/triangle counts, source immutability
check, surface diagnostic, signed/CAD volume, relative volume error and whether
the unchanged volume criterion was met. Measurements requiring a published
surface are `null` for failed builds. At most three attempts are recorded.

During sized-front preparation, Normal admission can rebuild a published shell
which misses the requested relative volume bound at a tighter effective chord,
up to three builds.
Every attempt consumes the shared resource/cancellation budget. The global
size, minimum size, compiled field and requested volume tolerance stay fixed.
Explicit construction limits remain active; failure to admit required fidelity
is a refusal. Bounded sized-front preparation retains its historical single
build; existing fallback retry policies are unchanged.

`initial_front_fallback` retains `volume_diagnostic` and the original failed
`source_shell`, in addition to `surface_diagnostic`, work and time. A successful
surface stage can therefore coexist with a source-volume refusal without its
original reason disappearing behind `surface_diagnostic="none"`.

`volume_sizing.applied` is `true`, `false` or `null` for reported construction
use, explicitly absent use, or unavailable evidence. A kernel fan without a
volume size field reports `false`. `sizing_completion_certified` remains false:
applying a field, accepting source geometry or exhausting a queue does not
verify every final element size. Use independent mesh inspection for achieved
size distributions. Diagnostic summaries bound source-attempt lists by
`max_issues`, with totals and truncation flags.

`constraint_recovery` also records `initial_topology_used`,
`initial_topology_diagnostic`, `initial_topology_cells`,
`initial_topology_import_work` and `initial_topology_peak_modeled_bytes`.
These identify whether recovery consumed the retained initial tetrahedral
topology and report import work and modeled peak storage. They are construction
evidence, not new geometry, sizing or solver certificates. A false `used` value
does not distinguish an absent carrier from refused import; read the diagnostic.
Historical reports may omit the object or these additive fields. Diagnostic
summaries preserve available fields without inventing a successful import.

### Construction receipts (development)

See [construction effort and caps](CONSTRUCTION_BUDGETS.md) for request scope and
replay behaviour. The native report's additive `construction` object records
`effort` (`legacy`, `auto` or `thorough`) and the four requested optional caps:
`max_surface_points`, `max_surface_chord_points`, `max_interior_candidates` and
`max_interior_sites`. An automatic cap is `null`; effective limits appear in
stage receipts. `diagnose()` copies this object to `summary.construction`.
Historical reports can omit these fields.

Each native body can provide `surface_spacing`, a list of rebuilt-face receipts:

| Fields | Meaning |
| --- | --- |
| `face_tag` | Imported CAD face identity for the reported stage. |
| `point_limit`, `chord_point_limit`, `chart_node_limit` | Effective stage ceilings within overall resources. |
| `added_points`, `chord_added_points` | Points actually inserted by the two stages. |
| `pending_entries`, `blocked_candidates` | Recorded unresolved triangle entries and blocked proposals. |
| `termination`, `queue_exhausted`, `budget_limited` | Stop diagnostic and queue/allowance state. |

Pending entries can increase as triangles subdivide. They are not a monotone
error measure, a count of bad final cells or evidence of complete achieved
sizing. Queue exhaustion does not certify CAD fidelity or cell quality. Routes
that do not rebuild charts need not provide face receipts.

A route can fail during shared-edge sampling before rebuilding any chart;
an empty `surface_spacing` list therefore does not certify surface completion.
The final body's route also need not describe every earlier failed trial.
Retain stage traces and their retry costs when diagnosing a fallback. A generic
`nonconforming_edge` surface refusal can describe a numerical sampling failure;
it does not by itself establish invalid source CAD.

Normal's curve integration measures a growth-prepared metric: endpoint feature
targets extend continuously by physical distance, and local integration
stencils lower sizes under the existing internal growth slope before comparing
coarse/fine estimates. Raw sampled requests and final global growth preparation
are retained separately. This is not an exact integration-error certificate
for raw curvature, a new public per-curve receipt, or a final sizing certificate.
Internal numerical probes can record the last refused parameter interval and
raw/prepared samples; those private fields are not promised by this JSON schema.

The body's `initial_front` adds `budget_reason`, effective `candidate_limit`,
`site_limit`, `query_limit`, `work_limit`, and allocation counters
`reserved_bulk_candidates`, `reserved_bulk_sites`, `source_candidates`,
`source_deferred_faces` and `source_phase_deferred`. Reservations record effort
allocation, not guaranteed accepted bulk points. A nonfinite
`maximum_unresolved_void_ratio` is emitted as `null`, indicating unavailable
completion evidence rather than zero unresolved space. This historical
construction observation is not a final covering-radius or achieved-sizing
certificate: later front closure can change the unresolved regions. Neither a
finite ratio nor queue exhaustion certifies that sizing targets were met.
The full native body
also exposes `native_surface_work`, `native_volume_work`,
`native_surface_seconds` and `native_volume_seconds`.

Progressive initial-front storage adds these receipts under `initial_front`:

| Fields | Meaning |
| --- | --- |
| `auxiliary_imported_cells` | Cells actually installed from the initial numerical seed. This is not final output Tet4 count. |
| `auxiliary_import_work` | Charged point/cell import operations; excludes the constructor's own work. |
| `auxiliary_import_modeled_bytes` | Auxiliary bootstrap estimate while the candidate and staged imported complex coexist. |
| `auxiliary_initial_capacity` | Admitted auxiliary cell capacity after successful bootstrap. |
| `auxiliary_hard_capacity` | Current cell ceiling derived from the explicit limit and the auxiliary's current memory share. It can change as the front and heap use that shared allowance. |
| `auxiliary_growth_attempts`, `auxiliary_growths` | Attempted and completed capacity-growth operations. An attempted operation can refuse without changing live cells. |
| `auxiliary_growth_work` | Charged copies of occupied auxiliary arrays during growth, rather than a charge for empty future slots. |
| `auxiliary_stored_cells`, `auxiliary_live_cells`, `auxiliary_stored_peak` | Current stored cells, current live cells and peak stored cells; stored arrays can include retired cells until compaction. |
| `auxiliary_storage_capacity`, `auxiliary_working_bytes` | Current admitted cell capacity and modeled byte share of the auxiliary stage. Neither is the final mesh count or whole-process memory. |
| `compactions`, `retired_cells_reclaimed`, `compaction_work` | Completed compactions, reclaimed retired cells and charged compaction work. |
| `compaction_peak_modeled_bytes` | Auxiliary estimate admitting old and replacement storage together during compaction. |
| `maximum_combined_modeled_bytes` | Largest recorded initial-placement admission estimate, including the enclosing boundary allowance and the front, bulk heap and auxiliary estimates. It is not a whole-job peak-memory measurement. |
| `maximum_heap_capacity` | Peak modeled bulk-heap reservation in entry slots. During reallocation this includes old and new slots simultaneously; it is not maximum queue length. |
| `resource_retries` | Source-face or bulk-candidate requeues for auxiliary compaction/capacity recovery. These are construction retries, not independent full meshing attempts. |
| `auxiliary_resource_diagnostic` | Latest auxiliary resource stop or state diagnostic; inspect it with `budget_reason`, overall status and the other stage receipts. |

Historical or inapplicable routes can omit these fields or leave their counters
at zero. A zero does not establish that construction, sizing or validation
completed. Byte fields use the native logical model; they do not measure
operating-system resident memory or all native-kernel allocations. Capacity
growth and compaction preserve geometric/certification checks. Completed work
returned by an independent cancelled stage is included before its cancellation
is propagated; incremental work callbacks check cancellation before another
operation.

Examples of auxiliary stop diagnostics are `bootstrap_fixed_memory`,
`bootstrap_construction_budget`, `bootstrap_import_memory`,
`bootstrap_actual_import_work`, `auxiliary_cell_limit` and
`auxiliary_memory_limit`. Initial-front `budget_reason` can also report
`combined_working_memory` or `auxiliary_resource_stalled`. These are native
stage strings, not replacements for summary issue codes. A stop identifies an
admission or progress limit; it does not prove that increasing it will resolve
the geometry. Preserve unknown additive diagnostic strings when displaying
reports.

`construction_budget_limited` carries `surface_stages`, an available
`initial_stage`, `sizing_completion_certified=False` and
`native_validity_gates_changed=False` in its evidence. Body/face identities come
from those receipts; no cell identity or cause of poor quality is inferred.
Diagnostics preserve native refusal status and do not automatically retry.
Diagnostic summaries bound nested face-stage and identity lists by
`max_issues`, with explicit totals and truncation flags. Consult the full native
report for the complete receipts; an omitted identity is not evidence that its
stage completed.

## Native resource intent

Development native mesh reports add `resources`, copied by `diagnose()` to
`summary.resources` without re-certifying the mesh:

| Field | Meaning |
| --- | --- |
| `policy` | Requested `normal` or `bounded` admission policy. |
| `requested_max_work` | Explicit requested ceiling, or `null` when omitted. |
| `effective_max_work` | Finite internal logical-work ceiling used by native admission. |
| `logical_work_cap_disabled` | True only for normal policy with omitted work; a finite overflow guard remains. |
| `max_working_bytes` | Effective modeled working allowance, not a hard process RSS cap. |
| `scope` | Explanatory accounting restrictions in Python reports; the native CLI emits the five policy/limit fields above. |

Explicit positive work ceilings apply to both policies. Historical reports can
omit this object; diagnostics then leave the summary field unavailable rather
than inventing a policy or completion certificate. Capability discovery adds
`resources.policies`, `default`, `explicit_limits_override_policy`,
`logical_work_is_cpu_instructions` and `modeled_memory_is_OS_hard_limit`.

New development Python meshing calls, new recipes and direct STEP CLI requests
default to `normal`. Saved recipes retain their policy; missing-policy older
recipes and the C++ options aggregate retain literal `bounded` defaults.
Normal with omitted work retains the native 512 MiB modeled allowance unless
the caller changes `max_working_bytes`. Explicit work, node, cell, output,
validation, cancellation and cooperative time limits still apply. Normal
admission does not select experimental `construction_effort="auto"` or
`"thorough"`, relax validity gates or guarantee arbitrary CAD compatibility.
The default effort remains `legacy`. These are resource-admission semantics,
not quality, requested-sizing, runtime or solver-accuracy guarantees.

## Mesh exploration contracts

See [mesh inspection](MESH_INSPECTION.md) for complete examples and metric
conventions. These are copied JSON objects; unavailable measurements are `null`,
not invented zeros. Nonfinite values are not emitted in JSON.

| Schema | Main evidence |
| --- | --- |
| `paacs-mesh-cell-page-v1` | Mesh identity, metric conventions, exact filtered `total`, `offset`, `limit`, `next_offset`, filters/order, optional CAD proximity header and `cells`. |
| `paacs-mesh-quality-summary-v1` | Whole-snapshot counts, metric distributions, requested/effective/achieved sizing, spatial bins and native timing/validation evidence. |
| `paacs-mesh-inspection-bundle-v1` | Exact archive/report hashes, source identity, units and native-ID schemes for the saved inspection bundle. |
| `paacs-mesh-cad-proximity-v1` | Matching source/kernel/units, complete selected-face centroid samples, evaluation count, limits and explicit non-conclusions. |
| `paacs-mesh-visualization-v1` | Embedded SVG evidence: original coordinates, selected cell IDs/values, identities, camera and cut rule. |
| `paacs-mesh-visualization-receipt-v1` | Output path/hash, selected IDs/counts, camera, cut and evidence identity. |
| `paacs-mesh-comparison-v1` | Two reported-valid snapshots, count/quality/sizing/CAD-error/timing comparisons and common-world spatial regions. |

Cell IDs are zero-based native archive cell indices, including intervening
surface cells. They remain stable for one mesh hash; different mesh hashes have
no implied cell-to-cell correspondence. Bare archive v1 files carry neither
units nor source hashes nor validation reports. Saved paired reports are
hash-checked, not re-audited or cryptographically authenticated.
CAD body-error deltas additionally require matching known kernel versions.
When either kernel is unknown or differs, measurements remain available as
separate unpaired body lists; population comparisons can still proceed.

## Refusal handling

`InspectionError` exposes `stage="mesh_inspection"`, `reason`, and `report`.
Resource reasons include `memory_budget`, `count_budget`, `work_budget`,
`time_budget`, `output_budget`, `association_budget`, `triangle_budget` and
`highlight_budget`. CAD distance adds `cad_work_budget`, `cad_time_budget`,
`cad_memory_budget` and `cad_kernel_allocation`. These refuse the operation
without silently dropping matching cells.

Identity/provenance failures include `bundle_hash`, `report_unavailable`,
`cad_identity_unavailable`, `cad_identity_mismatch`, `source_mismatch`,
`comparison_source`, `comparison_units` and `comparison_validation`. CAD query
completion failures include `cad_query_failed`, `cad_query_identity` and
`cad_query_receipt`. Archive structure and unsupported display/metric failures
have distinct reasons; they do not imply a native volume construction attempt.

CLI JSON errors use `paacs-mesh-cli-error-v1`, with `status`, `stage`, `reason`,
`message` and available bounded `report` evidence. Output is capped at 8 MiB;
oversized nested error reports are explicitly omitted. A completed operation
returns 0, an inspection/native refusal returns 1, and invalid ordinary command
arguments return 2. No inspection command retries meshing or relaxes validity.

Budgets bound admitted data, traversal and output; deadlines are cooperative.
They do not impose an operating-system RSS or CPU-time limit. A single native
CAD solve and its internal workspace cannot be interrupted by these counters.

## Planar-core receipts

Development Python and native CLI body reports add `planar`. Its `attempted`
flag is set immediately before the geometry-core call; skipped calls have a
null reason. An attempted call retains its outcome, analysis status/work,
stage seconds and up to 16 body rows. Rows include the resolved strategy and
options, count estimates, engine status/reason/phase codes, conversion outcome,
and available boundary/surface diagnostics. Estimates include surface cells;
they do not establish that the proposed strategy generated an output.

Engine codes use the selected strategy's enum definitions. Default zero
codes or timers do not prove that a substage was skipped. In particular,
analysis time is currently retained only after analysis returns, so a
cancellation during analysis can leave its timer at zero. Generation failures
retain their recorded body time. Use outer Native `seconds.generation` for
the complete generation-attempt cost and `seconds.meshing_total` for the full
native operation; nested times must not be added as independent costs.

`planar` describes core strategy/construction evidence. It does not replace
the separate original-CAD validity checks, certify achieved sizing or establish
solver adequacy.

### Method resolution and geometry fidelity (development)

Native reports add `fidelity`: nullable `requested_max_chord_error` and
`requested_max_relative_volume_error`, plus the effective relative-volume limit.
Null preserves the existing Native default. Body `fidelity` records
`effective_max_chord_error` and `chord_policy`; an unavailable limit is null.
Primitive radius defaults, source-preparation chords and exact planar geometry
remain distinguishable. These are acceptance limits, separate from observed CAD
errors, requested element size and independently measured achieved lengths.

`route_resolution` records `reason`, `candidate`, `preflight_attempted`,
`whole_body_verified`, preflight work/time, an optional retained fallback chord,
and at most two construction `attempts`. Each attempt records its method,
Native reason, publication flag, work and elapsed seconds. Failed attempts remain
present. Normal AUTO can select an existing whole-body primitive, affine or sweep
capability before general initial-point construction when controls are compatible.
Explicit initial-front choices and historical Bounded routing remain authoritative.

Capability recognition does not certify the final mesh. Native geometry,
orientation, connectivity, body and quality gates still apply after construction.
If a verified primitive attempt fails and general construction is tried, its chord
limit is retained. Top-level fidelity still records the original user request;
the fallback's source attempts report their actual effective chord. Logical
preflight and attempt work consumes the same overall allowance. Outer timings
include all attempts; nested stages and wrappers must not be counted twice.

Body `initial_topology_quality` reports the optional numerical-topology stage.
`complete` means that stage finished its bounded operation, not that every poor
cell disappeared or that the final volume mesh is complete. Read the remaining
candidate counts, diagnostics and final measured quality. `peak_modeled_bytes`
is helper-local and includes its retained input owner; it excludes the outer
held mesh envelope and is not process RSS. `delaunay_certified=false` remains
explicit. Native final geometry and validity gates still apply independently.

Quality interpretation requires a reported nonempty volume population. Empty
zero-filled summaries on a cancelled/unpublished operation do not trigger shape,
preferred-quality or strict-quality measurements; shape review remains
`unavailable` with reason `no_reported_volume_population`. The Native failure
status and requested requirements remain intact.

Python body reports also preserve the existing partition-stage evidence exposed
by the native CLI: `native_partition_seconds`, `native_assembly_seconds`, policy,
reason, attempts, piece/interface counts, logical work, volume-change admission,
trim-witness measurements and `partition_attempt_reports`. Each attempt keeps
its policy, diagnostic, failed internal piece tag, piece/interface counts, work,
seconds and surface deflection. Internal piece tags do not become original CAD
body identities merely by appearing here.

`diagnose()` limits the attempt list to `max_issues` rows per body and records
`partition_attempt_reports_total` and `partition_attempt_reports_truncated`.
The original full report is unchanged. Partition/assembly times are nested
construction evidence; they must not be added to outer generation/total times
as independent elapsed intervals. Zero/default fields do not prove stage entry,
and these receipts neither redo Native audits nor certify sizing completion.

### Native source-provider fallback (development)

`requested_native_source_fallback` is the strict Boolean request (default false),
also retained in `effective_policy`. Each source shell reports
`native_fallback_requested`, `native_fallback_attempted` and
`native_fallback_diagnostic`. Request and actual entry remain separate.
The experimental fallback is considered only after a remaining named geometric
chart/wire refusal; ordinary supported curve sampling is attempted first.

Source attempts add `provider`, `global_sizing_applied`,
`curvature_points_per_turn`, `angular_deflection` (radians, or null when
unavailable), and `witness`. A witness contains `stage`, `detail`, nullable
`edge_tag`, `face_tag` and `use_id`, plus `geometric_refusal`. Tags identify the
guarded source; an unavailable tag is null, and `use_id` may be zero. These are
construction witnesses, not repair instructions or quality-causation claims.

The same bounded attempt list counts chart and native-provider work together:
Normal admits at most three source attempts, Bounded explicit opt-in admits two,
and the default Bounded policy retains one. Existing failure, volume and source
immutability fields keep their meaning. All attempts consume cumulative work;
cancellation, budget, allocation, source-change and internal failures cannot
become a successful source switch.

Native provider angle/curvature and global volume sizing do not establish a
source edge-size certificate. `source_fallback_spacing_unverified` is info only
when the experimental fallback shell was accepted with
`global_sizing_applied=false`; it
remains available on a subsequently refused volume attempt. The summary never
uses it to recertify a mesh or hide a Native error. `max_issues` bounds retained
source-attempt lists and diagnostic evidence without changing the full report.

`requested_boundary_first` is another strict Boolean request, default false,
retained in `effective_policy` and diagnostic summaries. It is independent of
the source-provider switch. In a participating live spatial Tet4 route, true
explicitly selects boundary recovery before interior insertion, including when
the caller also explicitly enabled the initial front. The initial report then
uses `diagnostic="deferred_until_boundary_recovery"`; subsequent work is in
`interior_front`. Unavailable or deferred initial work must not be interpreted
as successful interior construction. Interior-front `nodes` and `candidates`
record attempted or privately staged work, including work retained in an
interrupted report. They do not count published nodes or certify a partial
mesh. `interior_front.work` includes charged CAD callbacks; `queries` is the
change in the body's CAD admission counter, excluding cache hits. These differ
from helper-local topology work and callback invocation counts. Preflight and
interrupted callback time remain included in this stage's elapsed time.
Both switches retain all Native validity checks, and neither certifies
requested sizing or solver adequacy.

### Material classification progress (development)

Body `material_classification` keeps cumulative returned-batch evidence:
`attempts`, `requested_points`, `attempted_points`, `native_queries`,
`strict_reuses`, `ray_admissions`, `work` and `seconds`. `diagnostic`, `complete`
and `source_model_verified_unchanged` describe the last batch. The scope is
`cumulative material classification attempts; last batch status; not final publication`.
Interrupted calls retain attempted work before cancellation propagates;
their partial admission masks are not published. Repeated batches may request
the same point, so these are not distinct node or cell counts.

`attempted_points` partitions into ordinary native queries, strict cached
reuses and positive original-CAD ray admissions. It must not exceed cumulative
requested points. An ordinary query is counted when it starts; interruption
may therefore leave an attempted point without a classification decision.
Classification completion does not mean every point was
inside: completed masks can contain outside decisions that drive cavity repair.
Neither classification completion nor a source guard certifies a final mesh.

Strict serial boundary-first batches additionally report `profiled_batches`,
`ray_witness_attempts`, `classifier_setup_seconds`, `native_query_seconds` and
`ray_witness_seconds`. These accumulate only profiled returned batches;
completion/source guards still describe the last batch, not the cumulative
profile. `profiled_batches=0` means the phase cost was not measured, not that it
was zero. Supplied phase times are finite and nonnegative; their sum is at most
material elapsed `seconds` within floating-point tolerance. Unprofiled calls
can still contribute to the broader material counts and elapsed time.

The explicit boundary-first mode can use positive ray witnesses against the
actual original trimmed CAD in strict serial classification. Ambiguous and
unadmitted points use the ordinary native classifier. This does not relax CAD
tolerances, allow a distance band, or replace final validity checks. Actual
`ray_admissions` record its use; a requested switch alone is not evidence of it.

Python also preserves the existing body `native_candidate_queries` and
`native_candidate_ray_admissions` counters exposed by the native CLI. Material
times/counters are nested within volume work and can overlap these broader
counters. Do not add them to generation, total time or other nested counts as
independent costs. Diagnostic summaries validate counter types and coherence,
retain the supplied progress, and do not rerun geometry classification.
