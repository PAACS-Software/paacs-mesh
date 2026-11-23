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
`validation`, `summary`, `bodies`, `timings`, `issues`, `issue_count`,
`issues_truncated` and `unavailable`. Status is `valid`, `valid_with_warnings`,
`failed` or `incomplete`. `summary.native_hint` copies an available native hint.
`unavailable` describes this report-only operation, even when a separate opt-in
inspector could supply additional evidence.

Each issue has `id`, `code`, `severity`, `scope`, `body_tags`, `cad_face_tags`,
`cell_ids`, `evidence` and `candidate_actions`. Its ID derives from source hash,
code, scope and body tags; it is not a cell ID or proof of correspondence between
different geometries. Actions are suggestions, not automatic operations.

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

The body's `initial_front` adds `budget_reason`, effective `candidate_limit`,
`site_limit`, `query_limit`, `work_limit`, and allocation counters
`reserved_bulk_candidates`, `reserved_bulk_sites`, `source_candidates`,
`source_deferred_faces` and `source_phase_deferred`. Reservations record effort
allocation, not guaranteed accepted bulk points. A nonfinite
`maximum_unresolved_void_ratio` is emitted as `null`, indicating unavailable
completion evidence rather than zero unresolved space. The full native body
also exposes `native_surface_work`, `native_volume_work`,
`native_surface_seconds` and `native_volume_seconds`.

`construction_budget_limited` carries `surface_stages`, an available
`initial_stage`, `sizing_completion_certified=False` and
`native_validity_gates_changed=False` in its evidence. Body/face identities come
from those receipts; no cell identity or cause of poor quality is inferred.
Diagnostics preserve native refusal status and do not automatically retry.
Diagnostic summaries bound nested face-stage and identity lists by
`max_issues`, with explicit totals and truncation flags. Consult the full native
report for the complete receipts; an omitted identity is not evidence that its
stage completed.

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
