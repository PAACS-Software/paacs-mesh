# Inspection and diagnostics

These discovery commands and diagnostic helpers are available starting with
`1.2.0rc3`. The earlier `1.2.0rc2` release predates them; installing that release
alone does not add these commands. See the [Python guide](PYTHON.md) for installation.

PAACS Mesh provides the same inspection, recipe and diagnostic operations to
Python programs, shell users and agents. Discover the supported operations,
inspect the CAD model, make an explicit selection, then evaluate the native
result. The interface does not autonomously choose a simulation resolution.

STEP coordinates, target sizes and geometric errors use millimetres; reported
volumes use cubic millimetres. Face and body tags belong to a particular imported
model. Save the source SHA-256 and guarded selections with a recipe instead of
assuming that a face number identifies the same feature in another STEP file.

## Command line

```sh
paacs-mesh capabilities --json
paacs-mesh inspect part.step --kind faces --limit 100 --topology --json
paacs-mesh run part.recipe.json --source part.step --output part-run-001 --json
paacs-mesh diagnose part-run-001/report.json --json
```

The simple conversion command remains available:

```sh
paacs-mesh part.step part.vtk --size 0.7
```

Omit `--json` on the subcommands for a readable summary. JSON mode writes one
JSON document to standard output; native transfer messages and other logs go to
standard error. Consumers should parse the complete standard-output document,
not look for a JSON fragment among progress messages.

Inspection is paginated. `--offset` selects the next page, and `--limit` must be
between 1 and 1000. Follow `next_offset` until it is `null` rather than assuming
the first page contains every face. The inspection schema is
`paacs-cad-inspection-v1`; its entity records include geometric properties and
body ownership. `--topology` includes connectivity information. A geometric
selection can be checked visually with the Python preview API below.

`run` replays the source-bound recipe without changing its controls or retry
policy. The source must be supplied explicitly. The output is a new local
directory containing the STEP snapshot, intent, recipe, tagged mesh exports and
full attempt report. Existing destinations are refused. This local publication
does not upload anything to a service.

Exit codes are 0 for a completed successful operation (including a valid mesh
with warnings), 1 for a failed or incomplete native/workflow result, and 2 for
invalid arguments or malformed input metadata. Human and JSON presentations use
the same diagnostic findings.

## Select by geometry and save the decision

```python
import paacs_mesh as pm

model = pm.read_step("part.step")
page = model.inspect("faces", limit=100, topology=True)
print(page)

# Match the intended feature and require exactly one candidate.
wall = model.select_faces(surface="cylinder", radius=1.0, expect=1)
model.preview(wall).write_svg("selected-wall.svg")
model.define_group("bore_wall", wall, group_id="part-bore-wall")

recipe = pm.create_recipe(
    model, .45,
    selected_face_sizes=[(wall, .22)],
    assembly="conforming", min_size=.18, growth=1.3,
    workers=1, validation_workers=1,
    max_work=100_000_000, max_attempts=1, max_seconds=90,
    warning_quality=.25,
)
pm.save_recipe(model, recipe, "part.recipe.json")
result = pm.run_saved_recipe("part.recipe.json", "part.step")
diagnostic = pm.diagnose(result.mesh)
print(pm.format_diagnostics(diagnostic))
pm.publish_result(result, "part-run-001")
```

This illustrative selection requires an input containing the specified
cylindrical face. Ambiguous or empty selections raise `SelectionError` and
include the expected count and candidate tags. Review those candidates or a
labelled preview before narrowing the query. Predicates can also use body, area,
normal, centroid region, adjacency and saved names/groups. Property guards and
the source hash protect saved selections when a recipe is replayed.

Spatial face and volume controls currently require one supported Tet4 solid
with `assembly="conforming"`. Multiple independent parts retain their body
identities under the default `preserve_parts` policy. Material interference is
reported; overlapping part meshes do not become one non-overlapping domain.

Rebuilt development wheels provide opt-in [construction effort and stage
caps](CONSTRUCTION_BUDGETS.md) for the compatible single-solid Tet4 route.
`construction_effort="legacy"` remains the default; `auto` and `thorough` can
allocate more effort and reserve a share for bulk interior placement. Optional
positive integer caps control surface points per rebuilt face chart, chord
points, interior proposals and admitted sites. The overall work, memory, node,
output and time budgets still apply. This is separate from sizing controls and
recovery policy, and does not extend assembly support.

Direct requests and `create_recipe()` accept the same construction options.
Saved recipes preserve those options through replay and retries; older recipes
without them resolve to `legacy` with automatic stage caps. A bounded retry
does not silently increase them. Use `capabilities()["controls"]` to discover
the effort names, default and cap names supported by the installed build.

## What a diagnostic means

The development interface adds [per-element mesh inspection](MESH_INSPECTION.md)
and [explicit Python request signatures](API_USABILITY.md). Those operations
are opt-in and separate from the report-only `diagnose()` operation below.
See the [schema and issue-code reference](DIAGNOSTIC_SCHEMAS.md) for machine
contracts. Published RC3 wheels retain their original interface.

`pm.diagnose(mesh_or_report, max_issues=32)` accepts a native mesh or its report,
or a bounded workflow report. It returns copied JSON data using
`paacs-mesh-diagnostics-v1`. `pm.format_diagnostics(diagnostic)` presents that
same data as text. Neither operation independently reruns or certifies the
underlying mesh checks.

The diagnostic contains:

| Field | Meaning |
| --- | --- |
| `status` | `valid`, `valid_with_warnings`, `failed` or `incomplete` |
| `source_sha256`, `units` | The reported CAD source identity and millimetre convention |
| `validation` | Native audit/quality completion and validity, mesh scope and interbody certification |
| `summary` | Available node/cell counts and aggregate native quality measurements |
| `bodies` | Available body-level routes, CAD deviation and volume measurements |
| `timings` | Reported timings, including workflow attempts when supplied |
| `issues` | A bounded list of measured findings and candidate next actions |
| `issue_count`, `issues_truncated` | Whether the displayed issue list omits findings |
| `unavailable` | Evidence the current report cannot supply |

Each issue has an `id`, a machine-readable `code`, a severity and a scope. Its
`body_tags`, `cad_face_tags` and `cell_ids` contain only identities supported by
the evidence. Empty identity lists mean that localization is unavailable. The
`evidence` object holds the measurements behind the finding;
`candidate_actions` are suggestions, not operations already performed or
guaranteed remedies. Use the source identity and issue ID together when
tracking a finding. Entity and cell identities are not persistent identifiers
across changed source files or remeshing.

Quality warnings remain distinct from invalidity. A low mean-ratio value can
produce a warning while all mandatory native checks pass. Incomplete checks
cannot establish a valid result. For assemblies, read the mesh scope and
interference information alongside the status: independently valid overlapping
parts are unsuitable as one non-overlapping material domain without further
modelling decisions.

Development reports can also produce `construction_budget_limited`, with the
limited stage's body and available face tags. Review its stop reason and actual
size/quality measurements before increasing limits. A blocked candidate need
not be fixed by more effort. Pending triangle-entry counts can rise under
subdivision; neither those counts nor queue exhaustion certify sizing, CAD
fidelity or element quality. Existing native failures retain their refusal
status regardless of this warning.

Current native reports provide aggregate quality quantiles and body-level CAD
errors. They do not identify which individual cells caused a global quality
minimum, associate that minimum with a specific bore face, measure achieved
local size/growth everywhere, or estimate a solver's discretization error.
Diagnostics state those limitations instead of inventing a location. A
candidate sizing change should therefore be measured and compared; smaller
elements alone do not establish better shape quality or simulation accuracy.

## Explicit retries and comparisons

The default recipe makes one attempt. To authorize bounded alternatives, set
`allow_advanced=True`, an explicit `finer_factors` sequence, `max_attempts` and a
total `max_seconds` budget when creating the recipe. The existing workflow
selects eligible alternatives from native failure diagnostics. Resource limits,
cancellation, allocation failures, unknown failures and unsupported controls
are not invitations to retry indefinitely. A quality warning alone does not
trigger another attempt.

Reports retain each attempt's exact options, trigger, native report and elapsed
time. Time spent in failed attempts counts toward the workflow total; export
and local bundle publication are recorded separately. Deadlines are
cooperative, and CAD calls or queue waits can overrun them. Metadata and report
inputs are bounded to 8 MiB; diagnostic output is also bounded by `max_issues`.

Run the complete, standard-library-only example from a checkout:

```sh
python examples/python/inspect_refine.py --output out/inspect-refine
```

It uses the project-owned `tests/fixtures/cad/cylinder.step`, inspects faces,
selects its cylindrical wall by geometry, writes a labelled SVG, and compares a
baseline recipe with a deliberately finer wall recipe. It reopens both recipes,
checks the source identity, mandatory native audits and face/body tags, then
saves two atomic bundles and `comparison.json`. Review the measured quality,
CAD errors and timings before selecting either result. The example does not
require NumPy, a GUI or an external solver.
