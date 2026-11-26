# Construction effort and caps

These controls are part of the RC4 source candidate and require matching native
bindings. The matching release uses the complete qualification gates; published
RC3 wheels retain their original interface.
The default remains `construction_effort="legacy"` while the alternative
policies are qualified against geometry, quality and runtime tests.

The [measured budget study](CONSTRUCTION_BUDGET_STUDY.md) records both
improvements and unresolved quality/resource limitations behind that default.

Construction allowances are separate from requested element sizes, recovery
mode and validity requirements. Increasing an allowance permits more bounded
construction work; it does not guarantee better quality, completion or a
particular runtime. Existing global work, memory, node, output and time budgets
still apply. No geometry or validity check is disabled.

```python
mesh = model.mesh(
    .5, assembly="conforming", workers=1, validation_workers=1,
    construction_effort="auto",
    max_surface_points=400,
    max_surface_chord_points=128,
    max_interior_candidates=20_000,
    max_interior_sites=256,
    max_work=100_000_000, max_seconds=90,
)
```

| Option | Meaning |
| --- | --- |
| `construction_effort` | `legacy`, `auto` or `thorough`. The initial default is `legacy`, preserving established construction allowances. The alternatives request adaptive bounded allowances. |
| `max_surface_points` | Optional cap on surface spacing point insertion per rebuilt face chart. |
| `max_surface_chord_points` | Optional cap on subsequent chord-refinement point insertion per rebuilt face chart. |
| `max_interior_candidates` | Optional cap on initial interior candidate proposals. |
| `max_interior_sites` | Optional cap on admitted initial interior sites. |

Each cap is a positive integer or `None`. Booleans, zero, negative values,
floating-point values and strings are rejected by the Python API. `None`
leaves the allowance to the selected policy and remaining global resources.
Reported effective limits can be lower than an explicit cap.

These controls currently affect chart rebuilding and initial interior placement
in the compatible single-solid native Tet4 route. They do not extend assembly,
structured-constructor or kernel-tessellation support. A cap is a ceiling on
added points or proposals, not a requested node count. Existing boundary points
are not counted as new spacing points.

With bounded admission, `auto` estimates allowances from boundary density,
requested sizes and available resources. `thorough` requests a larger estimate,
still bounded by those resources. Both reserve construction effort for bulk placement before resuming
deferred wall-front work. A reservation does not guarantee that candidates will
be accepted: the same geometric and cell-shape checks still apply. Neither mode
adds a volume sizing region. With normal admission, unset stage caps instead
derive their ceilings from the remaining overall resources. Explicit stage
caps remain authoritative under either resource policy.

Leaving all caps unset with `legacy` and `resource_policy="bounded"` retains
the established construction path.
Supplying an explicit cap also resolves its supporting query, storage and work
allowances within the global budgets. Consequently, explicitly supplying an old
point limit need not produce the same mesh as leaving every cap unset.

`legacy` here names a construction-allowance policy. It does not turn off the
advancing interior front. Use the separate `advancing_interior_front` and
`recovery` options to request their existing behaviour.

## Resource policy

`resource_policy` is independent of construction effort and recovery. New
Python meshing calls, newly created recipes and direct STEP CLI requests in
RC4 default to `normal`. Explicit `bounded` retains
historical admission. Saved recipes use their recorded policy; an older recipe
without this field retains literal `bounded` semantics.

| Policy | Omitted `max_work` | Explicit positive `max_work` |
| --- | --- | --- |
| `bounded` | Historical 50 million logical work allowance. | The supplied ceiling applies. |
| `normal` | No practical logical-work ceiling; a finite internal overflow guard remains. | The supplied ceiling applies. |

Normal admission uses the same validity requirements. It does not promote the
experimental `auto` or `thorough` construction effort policies. Other resource
allowances remain independent: omitted `max_working_bytes` retains the native
512 MiB modeled allowance. Byte admission is not a hard operating-system RSS
limit. Explicit node, cell, output, audit and cooperative time limits still
apply. Work units count implementation-specific charged operations, not CPU
instructions, elapsed seconds or elements.

```python
mesh = model.mesh(.5, resource_policy="normal", max_working_bytes=1024**3)
mesh = model.mesh(.5, resource_policy="normal", max_work=100_000_000,
                  max_working_bytes=1024**3, max_seconds=90)
```

The second request retains a finite work ceiling. Neither policy nor a larger
allowance is a certificate of achieved sizing, preferred quality or solver
adequacy. Inspect the actual stage receipts and independent measurements.

### Native admission under normal policy

Normal admission separates a future capacity ceiling from storage needed now.
The initial interior front constructs its auxiliary Delaunay seed once, then
admits the cells actually returned. It does not reserve or charge for importing
every cell the request might eventually permit. Live storage grows between
candidate batches, with headroom for current cells and bounded cavity work.
Growth and compaction account for old and new arrays while both coexist;
failed admission leaves the live topology intact. The front and bulk-candidate
heap share the remaining modeled allowance with this auxiliary mesh. Work
receipts include completed constructor attempts, imports, retries and storage
copies, including work returned by a cancelled stage.

Stages that do not coexist receive separate allowances. In particular, initial
point placement retains the actual boundary and its related data without
reserving a future volume front or apex cache. Later volume and recovery
stages retain their own input and workspace envelopes. Normal curve sampling
adapts its numerical integration to the curve and requested spacing; its
segment and face-rebuilding ceilings come from model counts, index ranges and
remaining node, work and memory limits. Native shared-edge, trim, orientation,
chord and volume checks still apply. Increasing capacity does not remove
geometric refusal or bounded local recovery/stall checks.

Normal curve quadrature prepares the sizing metric before estimating its
numerical error. Incident-edge feature sizes influence a continuous endpoint
neighbourhood through the existing internal slope `sqrt(1.8)-1`, using physical
distance to each endpoint. Each five-point integration stencil then applies
forward/backward growth bounds to copies of its requested sizes. Both steps
can only lower sizes. Euclidean distance is no greater than arclength, so these
local envelopes are conservative relative to the corresponding arclength
growth bound. Raw cached requests and the final global arclength growth pass
remain separate.

This preparation prevents isolated endpoint or curvature-size jumps from
forcing relative-error subdivision down to machine precision. It deliberately
integrates a prepared sizing metric, rather than claiming an error certificate
for a discontinuous raw curvature field. The internal slope is a construction
rule, not a certificate that the user's requested growth or every final edge
size has been achieved. Sampling work remains charged, Bounded sampling retains
its historical path, and no permissive numerical fallback or relaxed CAD check
is introduced. The [resource study](RESOURCE_POLICY_STUDY.md) separates the
curve-level regression evidence from full volume-mesh qualification.

These are logical stage-admission estimates, not a measurement or hard bound
on whole-process resident memory. Native kernel workspaces, allocator overhead
and other process state can differ from the model. A larger allowance also
cannot make an unsupported geometric construction succeed. Requests still
need measured CAD validity, element quality and achieved sizing; an admitted
mesh can miss a requested sizing target without becoming an invalid mesh.

The C++ aggregate `NativeCadMeshOptions{}` retains historical bounded defaults.
Use `normalNativeCadMeshOptions()` to opt into normal admission with its work
allowances, then assign any explicit ceilings. Setting only
`construction.resourcePolicy` on an existing options object does not erase
work limits already stored in that object. Python and the CLI resolve omitted
work from the selected policy. `auto` and `thorough` remain separate opt-ins.

## Saved recipes

`create_recipe()` accepts the same construction options and `resource_policy`.
Save/load and every permitted
retry preserve them. A finer retry changes the requested global size according
to the saved policy; it does not silently increase construction caps or global
resource budgets.

Older recipe-v1 files lack these fields. Each missing effort field resolves to
the literal historical value `legacy`; each missing cap resolves to `None`.
This normalization is independent of the default for new requests.
Source, kernel, engine, selection and other recipe guards remain active.
An old recipe lacking `resource_policy` resolves to literal `bounded`, although
new requests default to `normal`. An omitted work value remains
`None` in the recipe and is resolved by that saved policy. A recipe's workflow
deadline remains separate from native resource policy.

## Command line

```sh
paacs-mesh part.step out.vtk --size 0.5 --construction-effort auto --max-surface-points 400 --max-surface-chord-points 128 --max-interior-candidates 20000 --max-interior-sites 256 --max-work 100000000 --max-seconds 90
```

The direct positional STEP command accepts these flags. `paacs-mesh run`
replays the values from its saved recipe and accepts no ad hoc construction
override flags.

The direct command also accepts `--resource-policy normal|bounded`,
`--max-work`, `--max-working-bytes`, `--max-nodes` and `--max-cells`. Saved
`run` commands replay resource intent from the recipe; they do not accept
resource override flags. An explicit positive work cap wins regardless of
policy/cap flag order. Invalid work-cap tokens do not replace a previously valid
configuration in a native CLI error receipt.

## Diagnostics

The native report's `construction` block records the requested effort and
optional caps. Each body's `surface_spacing` records actual limits, inserted
points, pending/blocked work and termination diagnostics. The `initial_front`
report records actual candidate/site/query/work limits, budget reasons and
source/bulk allocation counters. Inspect those outcomes together with achieved
sizes and element quality. A configured allowance is not evidence that a stage
completed its targets.

`pending_entries` counts unresolved triangle entries in the construction stages.
It can increase when a larger allowance subdivides triangles, so it is not a
monotone error metric. `queue_exhausted` reports queue state; it does not certify
achieved sizes, CAD fidelity, cell quality or simulation accuracy. A blocked
candidate can indicate a geometric or construction limitation that a larger
budget will not resolve.

`maximum_unresolved_void_ratio` is a historical construction observation, not
a final covering-radius or achieved-size certificate. Later front closure can
change which regions remain unresolved. A finite value, zero or an exhausted
queue must not be interpreted as proof that every requested size was met.
The progressive storage receipts are documented in
[diagnostic contracts](DIAGNOSTIC_SCHEMAS.md#construction-receipts-rc4).

`diagnose()` emits `construction_budget_limited` when the native report records
a limited stage, with its actual body and available face tags. This warning is
separate from native validity failures. Inspect the recorded stop reason and
measured quality/achieved sizes before changing limits. Global refusals remain
refusals; raising a stage allowance never bypasses overall resource limits.
Diagnostic summaries can truncate nested face-stage and identity lists at their
configured issue limit, reporting totals and truncation explicitly. The full
native report retains the complete construction receipts.

The C++ counterparts are `NativeCadMeshOptions::construction` and
`CadConstructionBudget`. Their zero-valued stage fields select the policy's
automatic allowance; zero does not mean unlimited work. The Python API uses
`None` for that choice and rejects explicit zero.
