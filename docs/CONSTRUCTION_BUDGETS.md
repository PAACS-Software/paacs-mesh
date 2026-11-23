# Construction effort and caps

These are unreleased development controls requiring rebuilt native bindings.
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

`auto` estimates allowances from boundary density, requested sizes and available
resources. `thorough` requests a larger estimate, still bounded by those
resources. Both reserve construction effort for bulk placement before resuming
deferred wall-front work. A reservation does not guarantee that candidates will
be accepted: the same geometric and cell-shape checks still apply. Neither mode
adds a volume sizing region.

Leaving all caps unset with `legacy` retains the established construction path.
Supplying an explicit cap also resolves its supporting query, storage and work
allowances within the global budgets. Consequently, explicitly supplying an old
point limit need not produce the same mesh as leaving every cap unset.

`legacy` here names a construction-allowance policy. It does not turn off the
advancing interior front. Use the separate `advancing_interior_front` and
`recovery` options to request their existing behaviour.

Historical allocation is preserved by `legacy` with all four caps unset.
Setting a cap explicitly also scales its supporting query, work and storage
allowances within the overall limits. Even an explicit historical cap value
can therefore produce different points.

## Saved recipes

`create_recipe()` accepts the same five options. Save/load and every permitted
retry preserve them. A finer retry changes the requested global size according
to the saved policy; it does not silently increase construction caps or global
resource budgets.

Older recipe-v1 files lack these fields. Each missing effort field resolves to
the literal historical value `legacy`; each missing cap resolves to `None`.
This normalization is independent of any future default for new requests.
Source, kernel, engine, selection and other recipe guards remain active.

## Command line

```sh
paacs-mesh part.step out.vtk --size 0.5 --construction-effort auto --max-surface-points 400 --max-surface-chord-points 128 --max-interior-candidates 20000 --max-interior-sites 256 --max-work 100000000 --max-seconds 90
```

The direct positional STEP command accepts these flags. `paacs-mesh run`
replays the values from its saved recipe and accepts no ad hoc construction
override flags.

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
