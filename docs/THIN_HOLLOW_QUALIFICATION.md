# Thin hollow extrusion qualification

The owned circular, square and rectangular STEP fixtures isolate wall
thickness, length and placement from assembly contact and imported modelling
history. A rounded rectangular profile adds curved-to-planar transitions.
The generator and serial Python driver are in
[benchmarks/thin_hollow](../benchmarks/thin_hollow/README.md).

The matrix has 22 models: 1 mm and 4 mm walls, 40 mm and 160 mm lengths,
rotated/translated thin variants, and four rounded rectangular variants.
Construction and STEP reopening check one valid solid, analytic volume,
millimetre units, placements and an open bore. These are meshing fixtures,
with no solver-convergence claim.

## Experimental phase ordering

`thin_region_ordering` is a separate, default-off option. It asks the compatible
Normal Tet4 front to inspect its certified source surface before spending time
on initial interior sites. It does not force that constructor: verified tube
and polygonal-sweep routes remain available. An explicit `boundary_first=True`
takes precedence.

~~~python
import paacs_mesh

model = paacs_mesh.read_step("hollow-part.step")
mesh = model.mesh(
    4.0,
    assembly="preserve_parts",
    thin_region_ordering=True,
    max_chord_error=0.1,
    max_relative_volume_error=0.01,
    max_working_bytes=2 * 1024**3,
    workers=1,
    validation_workers=1,
)
~~~

The CLI exposes `--thin-region-ordering` and `--no-thin-region-ordering`.
Saved recipes retain the independent boolean. It enables no CAD repair,
defeaturing or source-provider fallback.

Leave constructor selection automatic for normal use. The benchmark's explicit
`advancing_interior_front=True` profiles isolate that path and require
`assembly="conforming"`; they are not a recommendation to replace verified
tube or polygonal-sweep construction.

The selector builds one bounded spatial index over the immutable source and
samples up to 64 area-weighted facet locations. Each sample uses the actual
size-field value and the initial constructor's face-based proposal reach.
Opposite source facets intersecting that reach indicate competing wall fronts.
At least 16 distinct facets and 80% blocked samples are required; an ambiguous
ray prevents automatic selection. This is scheduling advice, not proof that
every possible interior point is infeasible.

Selected cases recover the boundary first and then use the existing interior
front on the constrained material mesh. Thick parts, finer local fields and
unresolved samples retain the established order. No source node, face, size
request or validity limit changes. The selector's ceiling is five million
explicit work units and 32 MiB of modeled temporary storage, subordinate to
the enclosing budgets. Cancellation and allocation failures remain visible.

The initial qualification covers global sizing and analytic volume controls.
Selected-face controls can require native CAD projections during size-field
evaluation; the selector abstains with `native_face_sizing_not_qualified` for
these requests. Their existing sizing and construction paths still run.

Automatic boundary-first construction also retains a bounded, topology-only
cleanup for mean ratios below `1e-6`, unless initial topology quality is
explicitly disabled. It can make at most 64 accepted edits and does not move
nodes. Frozen source facets and edges remain protected. This is a best-effort
critical cleanup, not a quality certificate or general optimization pass;
`local_delaunay_quality_edits` records accepted edits and final quality remains
measured independently.

## Size-aware native sweep caps

A separate qualification finding concerns the fallback that sweeps a native
cap triangulation. Geometric tessellation controls chord error, but does not
set a physical triangle size on a plane. Applying the requested spacing only
along the sweep can leave large, thin cap triangles and poorly shaped Tet4s.
The rotated sharp profiles exercise this route without changing the CAD.

Normal single-body requests without local controls first try the existing
size-aware regular surface builder for this cap source. Cap trimming, holes,
original entity classification and the original-CAD sweep validation remain
mandatory. A refusal can make one historical cap attempt within the remaining
request allowance. Cancellation, allocation failure, source mutation and
internal or kernel failures are not retried. `bodies[].native_cap_sizing` records the requested cap size,
whether it was applied, source counts, diagnostic and all attempted cap work
and time. A coarse fallback must not be described as applied cap sizing.

This change is independent of `thin_region_ordering`. The private standalone
sweep API retains its historical source when `maximumCapSize` is zero. The
bounded Python resource policy (`resource_policy="bounded"`) also retains the
established cap dispatch.

## Evidence and diagnostics

Compare construction order at identical STEP bytes, size, chord and volume
limits, assembly policy, worker counts and resource allowances. Record every
attempt and retry, including configuration refusals and cooperative timeouts.
Keep route recognition, source generation, initial placement, recovery,
post-recovery placement, CAD checks and export costs separate.

`bodies[].initial_front.ordering` records selection, sampled/blocked/unknown
counts, work, storage, time and the measured clearance-to-proposal ratio.
`complete` describes the advisory query; it certifies neither mesh validity nor
achieved size. Initial and post-recovery reports also expose their existing
spacing, shape, geometry/cavity and Native refusal counters. A zero insertion
count alone is insufficient to identify the cause.

Measure all-cell mean ratio, dihedral angles, scaled Jacobian and achieved
edges independently. Preserve the original CAD body/face tags, positive
orientation, conforming exterior, open bore and CAD-error checks. Whole-element
cutaways retain real cell IDs. Compare tails as well as medians: a higher median
can coexist with a worse minimum.

## Explicit preparation qualification

The separate benchmark preparation tool defaults to `none`. Explicit
defeaturing operates on selected fillet faces of an owned rounded rectangle;
protected modified or deleted faces cause refusal. Repair is a distinct
operation with a caller-supplied tolerance ceiling. It currently qualifies
guarded copying of valid owned inputs, not repair of invalid raw CAD.

Reports retain original bytes/BRep guards, changed volume, bounds, tolerances,
open-bore probes and provider entity history where available. Repair history
and STEP-reimport face correspondence remain unqualified. Prepared variants
are new models with separate evidence; they do not count as unchanged-source
coverage. The standard importer and mesher do not run these operations.

## Local results: 25 November 2026

The Windows qualification used CPython 3.12, OCCT 8.0.1, MSVC 19.44 in the
C++23 compatibility build, and a Ryzen 9 8945HS with 31.29 GiB of RAM. Meshing
and native validation each used one worker. The main comparison requested
4 mm size, 0.1 mm chord error, 1% relative volume error and 2 GiB working
allowance. The final candidate enables `thin_region_ordering`; its default
remains off. The separate cap-source fix applies to eligible Normal requests
without that toggle.

The baseline completed 20 of 22 models within the explicit 180-second trial
allowance. Both remaining requests were cooperatively cancelled, rather than
diagnosed as unsupported geometry. The final build completed all 22, preserving
all 20 earlier successes and every original body and face tag. Fourteen mesh
archives are byte-identical to the baseline. This is the owned hollow matrix,
not a rerun of the original 51 downloaded models.

| Rounded section | Baseline process time | Final process time |
|---|---:|---:|
| 40 mm length | 53.27 s | 5.45 s |
| 40 mm, rotated | 46.02 s | 5.25 s |
| 160 mm length | Cancelled at allowance | 18.71 s |
| 160 mm, rotated | Cancelled at allowance | 22.33 s |

Total process time was 487.45 s before and 88.23 s afterwards, including the
two baseline cancellations. This is a portfolio cost, not a speedup between
matching completed meshes. Every earlier attempt and refusal remains in the
development ledger.

The four rotated sharp profiles expose a different defect: their historical
cap tessellation applied axial spacing but left oversized transverse edges.
Physical cap sizing raises their median mean ratio from roughly 0.04–0.05 to
0.49–0.50. In the short rotated rectangle, independently measured median cell
edge changes from 30 mm to 3.82 mm, minimum mean ratio from 0.0193 to 0.377,
and minimum dihedral from 0.971 to 13.8 degrees. The mesh grows from 264 to
2,970 Tet4s; the extra runtime buys actual cross-sectional resolution.

All five 2 mm trials also pass. The simple circular, square and rectangular
profiles have minimum mean ratios of approximately 0.61, 0.63 and 0.63, with
minimum dihedrals around 26.5 degrees. Rounded sections retain poor tails:
the long 4 mm case has minimum mean ratio 0.0118 and minimum dihedral
0.0777 degrees. The short rounded case's first percentile also decreases
slightly despite its nearly unchanged median. These are valid meshes under
the recorded checks, not solver qualification or uniformly good elements.

The finer long rounded case exhausts the selector's five-million-work advisory
allowance before sampling and completes through the established order. It
remains a source-index efficiency target; the advisory limit does not refuse
the mesh. Reusing a certified source index and improving curved transition
cells are further work.

Forty-eight focused Python tests and three native regression executables pass.
The latter include whole-request cancellation exactly after a completed cap
attempt, and a bounded source refusal followed by a valid coarse fallback.
Figures independently check every displayed archive's Tet4 orientation,
face incidence and tags, then measure all-cell quality and achieved edges.
Their cutaways retain whole elements and link displayed faces to archive IDs.

The preparation prototype passes two successful controls and five explicit
guard trials. Removing four selected outer fillets changes volume by
549.3807 mm³; that intentionally modified model meshes as a separate input.
Original STEP bytes remain unchanged. Valid-input repair also meshes, but
neither invalid-input repair nor STEP-reimport history is qualified.

The machine-readable measurements, settings, hashes, limitations and retained
cost ledgers are in [the qualification data](data/thin-hollow-qualification.json).
