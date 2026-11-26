# Mesh inspection and comparison

**RC4 candidate APIs.** These additions require RC4 bindings; published
`1.2.0rc3` wheels retain their original interface. Development wheels require
matching source receipts. Inspection itself does not change meshing algorithms
or recertify a supplied report.

Inspection provides independent Tet4 shape measurements, stable element
queries, optional CAD proximity, reproducible SVG views and comparisons. It
does not retry meshing, refine elements, change geometry or replace the native
validity checks. NumPy is not required.

## Start with a native mesh

```python
import paacs_mesh as pm

cad = pm.read_step("tests/fixtures/cad/cylinder.step")
mesh = cad.mesh(.5, assembly="conforming", workers=1, validation_workers=1)

with pm.inspect_mesh(mesh) as inspected:
    summary = inspected.summary()
    poor = inspected.query(
        metric="mean_ratio", quality_below=.2,
        sort_by="mean_ratio", limit=50,
    )
    print(poor["total"], "matching cells")
    for cell in poor["cells"]:
        print(cell["cell_id"], cell["node_ids"], cell["metrics"])

    inspected.save("cylinder-inspection")  # A new directory is required.
    inspected.write_vtk("cylinder-quality.vtk")
    inspected.visualize(
        "cylinder-quality.svg", metric="mean_ratio", highlight_below=.2,
        cut_axis="y", cut_value=0,
    )
```

`inspect_mesh(source, **limits)` accepts a native `Mesh`, an archive filename or
an inspection bundle directory. `open_mesh(path, **limits)` opens the same file
formats. Both return a `MeshInspection`; use its context manager or `close()` to
release the private snapshot. The input is copied before inspection so changes
to the caller's file cannot silently change measurements or element IDs.

An inspection bundle contains `mesh.paacsmesh`, `report.json` and
`inspection.json`. Mesh and report hashes, identities and counts are checked on
opening. `save()` requires a paired native report and never replaces an
existing destination. A bare archive cannot manufacture that report.

## Identity and validation provenance

Every query includes the mesh SHA-256 and a declared ID scheme:

- `cell_id` is the zero-based original archive cell index. Boundary triangles,
  lines and other cell families occupy indices too; it is not a Tet4-only row
  number.
- `node_ids` are the zero-based original archive node indices. Filtering,
  sorting, pagination and diagnostic export do not renumber them.
- An ID belongs to its mesh hash. Cell 100 in a second mesh is not assumed to
  correspond to cell 100 in the first.
- CAD body/face tags and node classifications are retained as reported. A
  positive dimension-3 tag provides a source-local body label; inspection does
  not independently rebind it to a different CAD revision.

A live native mesh retains its source hash, units, kernel version and native
report. A saved bundle retains a hash-paired copy of that information. Hashes
detect mismatched files; they are not signatures authenticating the producer.

Archive v1 itself stores no units, CAD source hash, kernel version or native
validation report. These remain unavailable for a bare archive. Its geometry
can still be inspected and rendered, with units explicitly unknown, but guarded
CAD proximity and source-bound comparisons need the missing provenance.

`inspection.validation` distinguishes a live native report from a paired saved
report and from unavailable validation. It records
`independent_native_reaudit_performed=False`: reopening checks archive
structure and provenance, not the original CAD embedding, body overlap or
native audit again. Shape statistics are separate from these native checks and
from solver accuracy or convergence.

## Metrics

Measurements currently apply to **linear Tet4 only**. Other cell families are
counted and remain queryable, but unsupported shape measurements are `None` in
JSON, not zero. A genuinely inverted or degenerate Tet4 has its own orientation
classification. Its defined zero mean ratio must not be confused with an
unsupported family.

Let `D = dot(cross(p1-p0, p2-p0), p3-p0)`, `V = D/6`, and `Lij` be the six edge
lengths. The implementation uses scale-aware arithmetic and an exact
binary-coordinate fallback for difficult orientation cases.

| Field | Convention |
| --- | --- |
| `mean_ratio` | `12 * (3*V)^(2/3) / sum(Lij^2)` for positive orientation; zero otherwise. A regular Tet4 has value 1. |
| `minimum_dihedral_deg` | Minimum of the six geometric interior dihedral angles, in degrees. A regular Tet4 has approximately 70.529°. Degenerate cases are unavailable. |
| `scaled_jacobian` | `sqrt(2) * D / max(L01*L02*L03, L01*L12*L13, L02*L12*L23, L03*L13*L23)`. Signed, dimensionless, regular Tet4 = 1; negative indicates inverted orientation. |
| `signed_volume` | `D/6`, in input length units cubed. |
| `edge_min`, `edge_max` | Smallest and largest of the six physical edge lengths. |
| `edge_median` | Arithmetic mean of the third and fourth sorted edge lengths. |
| `edge_rms` | `sqrt(sum(Lij^2)/6)`. |
| `volume_equivalent_edge` | `(6*sqrt(2)*abs(V))^(1/3)`, the edge of a regular Tet4 with the same absolute volume. |

The scaled Jacobian uses one common denominator, the largest incident
edge-length product. It is not an unsigned value or a minimum of separately
normalized corners. These linear Tet4 metrics do not evaluate curved Tet10 or
Hex27 Jacobians. Inverted cells can still have geometric edge and dihedral
measurements; those do not make their orientation valid.

`range_limited` identifies a nonzero measurement that cannot be represented as
a finite nonzero binary64 value. Some bounded quality values can round to zero
with that flag set; use the orientation and range flag to distinguish this
from degeneracy. Unavailable measurements are excluded from finite
distributions. Their `count` describes the available population; spatial bins
and comparisons also include explicit unavailability counts. Summary and
query results include the versioned `metric_conventions` object. See
[the implementation](../python/paacs_mesh/quality_metrics.py) for arithmetic
details. Small numerical differences from the retained native report can occur
for almost-degenerate cells.

Distributions provide count, minimum, p01, p05, median, p95, maximum and mean.
Quantiles use linear interpolation at rank `p*(n-1)` in the sorted finite
population. Histogram edges are left-closed: a value equal to an edge enters
the higher bin.

## Query elements and regions

```python
with pm.open_mesh("cylinder-inspection") as inspected:
    page = inspected.query(
        quality_below=.2, metric="mean_ratio",
        region=[[-1, -1, 0], [1, 1, 1]],
        spatial_mode="centroid", body_tags=[1],
        sort_by="mean_ratio", offset=0, limit=100,
    )
    if page["cells"]:
        cell_id = page["cells"][0]["cell_id"]
        detail = inspected.query(cell_id=cell_id)["cells"][0]
        print(detail["points"], detail["classification"])
```

`quality_below` is a strict `<` filter. The three selectable quality metrics are
`mean_ratio`, `minimum_dihedral_deg` and `scaled_jacobian`. Lower values are the
poor end of each convention; thresholds have the chosen metric's units.

Queries default to `cell_type="tet4"`, ascending `cell_id` and a page limit of
100. `limit` must be 1–1000; `offset` is nonnegative. `total` is the exact number
of matches, and `next_offset` identifies the next page. Sorting by a quality
metric resolves ties with the original cell ID and keeps unavailable values
last. `descending=True` reverses the requested ordering.

Use `cell_type=None` to include all supported archive cell families and
`orientation="positive"`, `"inverted"`, `"degenerate"` or `"unsupported"` to
filter orientation. A direct `cell_id` lookup bypasses the default Tet4 type
filter; other supplied filters still apply.

Region mode `centroid` tests inclusion of the arithmetic node centroid. For
Tet4 this is its geometric centroid. `aabb_overlap` tests overlap of the cell's
axis-aligned bounding box. Neither mode establishes exact intersection of the
entire cell with the requested region or CAD.

`cad_boundary_associations` are exact matches of node IDs to stored boundary
triangles. They are not guessed from coordinate proximity. Set
`associations=False` when that optional index is unnecessary.

## Optional trimmed-CAD proximity

```python
cad = pm.read_step("tests/fixtures/cad/cylinder.step")
wall = cad.select_faces(surface="cylinder", radius=1, expect=1)

with pm.open_mesh("cylinder-inspection") as inspected:
    near_wall = inspected.query(
        quality_below=.2, cad=cad, cad_face_tags=wall.tags,
        cad_distance=.1, max_cad_work=100_000, max_cad_seconds=10,
        limit=100,
    )
    for cell in near_wall["cells"]:
        print(cell["cell_id"], cell["nearby_cad_surfaces"])
```

Proximity requires the same known source SHA-256, kernel version and millimetre
units. Bare archives and older reports lacking kernel provenance are refused.
OCCT measures point-to-shape distances on the truly trimmed selected faces,
including their edges and vertices. It does not measure distance to an infinite
supporting plane or untrimmed surface.

The sample is the arithmetic mean of the cell's nodes. It is the Tet4 geometric
centroid, but a node mean for nonlinear and other families. Results include the
sample point, nearest selected face tags, closest point, distance and native
tolerance. Ties use the declared native-tolerance/numerical rule. This is not a
minimum over the whole cell, boundary ownership, a CAD-error bound, proof of a
quality defect's cause or a validity certificate.

`cad_distance` retains distances `<=` the supplied value. Distances are evaluated
for **all candidates surviving the other filters before pagination** so totals
remain exact. Reducing `limit` alone does not reduce that work. Restrict a region,
quality filter, cell ID or face list when needed. `cad_include_faces=True` adds
every selected-face measurement and increases output/storage requirements.

For an explicit cell list, the lower-level public function is:

```python
with pm.open_mesh("cylinder-inspection") as inspected:
    cells = inspected.query(quality_below=.2, limit=10)["cells"]
    measurements = pm.query_facial_proximity(
        inspected, [cell["cell_id"] for cell in cells], cad,
        face_tags=wall.tags, max_work=100_000, max_seconds=10,
        include_faces=True,
    )
```

Call it while the inspection context is open. The admitted work is bounded by
`candidate_count * selected_face_count`, with at most 100,000 CAD evaluations
and 1,000 selected faces. All selected measurements must finish or the operation
refuses; it does not return a partial nearest-face answer.

## Requested, effective and achieved sizing

`summary()` describes the entire snapshot. Its `sizing` block separates:

- `requested`: global/face/box requests, the optional minimum and growth
  overrides. `None` preserves an unset caller override.
- `effective`: field targets after the reported minimum clamp and growth
  inheritance, plus `resolved_routing` from native body reports.
- `achieved`: per-Tet4 edge/volume-proxy distributions and boundary-face edge
  distributions measured from the archived coordinates.

Boundary distributions use **triangle-edge incidences**, so a shared edge
appears once for each incident triangle. They differ from a unique-edge
population. Tet4 distributions use one measurement per cell. A small
volume-equivalent edge can result from a thin sliver and does not certify a
small maximum edge or good accuracy.

An effective target or an observed construction route does not certify every
achieved edge length. `sizing_completion_certificate` remains unavailable, and
these operations do not perform a neighbouring-cell growth audit.

Native `effective_policy` preserves wrapper sizing choices separately from
observed per-body routes. The existing `sizing.mode` is the wrapper policy
before automatic native dispatch; it can say `standard` while the actual body
route is `native_spatial_sizing`. Per-body `initial_front` and `interior_front`
retain actual termination diagnostics and counters. `not_attempted` with zero
remaining counts does not certify completion. See
[the refinement findings](REFINEMENT_FINDINGS.md) for the reproduced RC3 limits.

## Genuine SVG views

```python
with pm.open_mesh("cylinder-inspection") as inspected:
    receipt = inspected.visualize(
        "cutaway.svg", metric="mean_ratio", highlight_below=.2,
        cut_axis="y", cut_value=0,
        view_bounds=[[-1.1, -1.1, -.1], [1.1, 1.1, 2.1]],
        azimuth=35, elevation=25,
    )
    print(receipt["path"], receipt["highlighted_cell_ids"])
```

`pm.visualize_mesh(inspected, path, **options)` is the equivalent function.
Rendering uses the original coordinates, node IDs, Tet4 connectivity and
measured values. There is no illustrative replacement mesh, random jitter or
artificial interior reconstruction.

The context is translucent native boundary triangles, or boundary faces derived
from actual Tet4 cells when triangles are absent. `cad=cad` supplies optional
source-matched CAD display tessellation instead. That display approximation is
labelled with its deflection; it is not a rigorous mesh-error bound. Unsupported
complete boundary contexts require matched CAD rather than silently omitting
part of the object.

A cut keeps whole Tet4 cells with `centroid[axis] <= cut_value` and exposes their
original faces. It produces a stepped surface, not a fictitious planar slice.
All threshold-selected cells are drawn as an explicitly labelled **x-ray
overlay**, including interiors and cells removed from the context by the cut.

`highlight_cell_ids=[...]` isolates explicit native IDs and overrides threshold
selection. It requires supported Tet4 measurements. A close-up uses world
`view_bounds`; off-screen drawing is clipped to the viewport, with cropping
labelled. All selected IDs and values remain in SVG metadata. Reuse the same
bounds, angles, width and height for before/after camera matching.

Numeric legends show metric values and conventions, including negative/zero
scaled-Jacobian separation and dihedral degrees. Hover titles and
`data-cell-id` attributes identify actual cells. A single-cell caption also
prints its ID and measured value for PNG or phone viewing. SVG metadata includes
the exact mesh/source hashes, units, selected cell coordinates/values and camera
frame. The receipt records the SVG hash and byte count.

Defaults are 1200×900 pixels, 100,000 context/cut triangles, 10,000 highlights
and a 32 MiB SVG output limit. Exceeding a limit refuses the render; poor cells
are not silently sampled away. The API writes SVG, without a runtime plotting
or NumPy dependency. PNG conversion is an external display step.

## Compare two meshes

```python
with pm.open_mesh("before-inspection") as first, pm.open_mesh("after-inspection") as second:
    comparison = pm.compare_meshes(first, second, spatial_bins=4)
    frame = comparison["spatial_comparison"]["bounds"]
    first.visualize("before.svg", view_bounds=frame)
    second.visualize("after.svg", view_bounds=frame)
    print(comparison["quality_delta"], comparison["counts_delta"])
```

Comparison requires the same known CAD source hash, matching known units and
paired reports whose native status is reported valid or valid with warnings.
That requirement is provenance, not a fresh native re-audit. Source-local body
tags are retained without independent rebinding. Body-level CAD-error deltas
require matching known kernel versions too. Unknown or differing kernels leave
those measurements as separate unpaired lists, with body deltas unavailable;
overall quality and spatial population comparisons can still proceed.

Results include cell counts, quality/achieved-size distributions, reported CAD
errors, assembly/interference information and changes in common spatial bins.
Deltas are **second minus first**; `None` means unavailable or out of range.
Both snapshots share the union bounding box for centroid bins. The spatial
comparison's low mean-ratio threshold is 0.2. No element-to-element
correspondence is inferred between meshes.

Compare achieved sizes and reported CAD errors alongside quality and counts.
Finer elements alone do not establish better element shape or solver accuracy.
CAD errors in this report are retained native measurements, not independently
recomputed surface distances or a global Hausdorff certificate.

Native timings describe the reported final native meshing call. The comparison
does not accept a source/final-report-paired workflow history, so
`workflow_total_seconds`, `failed_retry_seconds` and `retry_count` remain
unavailable. `end_to_end_runtime_comparison_available` is false. Do not replace
these missing values with zero or claim a runtime speedup including retries
from this report alone.

## Command line

The RC4 interface provides these commands through `paacs-mesh` or
`python -m paacs_mesh`:

```sh
paacs-mesh inspect-mesh cylinder-inspection --summary --json
paacs-mesh inspect-mesh cylinder-inspection --quality-below 0.2 --sort-by mean_ratio --limit 50 --json
paacs-mesh inspect-mesh cylinder-inspection --cell 3295 --json
paacs-mesh inspect-mesh cylinder-inspection --region -1 -1 0 1 1 1 --spatial-mode centroid --json
paacs-mesh inspect-mesh cylinder-inspection --source tests/fixtures/cad/cylinder.step --cad-face 1 --cad-distance 0.1 --json
paacs-mesh inspect-mesh cylinder-inspection --vtk cylinder-quality.vtk --json
paacs-mesh visualize cylinder-inspection --output cutaway.svg --color-by mean_ratio --highlight-below 0.2 --cut-axis y --cut-value 0 --json
paacs-mesh visualize cylinder-inspection --output closeup.svg --cell 3295 --view-bounds 0.25 -0.99 -0.02 0.55 -0.82 0.11 --json
paacs-mesh compare before-inspection after-inspection --output comparison.json --json
```

Replace example IDs with IDs returned from the inspected mesh. In particular,
3295 is meaningful only in the corresponding refinement example, not every
cylinder mesh. `--body` and `--cad-face` may be repeated. Region/view bounds are
`x0 y0 z0 x1 y1 z1`. `--cell-type all` includes unsupported quality families;
their unavailable measurements stay explicit.

`--summary` always describes the entire snapshot and rejects query filters.
`--vtk` exports the **whole snapshot**, even when a filtered page is requested.
The VTK export preserves original cell/node IDs, coordinates, classifications
and measured fields; unavailable measurements are NaN fields rather than
fabricated zero quality. It is a diagnostic export, not a solver validation.

`--json` emits one JSON document to stdout. Native logging is routed to stderr.
Exit codes are 0 for a completed operation, 1 for a refused operation and 2 for
invalid input/usage. A successful inspection command is not a declaration that
an arbitrary input archive is a native-valid simulation mesh.

## Limits and output behaviour

Snapshot defaults are `max_bytes=256 MiB`, `max_cells=1_000_000`,
`max_nodes=1_000_000`, `max_work=20_000_000` and `max_seconds=60`. Query and
comparison operations admit their own scratch before allocating it and check
cooperative wall deadlines. Spatial bins are limited to 1–16 per axis. JSON
reports have an 8 MiB output ceiling; diagnostic VTK defaults to 128 MiB.

These are bounded PAACS storage/work allowances, not process RSS limits. A
comparison retains two input snapshots and has an additional scratch allowance
bounded by the smaller input `max_bytes`; it is not a single combined memory
ceiling for both inputs and Python. Exact-total filtering, sorting, requested
associations and CAD proximity can require substantial additional scratch.

An OCCT distance solve or CAD tessellation and its internal allocations cannot
be interrupted or strictly RSS-bounded by these options. Deadlines are checked
cooperatively, including native operation-gate waiting, and can overrun while
one kernel call finishes. Restrict candidate cells/faces or increase an explicit
allowance when a diagnostic is refused; no hidden retry is performed.

`InspectionError` carries a `reason` and structured report for bounded operation
refusals. Invalid arguments raise ordinary input errors. SVG/VTK and CLI
comparison-file outputs use temporary sibling files and atomic replacement;
failures preserve the previous destination. All known input paths, including
archive, snapshot, bundle provenance and source STEP paths, are protected from
replacement. Existing unprotected output files may be replaced atomically;
inspection bundle directories must be new.
