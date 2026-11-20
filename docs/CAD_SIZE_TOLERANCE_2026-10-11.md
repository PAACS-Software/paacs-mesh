# CAD sizing and geometric tolerance — 11 October 2026

This evaluation distinguishes the geometric approximation of a linear mesh
from the validity of its elements and its discrete topology. The fixed 51-file
rerun reproduces **39/51 accepted meshes**, all 39 archives byte for byte,
without regressions. Every one of the 12 refusals was retried at half and
quarter requested size: **neither size added an accepted original-file mesh**.
All 75 size-only attempts completed, with no timeouts or unexpected volume
publication. Three further explicit-chord trials also produced no accepted
volume. The current regression suite passes 77/77. No production meshing or
acceptance code has changed.

## What the previous refusals mean

The existing original-CAD centroid checks use `model.tolerance(3, bodyTag)`,
usually 1e-7 mm. The native surface path independently permits a sampled
chord error of `min(0.05 * targetSize, 0.0025 * bodyDiagonal)` by default.
Consequently, a candidate may satisfy its recorded linear-surface error bound
and still fail the much stricter original-solid centroid check.

Validation depth also differs by route. The native surface/partition fallbacks
check original-CAD tetrahedron centroids, whereas standalone recognised
primitives use their analytic construction, mapped CAD boundary nodes, sampled
boundary errors, volume agreement and discrete audits. Their zero reported
native-centroid check count must not be read as an exhaustive test with zero
outside centroids. The 39-file accepted baseline records the checks actually
performed; it is not a claim of continuously inscribed cells on every route.

For a concave boundary such as a cylindrical bore, a chord joining two CAD
points cuts across a small part of the exact hole. A linear tetrahedron next
to that chord can have a centroid outside the exact solid while being a
positive, conforming element of the approximating discrete domain. Exterior
convex surfaces can instead produce omitted material. Both directions of
geometric deviation matter; requiring an inscribed mesh alone introduces a
systematic geometric bias.

This corrects the earlier interpretation: an outside centroid is evidence
against strict CAD containment, not sufficient evidence that an element is
invalid or unusable in an engineering simulation. Conversely, simply being
near some CAD face is insufficient to establish correct geometry or body
ownership. Missing holes, merged features, overlapping regions, unmatched
interfaces, inverted/zero cells and incomplete bodies remain separate defects.

The diagnostic records geometric approximation, topology, element size and
element quality separately. Its sampled distances describe the evaluated
points. A continuous error bound requires a certificate over the full
boundary.

## Independent deviation measurements

The analytic bore control passes: a positive Tet4 has all four vertices ON
a valid native annular solid, but its centroid is OUT by **0.002046 mm**.
The corresponding circular chord sagitta is **0.006092 mm**. This local
control demonstrates the difference between strict containment and a bounded
linear approximation. It is not a complete mesh of that solid or a corpus pass.

On CAD 27, the initial and finer PAACS boundary meshes had sampled
mesh-to-CAD deviations of 0.053947 and 0.027777 mm. In the reverse direction
the sampled deviations were 0.047400 and 0.023875 mm. The CAD-to-mesh
witnesses were verified on their original trimmed faces. These surface
measurements are separate from whole-volume acceptance.

The generic body-wide boundary-ON check can refuse points whose measured
shell distance is effectively zero. The finer PAACS boundary has 87 such
body-ON refusals and a maximum refusal-point distance of about 1.8e-15 mm.
The owning-entity node checks provide a separate record of node placement.

The workspace-only `geometry_deviation.cpp` probe measures retained meshes
against the original STEP import. It records the following diagnostics:

- All original-body outside centroids, their distances to actual trimmed CAD
  shells, and whether their tetrahedron touches the discretized boundary.
- Four physical samples on each of up to 512 deterministically spread mesh
  triangles, measured against the CAD shells.
- A 9-by-9 interior parameter grid on each native trimmed face, filtered by
  that face's trim, measured against all mesh boundary triangles. Faces with
  no admitted grid samples are explicitly listed.
- Exact discretized surface embedding, volume boundary incidence, positive
  cells, native classification refusals and stable original volume comparison.
- A diagnostic 0.001–0.1 mm threshold grid. None of those thresholds is selected
  as engineering acceptance, and none adds corpus coverage.

The probe's nine analytic checks verify boundary-shell distance, including an
inside-box point whose distance to the solid itself would misleadingly be
zero, and triangle-interior/edge/vertex distances.

These finite distances are diagnostic estimates. Continuous containment,
unsampled extrema, small-feature preservation and exhaustive volume
intersection checks require separate evidence.
Engineering acceptance also depends on critical feature dimensions, intended
boundary conditions, solver requirements and solution convergence. No universal
engineering tolerance is inferred from these files.

### Independently measured accepted primitives

Two fresh accepted archives were also read by the same probe. Their original
STEP and archive hashes remain unchanged. These are controls within the 39
baseline successes, not additional successes among the 12 refusals.

| Accepted original CAD | Volume Tet4s | Boundary triangles | Outside centroids | Sampled mesh → CAD / mm | Sampled CAD → mesh / mm |
| --- | ---: | ---: | ---: | ---: | ---: |
| 7, cylinder | 3,456 | 1,056 | 0 | 0.194942 | 0.173275 |
| 8126, torus | 438,912 | 36,576 | 0 | 0.004255 | 0.003810 |

Both controls pass the probe's positive-cell, discrete surface-embedding and
volume-incidence checks, with no interior-node IN refusals or unknown native
states. All centroids were tested. This also shows why a zero outside-centroid
count does not certify geometric accuracy: an inscribed approximation can
omit material near a curved exterior. The two models have different dimensions
(diagonals approximately 258.77 and 42.57 mm), so their absolute errors are not
a direct comparison of algorithm quality. The finite sampling and engineering
scope limitations above still apply.

An eventual tolerance-aware acceptance mode should retain positive-cell,
connectivity, body-ownership and interface checks, then evaluate the mesh
boundary against its owning trimmed CAD entities in both distance directions.
Its report should identify the supplied physical tolerance and the validation
depth actually completed. A near-boundary OUT classification can then be
reported as geometric approximation when those other conditions hold; a
missing cavity or wrong body cannot be accepted merely because some CAD face
is nearby. The present probes establish why that distinction is needed, but
do not implement or certify that production acceptance mode.

## Controlled size campaign

### Completed results

| Phase | Original files attempted | Accepted | Geometry/recovery refusal | Resource refusal | Body-interface refusal | Invalid import |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Original size, 0.15 × model diagonal | 51 | 39 | 2 | 6 | 3 | 1 |
| Half size, 0.075 × model diagonal | 12 | 0 | 2 | 6 | 3 | 1 |
| Quarter size, 0.0375 × model diagonal | 12 | 0 | 2 | 6 | 3 | 1 |

Every failure retains its public outcome category at both smaller sizes:

| CAD | Original size | Half size | Quarter size |
| --- | --- | --- | --- |
| 1 | Geometry refusal | Geometry refusal | Geometry refusal |
| 10 | Body interface | Body interface | Body interface |
| 11 | Resource refusal | Resource refusal | Resource refusal |
| 12 | Resource refusal | Resource refusal | Resource refusal |
| 13 | Resource refusal | Resource refusal | Resource refusal |
| 14 | Body interface | Body interface | Body interface |
| 23 | Resource refusal | Resource refusal | Resource refusal |
| 25 | Invalid import | Invalid import | Invalid import |
| 26 | Resource refusal | Resource refusal | Resource refusal |
| 27 | Resource refusal | Resource refusal | Resource refusal |
| 28 | Meshing refusal | Meshing refusal | Meshing refusal |
| 3636 | Body interface | Body interface | Body interface |

The named stage reports provide more information than the public resource
category. CAD 11, 12 and 13 retain nonmanifold native surfaces; CAD 23 retains
a nonconforming edge. CAD 26 and 27 prepare surfaces but their volume fan
refuses nonzero genus, followed by unsuccessful partition attempts. CAD 10,
14 and 3636 retain unresolved body interfaces; CAD 25 fails import validity.
CAD 1 and 28 retain `geometry_error` and `meshing_failed` respectively. The
geometry/recovery failures require investigation at those stages; increasing
resolution did not repair them in this profile. These bounded trials do not
prove that every other size or resource configuration would also fail.

The figure �All 51 original outcomes and both smaller retries� is retained in private development evidence at `out/visual-validation-2026-10-11/size-tolerance/size-campaign.png`.

The frozen 39/51 executable is replayed on the identical 51 original files at
size factor 0.15. Every fresh baseline failure is then attempted at 0.075 and
0.0375, including invalid and overlapping inputs so that all baseline failures
receive the requested smaller-size trials. Size is the only changed meshing
setting. Both smaller phases use the full baseline failure set.
The preferred quality threshold is 0.25 as a warning; the strict minimum
quality option remains its default zero. Positive orientation and the other
validity checks are retained. These trials do not require reaching 0.25 in
order to publish an otherwise accepted mesh.

Each phase retains 300 million reported work, 200,000 nodes, 1 million cells,
128 MiB logical mesh/audit output budgets, 4,096 native pieces, a 5,400-second
direct-child timeout, one meshing thread, eight validation workers and four
independent CAD jobs. These logical budgets do not cap native allocations or
process RSS. Peak direct-child memory is sampled. Source SHA256, executable
SHA256, exact commands, stage reports and publication/refusal outcomes are saved.
Elapsed time and RSS are operational receipts. The baseline overlapped the
regression suite and diagnostic work; isolated performance benchmarking
remains a separate task.

The default surface error cap is important. The CLI target uses the bounding
box of all imported solids, while the native fallback's cap uses the current
body diagonal. At factors 0.15 and 0.075, `0.05 * targetSize` exceeds
`0.0025 * bodyDiagonal`, so that cap remains the same when size is first halved.
For a single-body model, factor 0.0375 changes it to 0.001875 of the diagonal:
quartering the target tightens the default chord allowance by only 25%, rather
than fourfold. For sufficiently smaller bodies in an assembly, even that
quarter-size trial leaves their default cap unchanged. Later fallback
preparation may tighten it further. Actual reported surface settings/counts
accompany the outcomes; smaller requested volume size must not be presented
as guaranteed uniform surface refinement.

The current full CTest suite passes **77/77**, with no failures or skipped
targets, in 845.73 seconds. Test executable hashes are recorded. These results
are separate from the frozen CLI corpus qualification; installed consumer
examples were not rerun during this size/tolerance evaluation.

The completed surface comparisons below show why requested element size and
actual surface refinement must be tracked separately. These are reported
native preparation counts, including refused nonmanifold surfaces; they do
not represent published volume meshes.

| CAD | Original triangles | Half-size triangles | Quarter-size triangles | Retained initial-stage refusal |
| --- | ---: | ---: | ---: | --- |
| 11 | 34,428 | 34,428 | 42,922 | Nonmanifold surface |
| 12 | 34,430 | 34,430 | 42,926 | Nonmanifold surface |
| 13 | 9,825 | 9,825 | 10,956 | Nonmanifold surface |
| 23 | 0 | 0 | 0 | Nonconforming edge |
| 26 | 27,760 | 27,760 | 28,098 | Nonzero-genus volume fan |
| 27 | 3,274 | 3,274 | 3,750 | Nonzero-genus volume fan |

CAD 26 and 27 complete native surface preparation in each trial. Their volume
fan refusal is followed by partition attempts, which retain `invalid_result`
and `volume_mismatch` aggregate diagnostics respectively. No accepted volume
is produced in these comparisons.

### Additional explicit boundary-resolution trials

Three isolated frozen-CLI trials also set size factor 0.0375 and an explicit
surface chord bound. Their direct-child timeout is 1,200 seconds, rather than
the corpus campaign's 5,400 seconds. None timed out. They are separate from
the size-only campaign because an additional parameter changes.

| Original CAD | Explicit chord bound / mm | Reported prepared surface triangles | Outcome |
| --- | ---: | ---: | --- |
| 27 | 0.025 | 11,444 | Resource refusal; no volume archive |
| 27 | 0.005 | 49,214 | Resource refusal; no volume archive |
| 28 | 0.025 | No published volume | Meshing refusal; no volume archive |

For CAD 27, native surface preparation completes, but the single-kernel fan
route refuses the nonzero-genus boundary. Partition attempts report original
volume-conservation mismatches, provenance loss or local partition admission.
Only about 568,824 and 1,896,105 total validation work units are reported for
the two trials, far below the nominal 300-million whole-call ceiling: the
resource category must not be interpreted as that entire ceiling being spent.
These outcomes are construction/partition gaps before an accepted volume,
independent of the stricter centroid interpretation. CAD 28 remains a separate
assembly/recovery refusal. No additional whole-source success is claimed.

## Evidence and reproduction

Evidence root: `evidence/local-2026-10-11/size-tolerance`.
The campaign is restartable under identical pinned inputs and settings:

```powershell
& 'C:\Program Files\Python312\python.exe' evidence/local-2026-10-11/size-tolerance/run_campaign.py
& 'C:\Program Files\Python312\python.exe' evidence/local-2026-10-11/size-tolerance/run_tests.py
```

The read-only diagnostic build and measurements use a fresh build directory;
do not rebuild over an existing receipt:

```powershell
# Build once if build/geometry-deviation-2026-10-11-v2 does not exist.
& 'C:\Program Files\Python312\python.exe' evidence/local-2026-10-11/size-tolerance/build_deviation_probe.py
& 'C:\Program Files\Python312\python.exe' evidence/local-2026-10-11/size-tolerance/run_accepted_deviations.py
& 'C:\Program Files\Python312\python.exe' evidence/local-2026-10-11/size-tolerance/summarize_campaign.py
& 'C:\Program Files\Python312\python.exe' evidence/local-2026-10-11/size-tolerance/render_size_campaign.py
```

`summary.json` verifies phase manifest/record hashes, the unchanged 51 original
STEP files, the recorded production source/binary pins, all 39 accepted archive
hashes, the regression executable hashes, and the retained deviation-measurement
inputs and outputs. Exact native commands and stdout/stderr accompany every
attempt. The complete phase wall times are 77.22, 47.10 and 47.33 minutes;
peak sampled direct-child working sets are approximately 3.88 GiB. These are
operational receipts under concurrent work, with the limitations stated above.
