# Local facet and volume recovery — 10 October 2026

**Tolerance-policy clarification (11 October):** the original-CAD material
refusals below use native CAD tolerances, not an engineering approximation
allowance for straight-sided curved-CAD elements. The fresh
[size and tolerance evaluation](CAD_SIZE_TOLERANCE_2026-10-11.md) reports
physical deviations separately. Historical successes/refusals are unchanged.

Two reusable development methods were added: CAD-evaluated facet insertion and
an optional native-admitted interior kernel point for fixed-PLC recovery. The
second resolves a reproduced small solver failure. **Qualified original-CAD
coverage remains 39/51.** CAD 27 still fails after local facet refinement.
Neither method is enabled in the frozen
qualified CLI. The 40-file milestone remains outstanding.

## Implemented methods

`src/geometry/OcctFacetRefinement.{hpp,cpp}` takes explicitly selected source
triangles. It evaluates their centroids in a canonical coordinate order,
projects onto the actual trimmed native face, checks native evaluation and a
caller-supplied displacement limit, and requires the projected point's planar
foot strictly inside the original triangle. It replaces each admitted triangle
with three directed triangles. Original nodes/classes and every unaffected
directed triangle retain their values and relative order.

Shared faces require every actual body coowner and bitwise matching original
triangles. The inserted coordinate is evaluated once and reused across owners;
local node indices may differ. Exact closed, oriented, embedded surface
certificates run before and after. Budget, cancellation, projection, unmatched
coowner or geometry failures discard all private copies; publication is a batch
of nonthrowing swaps. The caller retains CAD facet-error, material, whole-volume
and persistent cross-call identity gates. This helper does not create a general
trimmed-face triangulation or repair source cracks, seams or missing pcurves.

`OcctPartitionRecoveryOptions::nativeKernelSeed` is a private, default-off
construction option in `OcctRecoveryMesher`. It attempts the existing exact-PLC
star/kernel construction with native material admission. An admitted kernel
must additionally classify strictly inside the native piece before becoming
one new volume point. It changes no boundary coordinates or classifications.
The ordinary Delaunay construction, fixed-constraint recovery, eligible
diagonal retries, exact boundary audits, positive volume, and both original and
piece material gates still run. Native seed queries and explicit work are
included in the recovery receipt.

A geometric kernel-search refusal leaves the unseeded route available. This is
not a seed algorithm for arbitrary concave or multiply connected volumes.
Cancellation, resource, allocation, invalid-boundary, native-query and internal
failures propagate. The point must be retained by the final mesh. Work ceilings
are explicit-visit admission, not a total CPU or native-heap bound.

## A solver failure reproduced and fixed

The control is a unit cube with two coplanar native-face triangles subdivided
at CAD-evaluated centroids: 10 surface nodes and 16 triangles. That surface
passes closure, embedding, source-coordinate, CAD facet and unit-volume gates.
The standard star route can already mesh it; the experiment targets the
separate fixed-PLC recovery route.

| Identical input surface and budgets | Result |
| --- | --- |
| Direct recovery, kernel option disabled, up to four attempts | First constrained edge `[1,3]` reaches the general recovery limit; no volume or source mutation |
| One native-interior kernel point, up to four attempts | Two recovery attempts and an eligible surface-diagonal retry complete a mesh: **11 nodes, 16 Tet4**, unit volume |

The accepted trial passes exact discretized surface and volume incidence,
positive tetrahedra, owning native boundary-node tolerances, strict interior
point inclusion, every tetrahedron centroid, CAD facet samples, stable native
volume and unchanged source BREP. Relative volume error is approximately
`1.11e-16`. This is a regression fixture, not a new original-dataset success.
Preferred-quality optimization was not requested. Independent visualization
measurement reports volume CV 43.0%; uniformity was not the completion target.

The algorithm was first isolated in a diagnostic source copy, then moved to
the private optional adapter. The final `v6` trial and regression executable
use the actual adapter method; their only recovery source override prints
failure witnesses. With the option disabled, the identical four-attempt
control still refuses atomically.

## CAD 27: measured negative result

The selection came from the source surface triangles incident to recorded
native-OUT tetrahedron centroids. This produced 100 unique triangles.
Source IDs, input hashes and selection are recorded in
`facet-refinement-witness-selection.json`.

The all-or-nothing call refused at native face 46, triangle 2101, with
`native-trim-projection-failed`. An isolated controller allowed at most eight
atomic attempts within one 20-million explicit-work scope, recording every
omission. Its second attempt omitted that one unprojectable candidate and
inserted **99 points / 198 additional triangles**. All original coordinates and
classes remained bit-identical. Maximum admitted displacement was 0.036164 mm;
the resulting 1,730-node / 3,472-triangle surface passed exact embedding and
sampled native facet checks.

The PAACS recovery trial with local facet candidates and the optional
kernel method still fails its first constrained edge `[0,1]`, reason 26
(`InvalidCandidate`), after 17.391 seconds. It consumes 6,219,339 recorded
recovery work units, below the shared 20-million scope. CAD 27 has no admitted
star kernel in this trial. Increasing the retry count did not bypass that
ineligible edge failure. The gap remains recovery of the actual CAD constraint,
not an aspect-ratio threshold.

## Fresh validation and baseline scope

`paacs.geometry.native_facet_refinement` now has **158 checks**. They cover
planar/curved evaluation, unchanged source and node prefixes, unaffected
triangles, deterministic work and mesh, duplicate selection, coowner matching
and opposite winding, missing/mismatched coowners, late rollback, work/count/
payload admission, cancellation, displacement refusal, the reproduced recovery
failure, the admitted point/retry fix, and its exhausted-budget rollback.

All **five selected targets passed** on the current adapter, in 29.36 seconds:
native boundary, wire refinement, facet refinement, PLC recovery and contact
components. CTest now registers 77 targets; **the full 77-target suite was not
rerun**. The earlier full 75/75 tests and 15/15 installed examples apply to the
frozen 39/51 checkpoint. No new installed-example or complete-corpus claim is
made for the optional development strategies.

The production CLI and core library remain byte-identical to that checkpoint:

- CLI SHA256: `84fb5fc9ce71e7b4090d139e16b6177134e728f4cb4826bd6723959594439975`
- Core SHA256: `e864b4ad7052a44a01196724acaf3bb7023a029b29febc6db3c239d14f2306eb`

The optional adapter library changed and is separately pinned. The old source
ZIP predates these additions. The MSVC build still uses `/std:c++latest` in
C++23 compatibility mode; this is not C++26 qualification. The logical
256 MiB trial admission is not an RSS cap. Child stages are limited to 180
seconds and their sampled peak process memory is recorded.

## Next direction

The positive control shows that interior construction points and recovery
retries can fix a real constraint-recovery failure without relaxing validity.
It does not establish the correct insertion policy for CAD 27. The negative
facet result argues against adding points everywhere a previous fill fails.
The next broader method needs a coordinated failed-edge/facet cavity strategy
with new recovery points and native material admission during construction,
then common shared-face publication. It must preserve hole/trim topology and
strict final source geometry checks.

## Reproduction

Evidence lives under `evidence/local-2026-10-03/coverage-40/`. Relevant receipts
include `facet-kernel-strategy-{build.log,tests.xml,tests.log}`,
`facet-kernel-qualified-state.json`, `facet-refinement-witness-selection.json`
and `facet-refinement-trials/*/receipt.json`. Trial builders create fresh
directories and copy the current libraries before compiling. `v1`–`v5`
diagnostic snapshots remain separate from the final adapter-backed `v6`.

```powershell
cmake --build build/windows-msvc-release --config Release --target paacs-occt-facet-refinement paacs-occt-plc-recovery paacs-occt-wire-refinement paacs-occt-boundary paacs-occt-contact-components --parallel 2
ctest --test-dir build/windows-msvc-release -C Release -R "^paacs.geometry.(native_facet_refinement|native_wire_refinement|native_boundary|plc_recovery|contact_components)$" --output-on-failure
python evidence/local-2026-10-03/coverage-40/build_facet_refinement_trial.py --directory coverage-40-facet-kernel-reproduction
python evidence/local-2026-10-03/coverage-40/run_facet_refinement_trial.py --case box --mode core --refine-triangles 0,1 --star-interior-seed --core-attempts 4 --label reproduction --build-directory coverage-40-facet-kernel-reproduction
```

Native runtime paths come from `evidence/local-2026-10-02/occt-sdk-setup.json`.
The runner installs those paths only in its child environment and checks
executable, source, library and original STEP input hashes.
