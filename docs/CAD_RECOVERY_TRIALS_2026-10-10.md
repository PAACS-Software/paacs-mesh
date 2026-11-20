# CAD recovery trials — 10 October 2026

**Tolerance-policy clarification (11 October):** original-CAD centroid checks
in this report use the body's native tolerance, distinct from a geometric
approximation tolerance for a linear mesh. "Rejected" records that strict
policy; it does not by itself establish engineering unsuitability. The fresh
[size and tolerance evaluation](CAD_SIZE_TOLERANCE_2026-10-11.md) measures this
distinction without changing the historical corpus counts.

This records the wire-refinement checkpoint. The later
[facet/kernel iteration](CAD_FACET_KERNEL_TRIALS_2026-10-10.md) adds two optional
construction methods, a solved recovery control and fresh tests.
Original-CAD coverage remains 39/51.

The qualified whole-source baseline remains **39/51**. This iteration adds an
opt-in native-wire subdivision utility and bounded surface and volume trials.
It earns **no additional original-CAD success**. The current CLI and
core library remain byte-identical to the qualified baseline; the optional
adapter library now also contains the new utility. Default meshing dispatch
does not call it. The failed targeted remeshing policy was not enabled there.

## Implemented and tested

`src/geometry/OcctWireRefinement.{hpp,cpp}` subdivides an actual CAD wire by its
native parameters and splits all incident surface triangles consistently.
Original nodes and classifications remain bit-identical. New nodes retain
the owning CAD edge classification and are checked against actual incident
pcurves and trims at source tolerances. Every coowner of an affected face must
participate. Both bodies must have identical parameter chains; canonical
parameter evaluation makes shared coordinates bit-identical despite opposite
local endpoint ordering. All destinations publish together after exact closed,
oriented, embedded surface certificates. Failure preserves every destination.

The default admits native lines. Nonlinear, nonperiodic curve candidates
require explicit opt-in. Periodic curves are refused because their branch and
surface rebuilding need a different method. Fan subdivision is a candidate
construction; owning-CAD facet audits and whole-volume checks remain caller
requirements. This does not establish a general trimmed-spline surface mesher.
New cross-call persistent node keys also remain a consumer responsibility.

The new regression passed **63 checks**, including original source/node
preservation, deterministic output/work, actual shared-face matching and
opposite winding, missing coowner refusal, late second-body rollback, one-less
work, added-node/count/working budgets, cancellation and periodic refusal.
A subdivided CAD cube also produced an independently audited 18-tetrahedron,
unit-volume mesh through the established native-admitted star route.

The rebuilt native boundary, PLC recovery and contact-component regressions
also passed: **4 selected CTest targets passed**, including the new test.
The full suite now contains 76 targets; it was not rerun. The previous 75/75
and 15/15 installed-example qualification applies to the preceding frozen
baseline, not to an unrun new full suite.

## Quick PAACS recovery trials

Original CAD 27 has an embedded, closed prepared boundary with 1,631 nodes
and 3,274 triangles. One direct PLC recovery attempt refuses edge `[0,1]`
with `InvalidCandidate`. Its actual CAD wire is represented as a spline;
the original two boundary endpoints alone do not provide intermediate nodes.

Adding one native-curve midpoint produced a certified 1,632-node,
3,276-triangle surface. The first recovery refusal moved to `[0,4]` with
`ResourceBudgetReached`. Three new points produced 1,634 nodes and 3,280
triangles, but refused the same next edge. Neither produced a volume.
The baseline, midpoint and four-segment trials took about 15.7, 16.0 and
17.3 seconds respectively. Each admitted one core attempt, a 1-million-unit
general recovery reservation, and a 20-million-unit whole-solve work ceiling.
No larger repeated cavity-search campaign followed.

A cube control demonstrates why refinement is opt-in: its direct fixed-PLC
solve passes before subdivision, but the subdivided direct fixed-PLC solve
hits a recovery limit. The established star route handles that same refined
cube correctly. More boundary points can change recovery difficulty and do
not by themselves establish broader compatibility.

## Surface refinement diagnostics

On CAD 27, halving the global deflection increased the prepared surface
from 3,274 to 4,772 triangles. A separate trial ranked original trimmed
faces by outside-centroid witnesses and selected faces 70 and 71. Their
eight incident edges drove a joint patch rebuild at smaller deflection.
The resulting 4,780-triangle surface remained closed and embedded and
passed the same sampled facet-error limit.

CAD 26's separate quick trial hits the deliberately small 10,000-node /
20,000-triangle surface cap before volume construction. That is a trial
admission limit and does not diagnose its production whole-source failure.

## What the experiments change

The core recovery cannot complete CAD 27's prepared PLC under the bounded
trial. A discrete surface certificate establishes the embedded triangular
boundary. Volume publication also needs original-CAD material checks,
including where finer deflection changes the discrete geometry.

The next algorithm should use recovery failure witnesses to introduce
additional edge, face or interior points under a shared bounded controller.
It should evaluate candidate tetrahedra and their nearby boundary patches
against native material during construction, then rebuild only the affected
patches and their coowners. It must retain face trims, actual wire constraints
and the final original-CAD checks. Curved shared interfaces require one common
patch for both bodies.

No arbitrary-CAD guarantee follows from these trials. Continuous containment
against the exact CAD boundary and exhaustive global tetrahedron intersection
checks were not certified. Quality thresholds were not relaxed into validity
claims, and source geometry was not repaired or replaced.

## Local reproduction and evidence

The reproducible build driver, input-hash-checked runner and witness ranker
live under `evidence/local-2026-10-03/coverage-40/`:

- `build_wire_refinement_trial.py` builds a fresh isolated diagnostic pair
  from `wire-refinement-trial-sources/`, with copied native/core libraries.
- `run_wire_subdivision_trial.py` runs bounded controls or fixed original
  CADs. `--build-directory` and `--build-receipt` support a fresh reproduction.
- `rank_native_failure_faces.py` proposes bounded actual face selections
  from recorded outside-centroid distances.
- `wire-refinement-qualified-state.json` binds current source/test hashes
  and unchanged qualified CLI/core hashes.
- `wire-subdivision-trials/*/receipt.json` records inputs, controls,
  stage timings, sampled peak process memory, validation and refusal details.
- `wire-refinement-final-tests.xml` and the neighboring test receipts retain
  the fresh checks.

Native runs use the runtime paths in
`evidence/local-2026-10-02/occt-sdk-setup.json`. The build uses MSVC's
`/std:c++latest` C++23 compatibility mode; this is not C++26 qualification.
The 256 MiB working admission is logical payload admission, not a process
memory cap. Explicit work counts do not cap OCCT or predicate internals;
the runner enforces a 180-second timeout for each child stage.

```powershell
python evidence/local-2026-10-03/coverage-40/build_wire_refinement_trial.py --directory coverage-40-wire-subdivision-trial-reproduction
python evidence/local-2026-10-03/coverage-40/run_wire_subdivision_trial.py --case 27 --mode core --subdivisions 2 --label reproduced-midpoint --build-directory coverage-40-wire-subdivision-trial-reproduction --build-receipt build/coverage-40-wire-subdivision-trial-reproduction/build-receipt.json
cmake --build build/windows-msvc-release --config Release --target paacs-occt-wire-refinement paacs-occt-boundary paacs-occt-plc-recovery paacs-occt-contact-components --parallel 2
ctest --test-dir build/windows-msvc-release -C Release -R "^paacs.geometry.(native_wire_refinement|native_boundary|plc_recovery|contact_components)$" --output-on-failure
```
