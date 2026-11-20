# Historical readiness notes

These notes describe earlier development revisions. Current release scope is in
[RELEASE_POLICY.md](RELEASE_POLICY.md) and [RELEASE_READINESS.md](RELEASE_READINESS.md).

> Current candidate: **1.0.0-rc.8**. New scope: [local region/face sizing](LOCAL_SIZE_CONTROLS.md) and [planar prism/hex layers](PRISM_LAYER_CONTROLS.md). See [current validation](../VALIDATION.md) and [next priorities](STRUCTURAL_CFD_ROADMAP.md). Earlier release notes below are historical, not new qualification claims.

# Release readiness and development priorities

The historical rc.6 baseline was a modernized engineering candidate. A compiler
upgrade is useful for diagnostics, optimizations, and supported language tools;
it does not prove a new meshing algorithm or universal robustness.

The current prioritized [structural/CFD roadmap](STRUCTURAL_CFD_ROADMAP.md)
sets concrete acceptance gates for CAD, sizing, layers, higher-order elements,
solver accuracy and acceleration. rc.6 delivers the bounded
[uniform-cylinder route](UNIFORM_CYLINDERS.md), not those future milestones.

## Changes with direct correctness or efficiency value

1. Archive counts formerly reached reserve() before the payload was checked.
   The parser now checks lengths, counts, flags, classifications, cell types,
   and node references before allocating mesh storage. A span avoids copying
   caller buffers, fixed storage removes per-cell temporary allocation, and
   explicit budgets are available for untrusted files.
2. Mesh::addNode and addCell could leave parallel arrays/connectivity out of
   sync after allocation failure. They roll back logical additions now.
   Dedup insertion/rebuild also preserve registered nodes after failure.
3. Chunk-count ceiling division could wrap, chunk endpoints could wrap, and
   parallelReduce divided by zero for chunkSize=0. All now use checked-by-
   construction arithmetic. The shared thread count is atomic and worker
   context is restored by RAII.
4. PointGrid could convert NaN/infinite radius calculations to integers and
   enumerate enormous empty ranges. Finite validation and a populated-bucket
   fallback now bound that work; extreme distances use hypot.
5. C API boundary tag addition could overflow signed int32. Primitive requests
   also accepted infinite bounds/non-finite centers. They are rejected early.
6. CLI numeric conversion formerly accepted valid prefixes with trailing junk
   and unsigned negative values. Full-string from_chars validation replaces it.
7. Test coverage previously exposed only four CTest probes. The complete
   inherited unit suite and new failure-path tests are registered now.

## Required before calling this a stable general-purpose CAD mesher

| Priority | Work | Acceptance evidence |
|---|---|---|
| P0 | Keep the full suite passing as development continues | Reviewed changes, continuous CI, unchanged geometric contracts |
| P0 | Decide supported API/ABI surface and license | Versioned contract, owner-approved license, deprecation policy |
| P0 | Qualify supported platforms and sanitizers | Linux GCC/Clang; intended Windows/macOS builds; working TSan run |
| P0 | Expand hostile-input and allocation-failure coverage | Continuous parser fuzzing; failure injection at remaining pipeline boundaries |
| P1 | Extend the qualified planar STEP/BREP importer to curved CAD, seams, assemblies and healing | Independent industrial CAD corpus, curved-boundary fidelity, watertight boundaries |
| P1 | Robust constrained recovery/refinement on difficult domains | Cavities, slivers, near-coincident and thin features, multiple materials |
| P1 | Resource and cancellation guarantees across all meshing routes | Bounded memory/work and transactional cancellation on every route |
| P1 | Measure meshing performance | Equal geometry, target quality, boundary fidelity, CPU/thread budget |
| P2 | Remove remaining full-pipeline topology conversion costs | Profiled evidence, compact representation handoff and memory accounting |
| P2 | Improve parallel Delaunay conflict scheduling | Audited speedup vs SerialExact with stable deterministic outputs |
| P2 | General anisotropic and curved-layer support | Metric-quality and collision/junction corpus qualification |

The worker pool still has the documented one-external-dispatcher contract;
this change does not introduce concurrent independent meshing sessions.
Archive parsing validates structure, not geometric element quality. Run a
mesh audit before treating loaded elements as simulation-ready. File writes
report close failures but are not atomic replace transactions; use a temporary
path and application-level rename if preserving an existing file is required.

The new geometry route recognizes complete affine blocks and straight sweeps;
it does not automatically decompose arbitrary solids into conformal regions.
The controlled planar CAD corpus is a first qualification layer, not proof of
general industrial CAD coverage. The pyramid corpus exposes a very low-quality
tetrahedron; optional quality thresholds now refuse such output, but guaranteed
quality improvement remains work. See `GEOMETRY_GUIDE.md` and `../VALIDATION.md`.

SPD metrics currently feed a conservative scalar adapter; general anisotropic
Delaunay is not implemented. Prepared layers do not automatically solve CAD
corner/junction classification. Sliver removal is bounded. The pre-OCCT limits
in `docs/v0.10/KNOWN_LIMITATIONS.md` still apply unless explicitly superseded.

The rc.3 boundary audit keeps exact incidence validation while reducing
sorting and temporary storage. Per-body sizing and graded straight layers
improve control of existing routes; they do not add arbitrary curved CAD or
anisotropic Tet meshing. Coordinate-resolution regressions also exposed why
ideal affine/translation formulas need conservative emitted-cell fallback.
See `AUDIT_DESIGN.md` for the precise boundary-audit contract and cost model.

The rc.4 primitive route adds explicit cylinders, cones, frusta, tubes, tori
and elbows with accuracy-aware resolution, periodic connectivity and full
validation. It is not automatic curved-CAD import or general solid decomposition.
Cone apex fans need local coarsening/remeshing to maintain quality at high
resolution. The new Bernstein Hex8 test is a conservative floating-point
sufficient test; unresolved cases are refused without an exact/subdivision
fallback. See `PRIMITIVES.md` and the current scaling evidence in `../VALIDATION.md`.

## rc.5 primitive element choices

Tet4 and straight-sided Tet10 now work for all six explicit primitive families.
Graded tetrahedral refinement improves the standard cone-tip case; the swept
Hex/Pyramid route and extremely sharp or thin shapes can still have poor cells.
Quality warnings preserve valid meshes and identify affected cell ordinals;
strict rejection remains available. These are bounded, tested primitive routes,
not a general sharp-feature optimizer. Curved Tet10 primitive placement,
adaptive local transitions for hex-dominant cones, and imported CAD sharp-edge
recovery remain development work. See `ELEMENT_CHOICES.md`.
