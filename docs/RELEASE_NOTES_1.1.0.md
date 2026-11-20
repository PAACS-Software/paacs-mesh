# PAACS Mesh 1.1.0

Release date: 2026-10-02 UTC. C++26 / GCC 16 reference build.

This release adds a frontend-independent application layer. Meshing kernels,
numerical predicates, default mesh strategies and archive v1 remain unchanged.

## Added

- Bounded asynchronous `sdk::JobService` with polling, cancellation,
  serialized engine dispatch, structured refusal reports and immutable results.
- Owning adapters for primitives, adaptive boxes, geometry models, layered
  patches and mapped hexes; copied profiles/coordinate spans and retained
  immutable model/field/source ownership.
- Optional `PAACS::ui` shared library with C ABI v1, primitive/adaptive request
  discovery, 68 scalar controls, profile/field arrays and paged result reads.
- Dependency-free Python binding and JSON Lines protocol v1; explicit handle
  release, exact uint64 identifiers, archive export and strict input handling.
- Local reference web UI with discovered controls, recipe load/save, background
  jobs, cancellation, settings-revision checks and tagged surface previews.
- Per-job and retained-result logical byte budgets. Strict hashing is opt-in;
  disabling the extra hash pass never disables mesh validity checks.
- Concurrency/lifetime/cancellation/ownership tests, C/JSON contract tests, native
  output comparisons, installed examples and a reproducible job-overhead probe.

## Scope

The five native adapters preserve complete native options and reports. The
portable bridge currently covers primitives and adaptive boxes; it does not
claim exhaustive C/JSON coverage of geometry, layers, hybrid/prepared Tet or
feature analysis. Those native operations remain available to application code.
See [UI integration](UI_INTEGRATION.md) for contracts, recipes, limits and
the specific next steps before broader OCCT integration.

The reference UI is replaceable, and no frontend framework is required by the
library. It is a local example, not a deployed application. Cross-platform design
does not substitute for Windows/macOS qualification. Existing limitations in
general curved CAD, boundary refinement and curved CFD layers remain.

[Validation](../VALIDATION.md) separates fresh tests from inherited qualification.
Comparative superiority over external meshers has not been established.
The [1.0.0 release notes](RELEASE_NOTES_1.0.0.md) and original evidence remain
available. Distribution licensing is still recorded in [LICENSING_STATUS.md](../LICENSING_STATUS.md).

The [meshing roadmap](ROADMAP.md) records current element support and the
qualification milestones for broader native CAD, Prism6 wall layers, curved
Tet10, Hex8 and Hex27, local mesh editing, and structural and flow studies.
The [STEP cantilever example](../examples/elasticity_cantilever/README.md) includes
the CAD fixture, audited Tet4 meshes reproduced on four size levels, and a
three-dimensional elasticity solve with displacement and equilibrium checks.
The [STEP plate with a hole](../examples/elasticity_plate_hole/README.md) adds a
paired local-refinement and stress study. The [release status](RELEASE_READINESS.md)
records the latest fixed STEP corpus results and the remaining coverage gaps.

### MSH format API names

C++ callers use `ElementTraits::mshId`, `T5ExportFormat::Msh41` and
`writeT5MshExport` in the current source. Update those references and
rebuild and relink C++ consumers. The `writeMsh` entry point, MSH 4.1
wire format and serialized format IDs retain their existing behavior.
