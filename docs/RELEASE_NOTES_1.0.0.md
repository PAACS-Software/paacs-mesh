# PAACS Mesh 1.0.0

Release date: 2026-10-02 UTC.

This is the first official source release of the current PAACS meshing library.
It finalizes the tested implementation and documented controls with the C++26 /
GCC 16 reference toolchain. It is a version and packaging finalization; no meshing
algorithm, numerical predicate, element-quality policy or default strategy was
changed during this promotion.

## Included

- Exact-predicate tetrahedral infrastructure and bounded constrained meshing.
- Uniform primitive Tet4/straight Tet10 choices and mapped/swept Hex8/Prism6
  routes, with documented pyramid transitions and conforming element conversions.
- Reusable spatial and axial sizing, region/surface influences, growth controls,
  optional feature analysis and prepared planar wall layers.
- Geometry ownership, boundary/quality/size diagnostics, explicit budgets,
  cancellation and transactional publication on documented routes.
- Source, tests, examples, actual mesh figures and retained qualification evidence.

## Finalization changes

- Source, SDK, CLI, CMake package and CPack metadata identify version **1.0.0**.
- Current landing/readiness documents and control guide headings use the final
  version; historical test evidence keeps its original labels and measurements.
- Release scope and compatibility policy are explicit. CMake continues to require
  exact-version selection; archive v1 and the existing C symbol names are retained.
- Packaging verifies version agreement and produces a deterministic source ZIP,
  an internal per-file manifest and an external SHA256 checksum.

See [VALIDATION.md](../VALIDATION.md) for completed checks and
[docs/RELEASE_POLICY.md](../docs/RELEASE_POLICY.md) for the supported scope.

## Next stage

Matched comparisons with other meshers and additional OCCT work are deferred.
They are not prerequisites for assigning this source version, but the present
evidence does not establish universal CAD robustness or comparative superiority.
Curved CAD, coupled boundary/interior refinement, periodic/embedded constraints,
curved CFD layers, broader platform qualification and acceleration remain on the
[roadmap](../docs/GENERAL_MESHER_ROADMAP.md).

Distribution licensing remains unspecified; see [LICENSING_STATUS.md](../LICENSING_STATUS.md).
