# Changelog

## Unreleased

- Add a Normal resource policy with demand-driven construction storage and no
  practical logical-work ceiling when `max_work` is omitted. Explicit work,
  memory, count, output and cancellation limits remain authoritative. Bounded
  admission remains available and old recipes preserve its meaning.
- Default new Python and STEP CLI requests to Normal admission. Preserve
  Bounded C++ aggregates and old recipes, Legacy construction effort and the
  inherited 512 MiB modeled working allowance.
- Admit the actual initial numerical mesh before growing its live storage;
  account for queue growth, compaction and overlapping temporary allocations.
  Separate construction stages that do not coexist in memory.
- Replace Normal curve and constraint-source count cutoffs with admitted
  resource limits; retain native boundary, material and publication checks.
- Record effective resource settings and completed cancelled-stage work. Add
  resource-policy parity to installed-wheel qualification.
- Add opt-in, bounded Tet4 quality inspection, native-ID pagination, spatial
  distributions, trimmed CAD centroid distances and quality-field VTK export.
- Add deterministic SVG mesh views with whole-element cuts and visible interior
  highlights, plus common-region comparisons of reported-valid mesh snapshots.
- Expose Python meshing signatures, guarded selection sizing and native hints;
  preserve the distinction between requested policy and native body routes.
- Expose bounded surface/interior construction caps in C++, Python, recipes and
  the STEP CLI, with per-stage limit and termination receipts. Keep established
  allowances as the default; opt-in adaptive policies reserve interior effort
  for bulk placement and remain subordinate to overall resource limits.
- Report construction-limit warnings separately from geometric validity,
  achieved sizing and element quality. Preserve historical recipe defaults.
- Retain native validity gates and archive v1. These development changes are
  not part of the published RC3 wheels.

## 1.2.0-rc.3

- Add capability discovery, paginated CAD inspection, saved-recipe execution and
  shared human/JSON diagnostic summaries without changing native meshing checks.
- Add CPython 3.12–3.15 wheels on Windows and Linux x64/ARM64 and macOS ARM64,
  with macOS 15/26/27 installation checks and reusable verified CAD SDKs.
  Require the complete builder and installation matrix before publication.
- Preserve Unicode CAD names in JSON output on terminals with restricted encodings.
- Select the validated MSVC 14.44 toolset for Windows ARM64 and package the actual
  runtime dependency closure. Record native compiler and configuration evidence.
- Verify restored CAD SDK contents and create Linux cache directories on the host
  before container builds, so a cache hit can reuse the validated files.
- Keep native meshing algorithms, validity checks and mesh archive format unchanged.

## 1.2.0-rc.2

- Construct initial native Tet4 wall sites before filling interior gaps, with
  automatic selection for compatible single-solid requests and explicit on/off
  controls in C++, Python and the STEP command line.
- Reuse source-contact classifications and certify complete local cavities
  without repeating all cone-pair comparisons; retain exact admission checks.
- Record bounded automatic source-stage fallback separately from successful
  initial construction, including its failed work and elapsed time.
- Prepare three-platform PyPI publishing with matching source and qualification
  receipts, clean installed consumers and a fresh package-index installation check.
- Add an optional Python interface for native STEP import, audited Tet4 meshing,
  assembly interference reports and tagged VTK/archive export.
- Build and test a Windows x64 CPython 3.12 wheel with the pinned CAD runtime.
- Expose Python CAD face and box sizing, minimum field size and growth, with
  per-control overrides, request diagnostics and command-line controls.
- Add optional NumPy mesh snapshots with typed connectivity, CAD classifications
  and an extraction byte budget; preserve native mesh ownership and ordering.
- Build macOS 15 ARM64 and manylinux 2.28 x64 wheels alongside Windows, with
  repaired dependencies, native CLI parity and independent artifact installs.
- Inspect native CAD topology and geometry in Python, select faces by explicit
  predicates or representation names, and save guarded, source-bound groups.
- Add tagged CAD display previews and bounded, source-bound meshing recipes,
  with explicit diagnostic retries and atomic publication of project bundles.
- Expose optional Python quality warning and strict thresholds, reporting their
  effective policy and preserving the requested values through recipe retries.
- Preserve Unicode file paths and full VTK coordinate precision, export both
  classification dimensions and tags, and publish VTK files atomically.

## 1.2.0-rc.1

- Combine native CAD face and box sizing in the single-solid Tet4 route.
- Save and reopen STEP sizing recipes and verify face and body groups in C3D4.
- Build pinned OCCT on Windows, Linux and macOS and check the cylinder fixture.
- Support Apple SDKs that gate thread cancellation or lack floating-point `from_chars`.
- Document the native STEP preview, solver examples and candidate validation.

## 1.1.0

- Add asynchronous meshing jobs, a C and Python bridge and a local UI example.
- See the [1.1.0 release notes](docs/RELEASE_NOTES_1.1.0.md).

## 1.0.0 — 2026-10-02 UTC

- First official 1.0.0 source release, retaining the tested meshing implementation.
- Finalize source, CMake/CPack, SDK and CLI version identity as 1.0.0.
- Replace active candidate documentation with release notes, scope and compatibility policy.
- Add version-consistency checks to source packaging and an external ZIP SHA256 file.
- Preserve original qualification logs and historical version labels as evidence.
- Keep C++26/GCC 16, uniform defaults, archive v1, legacy C symbols and exact-version CMake selection.
- Defer external-mesher comparisons and additional OCCT work to the next development stage.


## 1.0.0-rc.11

- Add immutable composable spatial size fields: finite point/segment/triangle,
  box/ball distance, linear/smooth Threshold, Min/Max, Scale and Clamp. Validate
  graph/envelope/budgets; prune unreachable nodes; allocate nothing during queries.
- Share conservative whole-box sizing through BoundedSizeField and preserve
  RegionSizeField behavior. Add deterministic adaptive query-work accounting.
- Integrate global/per-body fields in the constrained geometry Tet4 route with
  explicit incompatible-control rejection and optional bounded interior refinement.
- Add independent linear-volume size auditing: conservative and sampled ratios,
  warning/strict policies, aggregate work budgets and callback-failure diagnostics.
- Add field oracles, scaled/rotated composition tests, dense cell-interior checks,
  rollback/fault injection, an installed SDK example and archive-backed visuals.
- Preserve uniform defaults and archive v1. Rebuild C++ clients after public
  layout/base-class changes. General protected-boundary refinement remains limited.

## 1.0.0-rc.10

- Reuse capped axial profiles in recognized polygonal geometry sweeps and mapped blocks;
  retain source-face selection, face/node/body classification and uniform defaults.
- Add Native/Tet4 body output and public convex Tet4/Prism6/Hex8-to-Tet4 conversion.
  Shared node-ID diagonals preserve conformity without inserting nodes; validate
  convexity, planar faces, exterior winding, positive tets, budgets and rollback.
- Accept explicit normalized mapped-axis coordinates in ordinary/workspace paths.
- Add optional geometry quality warnings and constrained-boundary failure details.
- Fix decimal layer/depth roundoff without permitting meaningful wall overruns.
- Add complex engineering profiles, converter permutations/scales, fault injection,
  non-sweep admission probes, installed example, benchmark scripts and exact mesh cuts.


## 1.0.0-rc.9

- Added immutable capped size profiles with bounded exact piecewise-linear
  lower-envelope compilation, overlapping controls, finite axial bands and
  selected cylinder/tube caps. Distribution independently checks whole-interval
  size bounds and actual adjacent axial growth, with plateau compaction.
- Reused those distributions for cylinder/tube Tet4/Tet10/Hex8 sweeps and
  prepared Tri3/Quad4/mixed prism/hex layers.
- Added opt-in varying-section solid-cylinder Tet4/Tet10 meshes and generic
  conforming Tet4-to-Hex8 conversion with shared nodes, semantic classification,
  quality policies, complete Jacobian certificates and boundary audits.
- Reuse consecutive identical slab sections to reduce Delaunay work. Added an
  uncached comparison switch and built/reused slab diagnostics. Resource
  admission remains conservative with reuse enabled.
- Added independent field, topology, geometric-volume and high-order checks,
  rollback and allocation-failure coverage, benchmarks, installed SDK example,
  reproducible archive-based visuals and a detailed sizing contract.
- C++ request/report layouts changed: rebuild consumers. Defaults and the mesh
  archive format remain compatible. This is still a release candidate.

# Changes

## 1.0.0-rc.8 — local sizing and planar wall layers

- Add analytic distance-based Box/Ball/Rectangle size controls, finite face
  selection, minimum-size floor, per-region growth and conservative box minima.
- Add bounded adaptive box Tet4 generation with strong 2:1 octree balance,
  shared interface triangulations, all-edge size caps, quality/growth acceptance
  gates, deterministic output and transactional publication.
- Add optional adjacent-volume growth measurement to the existing bucketed
  boundary audit without a second adjacency build. The default audit remains
  on its original compact path.
- Add explicit mixed-sweep layer fractions, prescribed first-wall thickness,
  geometric wall layers and a capped graded interior schedule. Planar
  structured/unstructured Tri3 patches produce Prism6/Wedge6; convex Quad4
  patches produce Hex8, including conforming mixed patches.
- Validate source intersections/T-junctions and unused nodes, retain complete
  boundary auditing, and test invalid input, budget limits and allocation faults.
- Include actual geometric section cuts, reproducible benchmark scripts and
  an installed SDK consumer. C++26/GCC 16 remains the reference.
- New box refinement currently produces Tet4. Curved layer fronts, independent
  tetrahedral-core attachment, spherical outer domains and automatic CAD-owned
  refinement remain outside this new route's qualification.
- Public C++ request/report structs have new fields: rebuild C++ consumers.
  The mesh archive format remains version 1.

## 1.0.0-rc.7 — balanced circular sweeps and optional surface features

- Extend balanced solid-disk Tet4/Tet10 generation to solid tori and elbows;
  share the torus seam, preflight topology-specific quadratic counts, and retain
  full validation, rollback, quality policies and explicit legacy selection.
- Add `SurfaceFeatures.hpp`: default-off, bounded, deterministic read-only
  Tri3/Tri6/Quad4 feature analysis with angle/tolerance/length controls, optional
  tag seams, topology diagnostics, cancellation and exact output allocation.
- Add CLI controls and counts, installed SDK example, independent circular
  quality and feature tests, allocation fault cases, actual cutaways, and
  repeat/thread/size scaling measurements.
- Preserve cylinder meshes and routes not selected for this change. C++26,
  GCC 16 and the existing strict floating-point contract remain the reference.
- Hollow/tapered uniform layouts, local feature-driven refinement, curved
  general CAD import, projected primitive Tet10, and cross-platform/GPU
  qualification remain future work. See VALIDATION.md for measured scope.


## 1.0.0-rc.6 — uniform cylinder tetrahedra

- Direct concentric triangular sections and balanced axial spacing for cylinder Tet4/Tet10; automatic selection, explicit Uniform/Swept controls, geometry-limited sizing diagnostics.
- Preserved exact orientation checks, full boundary auditing, certified straight Tet10, transactional publication, deterministic threading and checked memory budgets.
- Independent uniformity/shape tests, thin-disk warnings, strategy/CLI checks, scaling comparisons and actual mesh visuals.
- Structural/CFD roadmap with measurable acceptance gates and primary references.
- API source structs extended: rebuild C++ consumers; default cylinder tetrahedral connectivity changes. Select Swept to retain the previous route. Archive format unchanged.

## 1.0.0-rc.3

- Replace global face sorting with exact minimum-node buckets, compact temporary
  indices, checked temporary-memory budgets, cancellation and allocation diagnostics.
- Retain the rc.2 audit in tests and compare all diagnostics on mixed cells,
  renumbered/reordered meshes, cyclic quads, high-valence nodes and seeded damage.
- Add `GeometryBodyControl`: target size/strategy per solid, explicit source face,
  exact sweep layers and geometric grading towards either end.
- Support decreasing prepared-sweep widths with a symmetric adjacent-growth bound.
- Move direct mapped/sweep resource checks ahead of schedule allocation.
- Fall back to emitted-cell checks when ideal construction is too close to
  coordinate-resolution limits; reject nonpositive constituent sweep tetrahedra.
- Expand engineering-profile, STEP/BREP and graded CAD qualification, with
  reproducible reference-versus-new timing and process-RSS measurements.

C++26/GCC 16 remain the reference. Rebuild C++ clients: this candidate extends
request/report structures and does not promise ABI compatibility with rc.2.

## 1.0.0-rc.2

- Retain C++26/GCC 16 and add geometry-driven Hex8, straight Prism6 sweep,
  and constrained Tet4 selection, with per-body diagnostics and atomic output.
- Add an immutable planar polyhedral backend with loop, hole, shell and volume
  validation, and exact straight-edge projection.
- Add optional Open CASCADE 8+ single-solid planar STEP/BREP import and
  `paacs-cad-mesh`; refuse unsupported curved surfaces and assemblies.
- Add bounded deterministic planar constraint insertion and hole-loop support
  in the production surface bridge; preserve side/cap entity classifications.
- Add independent mixed linear boundary incidence audits, explicit limits,
  cancellation checkpoints and optional minimum element-quality enforcement.
- Fix undersized implicit-target sweep classification arrays, inverted-cell
  quality acceptance, and surface-recovery refusals misreported as allocation failures.
- Add geometry, CAD, failure-injection and seeded planar-recovery qualification,
  independent analytic-volume checks, thread determinism and isolated RSS/timing tooling.

This remains a release candidate. New virtual capabilities require recompiling
C++ clients. Archive format version 1 is unchanged; no general CAD or competitor
performance claim is made.

## 1.0.0-rc.1

### Toolchain and packaging

- C++26 required by default; GCC 16 reference preset and explicit C++23 mode.
- CMake 4.4.3 reference, 3.30 minimum; root build, explicit source list,
  compile commands, optional IPO/native tuning, ASan/UBSan and TSan presets.
- Version-independent `paacs-mesh` CLI and `PAACSMesh` installable package.
- Full inherited unit suite, retained probes, new failure-path regressions,
  installed C/C++ consumer examples and future GitHub CI workflow.
- Strict floating-point compilation for predicates and inline geometry.

### Correctness

- Avoid integer overflow in parallel chunk counts, worker narrowing and endpoints;
  define zero chunk size consistently in map-reduce.
- Make the global compatibility thread setting atomic; restore pool context after
  normal returns and exceptions.
- Preserve mesh node/classification and cell/connectivity consistency on allocation
  failure; make dedup insertion and rebuild transactional for logical mesh state.
- Reject invalid cell types, invalid dedup tolerances, non-finite spatial inputs,
  overflowing C API tags, and non-finite primitive geometry.
- Parse CLI numeric arguments completely and reject invalid/negative unsigned values.
- Validate archive header/payload/counts/reserved fields/connectivity before mesh
  allocation; expose configurable reader budgets and preserve outputs on failure.
- Detect delayed file write errors at stream close.
- Repair an inherited lattice test's instrumentation assumption and refresh two
  stale work-count goldens after reproducing them against the unmodified input;
  retain the geometric, quality, and mesh-identity checks.

### Efficiency and maintainability

- Span and initializer-list connectivity APIs avoid temporary vectors in existing
  brace-list cell insertion sites. Self-referential connectivity spans are safe.
- Allocate archive buffers to exact size, use endian-aware fixed-width stores,
  and reuse bounded decode connectivity storage instead of allocating per cell.
- Use a tested small integer sorting network for three/four-node face keys.
- PointGrid size() is constant-time; extreme-radius queries scan occupied buckets
  rather than enumerating unbounded empty cell ranges.
- Reformat compact multi-statement implementations and remove GCC 16 diagnostics
  without disabling the warning categories.

### Compatibility and scope

- Archive format version 1 and existing C symbol names are retained.
- C++ consumers must rebuild; final stable API/ABI and licensing remain decisions.
- Default archive-read budgets are new; explicitly increase them for larger files.
- The release scope covers the documented native core and its validated routes.
- Validation and measured performance are recorded in VALIDATION.md.
