# Local qualification — 1.0.0-rc.3

Date: 2026-09-29. C++26 / GCC 16 remain the reference. This candidate improves
boundary-audit cost and control of verified block/sweep routes.
The prior report is preserved in [docs/VALIDATION_RC2.md](VALIDATION_RC2.md).
Current raw evidence is under `evidence/rc3/`.

## Executed results

| Check | Result and scope |
|---|---|
| Release build | Core, OCCT adapter, CLIs, tests and retained benchmarks built with GCC 16.2.0 / C++26; production warnings are errors |
| Complete release CTest | **18/18 passed**, including **all 475 inherited cases**; 308.72 seconds |
| Geometry corpus | **42/42** published; positive independently computed cell volumes, analytic total volumes, classification/trim/face coverage checks, valid boundary incidence, identical meshes at 1/4 threads |
| Retained geometry output | All **32 rc.2 mesh hashes unchanged** |
| Audit differential/resource suite | **2,015 cases passed**, including 2,000 seeded mutation/renumbering trials, mixed element types, crossed quads, high-valence buckets, invalid rows and exact budget boundaries; all old diagnostic counters matched the reference |
| New body/layer/numerical suite | **58 cases passed**: grading in both directions, reversed caps, exact layer counts, per-body size/strategy, refusal rollback, huge schedule preflight, coordinate collapse and folded prism rejection |
| Existing geometry hardening | **56 checks passed** |
| Planar recovery | **40 seeded radial PSLGs + 3 refusals**, 3,498 flips, deterministic replay and work-budget rollback passed |
| CAD round trips | **32 STEP/BREP files**: 26 successful planar meshes and 6 expected curved/assembly refusals; **6 additional graded/reversed CAD cases** passed |
| Allocation failures | 17 retained SDK cases, **178 model + 242 pipeline allocation points**, plus separate checks for both audit allocations passed |
| ASan + UBSan | **6/6 selected suites passed**, 27.95 seconds: geometry corpus, differential audit, hardening, planar recovery, allocation failure and body/layer controls |
| Installed SDK | **5/5 consumers passed**: C, C++, geometry, installed CAD target and new graded-layer API |
| CAD CLI | **20/20 cases passed**; every failed request preserved existing output bytes |
| Corpus measurements | **264/264 fresh-process runs**, 42 cases plus 2 dense block sizes, 1/4 threads, 3 repetitions; stable hashes throughout |
| Source archive | Fresh extracted source configured successfully with the reference tools and OCCT; ZIP CRC and every manifest checksum verified |
| Paired audit measurements | **60/60 fresh-process runs**, 3 sizes × 2 numberings × 2 algorithms × 5 repetitions; matching face/count diagnostics |

New engineering profiles include I/T beams, a thin channel, stepped bracket,
slotted plate, thin-walled duct, nine-hole plate, polygonal tube, near-wall hole
and acute triangle. These are controlled polygonal fixtures, not a general
industrial CAD dataset. A polygonal tube is not a curved-cylinder capability.

## Audit runtime and temporary storage

Both algorithms were compiled into the same executable with the same flags.
The reference preserves rc.2's global sort; the new audit groups by minimum
node and sorts complete cyclic keys inside each group. No validation is
silently disabled. Records use compact temporary node IDs where possible;
the published mesh keeps 64-bit indices.

Medians of five fresh processes; MB in this table means 1,000,000 bytes.

| Hex8 volume cells | Numbering | rc.2 audit ms | rc.3 audit ms | Speedup | rc.2 temporary MB | rc.3 temporary MB |
|---:|---|---:|---:|---:|---:|---:|
| 15,625 | Natural | 10.523 | 2.469 | 4.26x | 3.90 | 1.70 |
| 15,625 | Shuffled | 12.790 | 4.622 | 2.77x | 3.90 | 1.70 |
| 250,047 | Natural | 175.886 | 37.990 | 4.63x | 60.96 | 26.48 |
| 250,047 | Shuffled | 227.617 | 89.751 | 2.54x | 60.96 | 26.48 |
| 1,000,000 | Natural | 714.809 | 146.222 | 4.89x | 242.40 | 105.20 |
| 1,000,000 | Shuffled | 956.670 | 434.194 | 2.20x | 242.40 | 105.20 |

On the 250,047-cell naturally numbered block, the audit is **4.63x faster** and
uses **56.6% less temporary array storage**. On one million cells it is **4.89x
faster**; randomly renumbering those same nodes reduces the gain to **2.20x**.
Renumbering also moves coordinates/classifications consistently. The input
geometry and connectivity relationships remain unchanged.

For the naturally numbered 250,047-cell benchmark, maximum observed whole
process peak RSS fell from **90.86 MiB to 57.88 MiB**. This RSS includes mesh
construction; shuffled runs additionally include renumbering buffers. The
reported temporary bytes cover requested arrays, excluding allocator
bookkeeping and the input mesh. They are not a process-memory cap. Small RSS
measurements have an approximately 12 MiB process-launch floor.

The new audit remains O(N + F + sum b_i log b_i), apart from the linear cell
preflight. A high-valence bucket can still require O(F log F) work. The
12,000-use hub fixture checks correctness of this case. Cancellation polls
scans and bucket boundaries; an individual bucket sort is not interruptible.
See [audit design](AUDIT_DESIGN.md) for exact guarantees.

## Full geometry qualification cost

The dense unit block contains 250,047 volume cells, 262,144 nodes and 273,861
cells including boundary quads. The following are medians of three fresh
processes; peak RSS is the maximum of the three. MiB means 1,048,576 bytes.

| Threads | Generation ms | Boundary audit ms | Whole qualification ms | Logical output MiB | Process peak MiB |
|---:|---:|---:|---:|---:|---:|
| 1 | 20.409 | 45.439 | 218.474 | 30.26 | 62.75 |
| 4 | 13.950 | 38.222 | 209.616 | 30.26 | 62.93 |

Whole qualification includes geometry construction, meshing, independent
volume/classification checks, quality analysis and archive hashing. It is
neither generation-only time nor only the `meshGeometryModel` call. Stage
numbers differ from the isolated audit microbenchmark because their preceding
work and allocation/cache state differ. Do not present the audit's 4.63x
speedup as the speedup of the entire program. The boundary audit is currently
serial; the thread setting affects eligible generation kernels.

## Robustness defects fixed

1. Global face sorting was avoidable. Exact minimum-node groups preserve quad
   cycles, ownership counts and interior orientation while using smaller
   records. Byte limits are checked before both temporary allocations;
   allocation failure and cancellation are explicit outcomes.
2. Layer schedules could allocate before direct mapped/sweep count budgets
   were checked. Huge-count refusal now precedes schedule allocation.
3. An ideal affine frame could certify cells that collapsed after coordinates
   were rounded far from the origin. The reproduced block at z=1e15 with 12
   layers was accepted before the fix. Conservative resolution guards now
   request emitted-cell checks for mapped blocks and straight sweeps; the
   regression refuses the collapsed meshes and leaves the destination intact.
4. A folded prism could have a positive total signed volume while one
   constituent tetrahedron was inverted. Prepared sweep checks now require
   every fixed-decomposition tetrahedron to have positive finite volume.
5. Prepared-sweep connectivity row offsets are checked before pointer access.

These are targeted safeguards, not a general geometric intersection or
continuous-Jacobian certificate for arbitrary warped elements. The separate
quality and topology contracts remain relevant.

## Controls and practical scope

`GeometryBodyControl` selects size and strategy per solid and, for verified
blocks/sweeps, a starting face, exact layer count and last/first layer-width
ratio. Ratios below one refine towards the target. Adjacent growth is bounded
symmetrically; impossible, invalid or numerically collapsed requests refuse
publication. The starting face can reverse a sweep. These are whole-body
controls, not automatic decomposition of arbitrary solids or general
anisotropic Tet meshing. [The guide](GEOMETRY_GUIDE.md) and
[compiled example](../examples/graded_geometry.cpp) show the API and CLI.

The pyramid fixture still has minimum SubTetMeanRatio approximately
**3.71e-10**. A valid boundary and correct total volume do not make that mesh
solver-ready. Optional quality floors can refuse output; guaranteed sliver
removal remains unimplemented. Thin valid cells can also have poor isotropic
quality. Full boundary auditing stays enabled by default.

CAD qualification still covers single valid planar, straight-edge solids.
Curved/periodic CAD, healing, assemblies, general solid decomposition and
conformal interfaces between touching bodies remain outside this route.
Separated bodies require disjoint bounding boxes. The polyhedral builder
does not prove absence of arbitrary inter-face intersections. Input/output
limits are not global process-RSS bounds; the OCCT parser is not time/memory
sandboxed. Archive writes retain the existing non-atomic file-write contract.

The optional OCCT binaries and all 475 inherited cases were not sanitizer-run
in this pass. Leak detection remained disabled for the known environment
restriction; leak freedom is not established. An initial concurrent sanitizer
compilation exceeded the 8 GiB workspace limit and its compiler was killed;
a serial rebuild succeeded before the six sanitizer suites passed. Both logs
are retained. This was a compilation resource failure, not a mesh-test failure.

No Windows/macOS/Clang qualification, GPU acceleration, new TSan result,
external industrial corpus, competitor parity, stable ABI or final 1.0
readiness is claimed. The one-external-dispatcher worker-pool contract remains.
Recompile C++ consumers for the extended request/report structures.

## Reproduction and contents

Linux x86-64, AMD EPYC 9V74 host, 8 CPU-equivalent quota, 8 GiB memory limit.
GCC 16.2.0, C++26, CMake 4.4.3, Ninja 1.13.2 and optional OCCT 8.0.1.
Release uses `-O3 -DNDEBUG -std=c++26 -fno-fast-math -ffp-contract=off`; IPO and
native-architecture flags are off. Timings are shared-host observations.
Build paths retain the working folder's rc.1 name; release metadata is rc.3.
The final 64-bit audit build was also rerun through the differential suite
after relaxing an unnecessary size assertion for 32-bit ABI compatibility.
Actual 32-bit execution remains unqualified.

```sh
ctest --test-dir build/gcc16-release --output-on-failure
python3 benchmarks/run_audit_comparison.py \
  build/gcc16-release/tests/paacs-audit-benchmark --output audit-results
python3 benchmarks/run_geometry_corpus.py \
  build/gcc16-release/tests/paacs-geometry-corpus --output geometry-results
python3 tests/cad_cli_checks.py build/gcc16-release/paacs-cad-mesh
```

Build commands are in [docs/BUILDING.md](BUILDING.md). Raw evidence includes
all qualification logs, JSON/CSV timings, compiler commands, the environment,
and unchanged-mesh comparisons. [The gallery](mesh_gallery_rc3.png) renders
actual CLI-generated boundary meshes from retained STEP fixtures. All **68
production C++ translation units** remain present, including all 63 original
ones. The ZIP contains full source and evidence; compiler and OCCT binaries
are not bundled.
