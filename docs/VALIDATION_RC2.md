# Local qualification — 1.0.0-rc.2

Date: 2026-09-29. C++26 remains the default. This is a tested release candidate,
not a stable general-purpose CAD release. The previous modernization report is
preserved in [docs/VALIDATION_RC1.md](VALIDATION_RC1.md); its raw evidence
remains under `evidence/`. New evidence is under `evidence/rc2/`.

## Executed results

| Check | Result and exact scope |
|---|---|
| GCC 16 / C++26 full release build | Core, optional OCCT adapter, both CLIs, tests and retained benchmarks built; production libraries use warnings as errors |
| Complete CTest | **16/16 passed**, including **all 475 inherited unit cases**; 277.81 seconds |
| Geometry corpus | **32/32 shapes published**, each with matching output at 1 and 4 threads, independent analytic volume, positive oriented cell volumes, boundary incidence, entity position/trim checks and complete face coverage |
| CAD files | **24 round trips**: 18 successful planar STEP/BREP meshes and 6 expected refusals for cylinders, spheres and assemblies; file-size, malformed-file, missing-file and format checks also passed |
| Geometry hardening | **56 checks passed**: invalid shells/loops, nonfinite geometry, count/work/byte budgets, quality thresholds, cancellation, forced-strategy mismatch, boundary defects, classification arrays, separated-body strategy selection and enumeration determinism |
| Planar constraint recovery | **40 seeded radial PSLGs**, 3 refusal cases, deterministic replay and budget rollback passed; 3,498 edge flips exercised |
| Allocation failure | Retained 17 SDK injection failures passed; new builder and mapped-meshing rollback passed across **178 + 241 allocation points** before first success |
| ASan + UBSan | **4/4 selected CTest entries passed**: the complete new geometry corpus, geometry hardening, planar recovery and allocation failure; 16.83 seconds |
| Installed SDK | **4/4 consumers passed**: C, C++, new geometry API and optional installed CAD target |
| CAD CLI | **8 publication/refusal scenarios passed**; failed requests preserved an existing output archive |
| Isolated measurements | **204/204 runs passed**: 32 cases plus 2 dense block sizes, 1/4 threads, 3 fresh processes per combination; hashes stable across repetitions and thread counts |

The 32 shapes include axis-aligned, rotated, skewed, translated, small/large,
thin and slender blocks; triangular and regular polygonal prisms; L/U/C-shaped
profiles; a narrow ligament; one/two holes; oblique sweeps; a tetrahedron,
pyramid and frustum; and explicit sweep/Tet overrides. The original five
surface-recovery failures at the same target sizes are retained in
`geometry-before-recovery.jsonl`. They now pass without relaxing the topology,
volume or classification assertions.

## What the tests found and changed

1. The production surface bridge exposed only the outer loop. It now carries
   all ordered holes into the existing audited surface adapter.
2. Some concave/multiple-hole profiles required boundary recovery beyond the
   existing no-Steiner route. A bounded deterministic edge-flip pass now inserts
   literal constraints before the audited publication adapter. It preserves
   source points and constraint endpoints. It is a constrained triangulation,
   not a guarantee of a constrained-Delaunay or optimal-quality triangulation.
3. An implicit-translation sweep could read a short target-classification
   array out of bounds. Array lengths and classifications are checked before
   access, including the new intermediate-node and side-edge maps.
4. Mixed quality evaluation converted negative tetrahedral volume to positive,
   allowing an inverted cell to appear valid. It now requires a positive exact
   orientation and positive finite signed volume. A reversed-tet regression
   pins this behavior.
5. An unattempted suffix after a geometric recovery refusal was reported as an
   allocation failure. The first actual failing outcome is preserved now.
6. The planar geometry backend uses direct analytic line projection; it avoids
   iterative nearest-point approximation for straight edges.

## Measured runtime and memory

| Geometry | Threads | Nodes | Volume cells | Generation ms | Boundary audit ms | Logical mesh MiB | Peak process MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| concave_l | 1 | 235 | 264 | 2.968 | 0.155 | 0.037 | 12.00 |
| holes_two | 1 | 755 | 936 | 12.705 | 0.581 | 0.124 | 12.00 |
| pyramid | 1 | 56 | 116 | 72.658 | 0.056 | 0.014 | 12.00 |
| block_forced_tet | 1 | 80 | 200 | 68.156 | 0.121 | 0.022 | 12.00 |
| block_unit (dense) | 1 | 17,576 | 15,625 | 1.150 | 10.495 | 2.048 | 12.00 |
| block_unit (dense) | 1 | 262,144 | 250,047 | 18.966 | 172.540 | 30.257 | 91.88 |
| block_unit (dense) | 4 | 262,144 | 250,047 | 13.333 | 172.079 | 30.257 | 91.91 |


Numbers are medians of 3 fresh processes; RSS is the maximum observed peak
among those 3. `volume cells` excludes boundary/curve cells. Logical storage
includes nodes, classifications, all cells and connectivity, and excludes
container capacity and workspaces. Peak RSS covers the whole qualification
process: geometry construction, meshing, independent checks, quality analysis
and archive hashing. A fresh Python wrapper obtains child `ru_maxrss`; this
includes process-launch overhead, with a roughly 12 MiB floor in small runs.
It is not an isolated allocator measurement or a guaranteed mesher memory cap.

The 250,047-Hex8 case had 262,144 nodes and 30.26 MiB
logical output. Four-thread generation measured 13.33 ms,
versus 18.97 ms on one thread
(1.42x for generation only).
Its independent boundary audit added 172.08 ms. The full
qualification process took a median 339.91 ms and peaked
at 91.91 MiB. The small 15,625-cell block did not
benefit from four threads, illustrating thread-dispatch overhead. Auditing is
the main cost at the larger block size; generation-only throughput must not be
presented as full-pipeline throughput.

These are different cell/quality workloads, not equal-quality comparisons
between strategies.
A prior measurement pass lacked RSS because `/usr/bin/time` was unavailable;
the retained final matrix was rerun with the fresh-process resource wrapper.

## Quality and support limits

The pyramid fixture has a minimum SubTetMeanRatio of approximately
**3.71e-10**, despite a valid boundary and correct total volume. This is a
real quality limitation, not a broad simulation-readiness pass. The new
`minimumElementQuality` option can refuse such output before publication;
it does not implement guaranteed sliver removal. Thin/slender geometries may
also produce low isotropic quality when the requested size is coarse.

The CAD importer is qualified for one valid solid with planar faces and
straight edges. The fixtures are controlled models generated independently by
OCCT; they are not an industrial customer corpus. Curved/periodic CAD, healing,
assemblies, automatic decomposition within arbitrary solids, and automatic
conformal interfaces between touching bodies remain outside this qualification.
Separated bodies can select different strategies. The polyhedral builder does
not certify arbitrary inter-face self-intersections.

No Windows/macOS/Clang run, GPU implementation, new full-pipeline TSan run,
competitor comparison, or stable C++ ABI is claimed. The prior focused TSan
parallel-helper evidence remains historical. The single external dispatcher
contract remains. C++ clients must be rebuilt; archive format 1 is unchanged.

ASan/UBSan used `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1` and
`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`. The environment's earlier
LeakSanitizer restriction remains; leak freedom is not established. The optional
OCCT binaries and the full inherited 475-case suite were not sanitizer-qualified
in this pass. Cancellation is at analysis/stage checkpoints, and output byte
limits are logical mesh limits, not whole-process RSS limits.

## Environment and reproduction

Linux x86-64, AMD EPYC 9V74 host, 8 CPU-equivalent quota and 8 GiB memory limit.
GCC 16.2.0, C++26, CMake 4.4.3, Ninja 1.13.2, optional OCCT 8.0.1 (novtk).
Release flags: `-O3 -DNDEBUG -std=c++26 -fno-fast-math -ffp-contract=off`.
IPO and host-specific architecture tuning were off. Shared-host timings are
observations, not performance guarantees. Raw build paths retain the working
folder's rc.1 name; build IDs, source changes and shipped release metadata are rc.2.

```sh
cmake --preset gcc16-release -DPAACS_WITH_OCCT=ON \
  -DOpenCASCADE_DIR=/path/to/lib/cmake/opencascade
cmake --build --preset gcc16-release --parallel 4
ctest --preset gcc16-release --output-on-failure
python3 benchmarks/run_geometry_corpus.py \
  build/gcc16-release/tests/paacs-geometry-corpus --output geometry-results
```

See [geometry guide](GEOMETRY_GUIDE.md), [build guide](BUILDING.md),
[actual mesh gallery](mesh_gallery.png), and `evidence/rc2/` for per-case
JSON/CSV, complete test logs, installed-consumer checks, CLI scenarios,
compiler commands and the machine-readable environment. All original
production source files are retained; new production translation units bring
the total to 68. This package contains source and evidence, not bundled compiler
or CAD-kernel binaries.
