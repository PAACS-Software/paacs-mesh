# Local CAD coverage iteration — 2026-10-03

The fixed 51-file ABC selection now produces **25 accepted whole CAD volume
meshes**, up from 12. All twelve earlier archives remain byte-for-byte identical,
and their geometric receipts agree. The original CAD files are independently
rehashed; no input was replaced or repaired to increase the score. This is a
local development qualification of the optional adapter, separate from the
original 1.1.0 GCC/C++26 release qualification.

The successful meshes pass the adapter's recorded validity and sampled native
geometry checks. Some have extremely poor element quality. Quality optimization
was not the acceptance target for this iteration, and these results should not
be described as uniformly sized or analysis-ready.

## Scope and evidence

The selection is the same 51 files from the downloaded 10,000-file ABC STEP
directory, rather than the complete dataset or an unseen test set. The final
same-binary screen is `evidence/local-2026-10-03/coverage-25/final-complete-screen`.
Its manifest, normal CLI commands, JSON receipts, stdout/stderr and successful
archives are retained locally. The compact dataset/result record is
[cad-coverage-2026-10-03.json](data/cad-coverage-2026-10-03.json).

Each counted success comes from the ordinary production CLI using default
bounded preparation. Private partition experiments, individual successful
pieces and diagnostic/display surfaces are excluded. The final CLI SHA256 is
`161d0dc0a4eb228df1fa58bd341f46b9d22e885bac4d1decad151f9a8abae686`.

| Outcome | Files |
|---|---:|
| Accepted whole mesh | 25 |
| Unsupported contact/body arrangement | 7 |
| Refused with a resource-limit diagnostic | 11 |
| External per-file timeout | 5 |
| Refused invalid/non-solid input | 3 |

The newly accepted source IDs are **4, 5, 6, 8, 9, 17, 18, 19, 20, 29, 1818,
5454 and 7272**. Preserved source IDs are **7, 22, 383, 469, 638, 691, 850,
4149, 5956, 6566, 6807 and 8126**. These are file directory IDs; full paths and
source hashes are in the dataset record. No failing run publishes an archive.

The canonical earlier 12-success screen used **2,000,000 work units and a
60-second timeout**. The new screen uses **20,000,000 work units and a
240-second timeout**, enables explicit separate-contact policies and requests
quality warnings at 0.25. Its size factor remains 0.15 of the CAD bounding-box
diagonal, with one PAACS thread. The improvement therefore measures the new
implementation and the expanded run profile together; it is not an isolated
speed or same-budget comparison.

Other limits are unchanged in this screen: 64 MiB input, 200,000 nodes,
1,000,000 cells, 128 MiB logical output and 128 MiB audit storage. Native
partition attempts use the default 512-piece limit. No per-file face overrides,
hand-picked partition planes or private executables are used.

The eleven `resource-limit` labels do not establish exhaustion of the full
20-million global allowance or OS memory. The wrapper can retain an earlier
local admission refusal while later bounded attempts fail geometry or kernel
search. The mechanism table below preserves those later diagnoses.

## What changed and why it helps

1. **Native surface projection and assembly.** Closest-point validation checks
   all applicable finite native extrema rather than assuming the first result
   is the nearest. Shared native edge identities, seam uses, exact exterior
   embedding and opposite cut-interface triangles continue to determine
   conformity. Native source-state checks are batched without removing the
   authoritative CAD queries or source-immutability check.
2. **More useful copied-body partitions.** Feature planes can use actual native
   cylinder, cone and circle axes. A bounded finer circular policy handles
   curved pieces whose coarser surface representation places a generated point
   outside the original solid. Every partition still must conserve native
   volume and retain OCCT history back to original entities. Fallbacks consume
   cumulative work and record their individual diagnostics.
3. **More reliable star-kernel proposals.** After the existing search fails,
   bounded halfspace proposals and near-parallel numerical retries can find
   an interior kernel. Every proposed point must still pass the exact side
   test for every boundary facet. This corrected a numerical/search gap without
   accepting an uncertified center or changing boundary coordinates. The normal
   CLI now meshes the curved-hole bracket, source 4.
4. **Explicit separate-body contact handling.** Point/edge contacts and a native
   circular rim can preserve separate body nodes. The new opt-in
   `--separate-matched-box-faces` handles complete native axis-aligned box faces
   at eight-corner Tet4 resolution. It proves each box, evaluates endpoints from
   actual CAD vertex identities, and requires unique exact opposite triangles
   covering the complete contact area. This admits the 64-block source 17.

The contact route keeps different body IDs and different node IDs even where
coordinates match. It is a separate contact boundary model. A welded/bonded
multimaterial interface, partial-face contact and general curved shared faces
still need further implementation. Finer box resolutions are refused when
their triangles do not satisfy this supported policy.

No external tetrahedral mesher or retired topology front door was added. The
core target remains independent of OCCT. Experimental general cavity recovery
and local adaptive face/cut proposals are retained as research evidence, but
were not promoted or counted as production coverage.

## Validity, geometry fidelity and quality

Before publication, the native route checks finite coordinates, positive Tet4
orientation, complete face incidence, an exact embedded exterior, original CAD
node classifications/projections, tetrahedron centroid inclusion, sampled
boundary deviation, mesh/native volume agreement and unchanged source state.
Partition interfaces cancel only when keys, coordinates/classes and opposite
winding match. Cross-body contact checks enforce the selected separate-body
policy. Failure, cancellation and exhausted bounds preserve previous results.

This validation level does **not** establish continuous containment of every
point of every tetrahedron in the curved CAD solid, a global tetrahedron-pair
intersection certificate or a Hausdorff distance certificate. The retained
receipts state measured errors and deflection limits. Coarse CAD fidelity must
be evaluated separately from validity of the discretized mesh.

For source 4, the accepted ordinary CLI mesh contains **4,372 Tet4, 3,264
boundary Tri3 and 1,677 nodes**. CAD volume difference is **0.3924%**; maximum
sampled boundary deviation is **0.1084 mm**, below the recorded **0.2828 mm**
sampled limit. Its minimum mean-ratio quality is **4.11e-11**, and **74.9%** of
tetrahedra fall below 0.25. The volume coefficient of variation is **236.9%**.
These warnings are substantial even though the acceptance checks pass.

`--warning-quality .25` retains an accepted mesh and reports low quality.
`--min-quality .25` remains a strict publication gate. Refinement needed to
represent the boundary remains necessary; optional quality optimization cannot
replace geometry checks. General native fallbacks do not yet enforce uniform
volume edge spacing or support arbitrary selected-face/spatial-field controls,
Tet10 or mixed elements.

## Mesh figures and measured distributions

These plots read the published archives. Cyan shows exterior mesh faces; gold
shows intersections through retained cells. The section polygons, element
counts, volume distributions and quality warnings come from actual mesh data.
No replacement/display mesh is substituted. Figure provenance records input,
archive, renderer and image hashes.

The [measured distributions](images/paacs-cad-coverage-25-distributions.png)
and [figure provenance](data/cad-coverage-25-figure-provenance.json) are in the
source package. The [asset provenance](ASSET_PROVENANCE.md) lists the original
ABC inputs for the archived gallery and individual mesh views. See
[the renderer guide](../scripts/mesh_visual_report.md)
for reproducible sections and distribution definitions. Low volume CV alone
does not prove correct geometry or requested sizing.

## Regression and installed-package checks

| Check | Result |
|---|---|
| Fresh complete native Release CTest suite | **64/64 pass**, 818.94 s |
| Final freshly installed consumer examples | **15/15 pass**, 3.00 s |
| Installed component modes | Core, optional/required OCCT available, optional unavailable, required unavailable refusal all pass |
| Fresh OCCT-disabled configure/build/install/components | Pass |
| Original successful meshes | **12/12 identical archive SHA256 and geometric receipts** |
| Source inputs and failed publication | All 51 source hashes preserved; no failure publishes an archive |

The full suite ran after the final native build. Installed core/native libraries
and CLI hashes match that build. New regression coverage includes bounded kernel
proposals, near-parallel constraints, exact-work boundaries, cancellation,
allocation admissions, circular feature axes and explicit box contact rollback
and repeated output. A Windows UI test loader issue was corrected by placing
the test executable beside its actual bridge DLL; UI acceptance/timeouts were
not weakened.

The tests and compact evidence are separate from the historical 58/62-test and
12-success reports. No sanitizer run or new C++26 qualification is claimed.

The final evidence checker additionally passes 40 synthetic accounting and
publication checks; those are separate from the 64 native tests. It binds the
full CTest command log to current test artifacts, the installed libraries/CLI
to the exact Release build, and every profile control to the recorded commands.
The five killed timeouts have no final JSON, so their process publication flags
remain unknown; their actual output archives are confirmed absent. Missing
flags never admit a success. The complete local qualification is retained in
`evidence/local-2026-10-03/coverage-25/final-qualification`.

## Local toolchain and reproduction

The machine is Windows 11 Home, Ryzen 9 8945HS (8 cores/16 logical processors),
32 GiB RAM. The build uses MSVC 19.44.35207.1 / VS 17.14, MSBuild 17.14.8,
Windows SDK 10.0.26100, CMake 4.0.2, Boost 1.89 headers, OCCT 8.0.1 x64 and
TBB 2021.13. The downloaded SDK contains headers, import libraries and runtime
binaries; these remain outside the clean source ZIP.

This is the explicit **C++23 compatibility build** using `/std:c++latest`,
`/MP2 /EHsc`, strict floating point and IPO off. CMake with this MSVC does not
accept `cxx_std_26`; the existing reference C++26 requirements remain separate.
The complete configure command is recorded in the earlier local validation
evidence. With SDK CMake directories and runtime PATH set:

```powershell
cmake -S . -B build/windows-msvc-release -G 'Visual Studio 17 2022' -A x64 `
  -DPAACS_CXX_STANDARD=23 -DPAACS_WITH_OCCT=ON `
  -DPAACS_BUILD_TESTS=ON -DPAACS_BUILD_EXTENDED_PROBES=ON `
  -DPAACS_BUILD_CLI=ON -DPAACS_BUILD_UI_BRIDGE=ON -DPAACS_BUILD_BENCHMARKS=ON `
  '-DCMAKE_CXX_FLAGS=/MP2 /EHsc' -DPAACS_ENABLE_IPO=OFF `
  '-DCMAKE_VS_GLOBALS=VcpkgEnabled=false;VcpkgApplocalDeps=false' `
  -DOpenCASCADE_DIR='<SDK>/cmake' -DTBB_DIR='<TBB>/lib/cmake/TBB'
cmake --build build/windows-msvc-release --config Release --parallel 2
ctest --test-dir build/windows-msvc-release -C Release --parallel 2 --output-on-failure
cmake --install build/windows-msvc-release --config Release --prefix out/native-install
```

The corpus runner is `scripts/run_cad_corpus.py`. The compact dataset record
contains the exact selected paths and runner command. Its normal per-file
flags include `--native-cad --planar-recovery --general-recovery
--separate-contact-bodies --separate-matched-box-faces --warning-quality .25`;
source alias checks and transactional archive publication remain enabled.
The frozen executable and all per-case commands are retained in local evidence.

PAACS explicit visits/admissions are bounded; OCCT parsing, Boolean/projection
internals and allocations are outside those counters. Cancellation is
cooperative. The runner kills its direct child at the external wall timeout
and measures Windows peak working set; it does not enforce a process RSS cap.

## Performance and source footprint

Source 17 was run three times through the final normal CLI: **3/3 pass with an
identical archive**, 384 Tet4 / 512 nodes / 64 bodies / 112 complete face-contact
pairs / 224 matched triangles. Median wall time is **17.63 s** (17.53–18.02 s),
maximum measured peak working set **37.86 MiB**. Stage medians are import 0.371 s,
preparation 11.832 s, generation 0.010 s, native validation 5.291 s, final audit
0.000383 s and quality 0.000078 s. Preparation/recertification dominates this
case. Concurrent regression/prototype activity was present: this is a repeated
functional timing measurement, not an isolated scaling benchmark. Other corpus
timings are single-run coverage measurements. Generation combines surface and
volume work rather than measuring those engines separately.

Comparable physical C/C++ source inventories exclude build, SDK, output,
evidence, docs and vendor files:

| Scope | Previous checkout | Current source |
|---|---:|---:|
| Library `include` + `src` | 167 files / 4,318,590 bytes / 96,188 lines | 216 files / 4,477,375 bytes / 98,009 lines |
| Library + apps/tests/benchmarks/examples | 240 files / 6,838,729 bytes / 152,533 lines | 367 files / 7,481,938 bytes / 163,431 lines |

The library is **3.68% larger in bytes**, and all measured directories are
**9.41% larger**. The optional CAD implementation and its regression coverage
add capability, so physical source size has not decreased relative to that old
checkout. The clean source ZIP excludes local research/evidence, SDKs, runtime
binaries, build trees and mesh outputs; a separate review ZIP retains text
evidence when wanted. The original source checksum manifest and earlier ZIP
are preserved, and each new ZIP has its own verified manifest.

## Remaining gaps and next work

The per-case record and retained attempt diagnostics distinguish failed native
surface conformity, absent strict star kernels, rejected native partition volume
conservation, unsupported contact arrangements, exhausted allowances and timeouts.
A timeout alone does not establish a specific geometry defect. Three refused
inputs are treated separately rather than silently repaired. Their native
validity/topology evidence identifies an invalid solid and two inputs containing
free faces outside the owned solids.

| Most useful observed failure guard | Files | Source IDs |
|---|---:|---|
| Other unresolved original-body interfaces | 3 | 0, 1, 28 |
| Timeout without a final phase/mechanism receipt | 5 | 2, 3, 11, 12, 21 |
| Proven original-body material overlap | 1 | 10 |
| Derived native split fails volume conservation | 4 | 13, 26, 27, 9090 |
| Mixed curved contacts; exact interface type unqualified | 3 | 14, 15, 16 |
| Surface stitching/derived split guards | 1 | 23 |
| Bounded split-star search finds no certified kernel | 6 | 24, 2727, 4545, 6363, 8181, 9999 |
| Original BRep fails native validity | 1 | 25 |
| Raw-valid BRep includes faces outside solids | 2 | 909, 3636 |

These groups choose a useful recorded guard per case, not a proof of exclusive
cause. Source 13 also fails the original surface nonmanifold check. Source 25
has a self-intersecting wire/unorientable shape; 909 and 3636 have three and
four unowned free faces respectively. They are import-contract refusals even
though their raw BReps pass OCCT validity. Source 10 has independently measured
positive body overlap, so a separate contact policy cannot accept it unchanged.
A missing strict kernel is a bounded search result, not proof that no valid
general tetrahedral mesh exists. The full
[per-case failure table](LOCAL_CAD_COVERAGE_FAILURES_2026-10-03.md) retains
uncertainties and links the recorded attempt evidence.

The next work should expand shared-face handling beyond complete boxes and
provide a general constrained Tet4 recovery route for pieces where the bounded
star search fails. Reproduced private cavity experiments still hit bounded recovery
or non-cycle/adjacency failures; increasing a bound alone did not make them
publishable. Native split operations that lose source volume must be corrected
before any mesh is counted. The general surface route also needs robust thin
feature/trim handling and cheaper repeated native queries. These are capability
and geometry-validity gaps, distinct from optional quality improvement.

A separate `paacs-studio` repository and a general arbitrary-CAD guarantee were
not part of this completed 25/51 iteration. Existing frontend/job integration
remains tested and available for later GUI work.
