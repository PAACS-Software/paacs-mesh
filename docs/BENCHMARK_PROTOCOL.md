# Performance methodology

The archive benchmark uses repeated Tet4 connectivity to isolate serialization
and reconstruction; it is not a geometric-quality or meshing-throughput test.
Compare old/new code with the same compiler, optimization, FP flags, allocator,
and hardware. Warm once and report medians plus byte/hash equality. Include
preflight validation in decoding time. Record peak memory separately.

The retained Delaunay benchmark verifies lattice, random, and near-coplanar
goldens before timing seeded point clouds. Its canonical topology hash must
match across runs. Keep raw unconstrained Delaunay timing distinct from the
prepared-boundary pipeline, which includes recovery, classification, and audit.

```sh
build/gcc16-release/paacs-bench-delaunay 50000 7
build/gcc16-release/paacs-bench-delaunay 200000 7
build/gcc16-release/paacs-bench-archive 200000
build/gcc16-release/paacs-bench-pre-occt
```

For comparable meshing measurements, freeze a shared geometry corpus and exact
mesher versions. Include boxes, curved cylinders, holes/cavities, multiregion
interfaces, thin walls, near-coincident surfaces, and CAD fixtures. Measure
success rate, boundary/Hausdorff error, inverted elements, minimum and percentile
quality, node/cell counts, memory, and wall time. Compare equal quality and
geometric tolerance, not just a shared nominal size setting. Publish failed
cases and workload scope. Use the same core count and separate serial/parallel
results. No competing-mesher benchmark was run in this upgrade.

The plausible differentiators to develop are deterministic exact topology,
compact storage, transactional edits, reusable prepared workspaces, and
identity-based hybrid assembly. Novelty and comparative speed remain hypotheses
until measured against a representative corpus and existing literature.
