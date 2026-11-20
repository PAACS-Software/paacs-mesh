# Exact boundary audit: design and historical comparison

The published mesh still receives the full `auditLinearVolumeBoundary` check
by default. Faster validation does not depend on trusting a geometry receipt,
suppressing checks, or accepting hash collisions.

## Contract

Supported volume cells: Tet4, Hex8, Prism6 and Pyramid5. Their exterior faces
are Tri3/Quad4. Line2 is ignored as before. Every node must be finite; every
supported non-line row must have in-range connectivity and distinct nodes.
Unsupported types are reported and refused. The audit checks:

- exactly one supplied boundary cell for each singly used volume face;
- opposite cyclic orientation for a face shared by two volumes;
- no supplied exterior boundary on a volume-interior face;
- no duplicate boundary cells, unexpected boundary or nonmanifold face incidence.

As in rc.2, supplied exterior normal sense is the geometry-ownership provider's
responsibility. This is a topological incidence test, not a geometric
self-intersection test, element Jacobian proof, CAD conformity proof, or
high-order audit. Two geometrically overlapping disconnected meshes could
pass this incidence test. Engine geometry checks and ownership checks remain
separate. Positive quality is also a separate requirement.

## Why the grouping is exact

Each face has a unique smallest node because repeated corners are rejected
first. Equal cyclic faces necessarily have the same smallest node, so they
must enter the same bucket. The smallest node is implicit in the bucket;
three remaining corners, face size, owner kind and orientation are stored.
Rotation to the smallest corner and choice of the smaller adjacent corner
canonicalize the cycle. Quad adjacency is retained: sorting just a vertex
set would incorrectly match a crossed quad. Full canonical keys, not hashes,
are compared inside every bucket. Group counts and orientation diagnostics
therefore match the former global-sort algorithm.

Dense minimum-node counters become exclusive starts; scattering advances them
into exclusive ends. This avoids a second cursor array. The face array is
allocated without redundant initialization; every record is assigned before
any sort or read. Nodes below the 32-bit limit use 16-byte temporary records;
otherwise records use 64-bit IDs. The public Mesh remains 64-bit. The actual
large-node-count fallback is not exercised by the local corpus.

Let N be nodes, F face uses and b_i the use count in bucket i. Cost is
O(N + cells + F + sum b_i log b_i); storage is N*sizeof(size_t) plus
F*sizeof(record). Low-valence meshes have small groups. A high-valence hub can
still have an O(F log F) sort; there is no universal linear-time claim. Many
unused nodes also increase the counter array. The declared byte budget can
refuse such a case before allocation. No persistent workspace is retained.

`temporaryBytes` reports requested array storage, not allocator overhead or
process RSS. `maximumFaceUses` and `maximumTemporaryBytes` are checked before
allocating. Failure returns incomplete/invalid; allocation and cancellation
have separate flags. Cancellation is polled during scans and between groups,
not within a single `std::sort`. The caller's mesh is never modified.

## Validation and benchmarking

`tests/audit_reference.cpp` preserves the rc.2 global-sort algorithm; its only
new instrumentation reports temporary array bytes. The differential suite
compares every pre-existing diagnostic field and verifies unchanged input.
Cases cover all supported mixed volume types, shared Hex8/Prism6 interfaces,
node and cell permutations, rotations/reversals, crossed quads, missing and
duplicate faces, nonmanifold uses, invalid rows, nonfinite points, a 12,000-use
hub bucket, exact resource limits and cancellation. Fault injection exercises
both temporary allocations separately.

`run_audit_comparison.py` runs the old and new algorithm in the same compiler
build on identical meshes, alternating order and launching fresh processes.
It includes natural and shuffled numbering at three sizes. Stage timings are
separate from process peak RSS, which also includes mesh generation and any
renumbering. See `../VALIDATION.md` and `../evidence/rc3/audit/` for measured
results. The standalone test executable also supports manual profiling.

## Comparing audit cost

PAACS's old audit dominated a particularly fast affine generator because it
rebuilt and globally sorted all face uses. For performance comparisons, fix
the geometry, boundary fidelity, element type, hardware and thread count.
Record which validity checks each timing includes.
