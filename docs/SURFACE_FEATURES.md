# Optional surface feature analysis (PAACS 1.0.0)

`detectSurfaceFeatures` identifies sharp edges on explicit surface cells. Its
default `enabled=false` returns immediately without allocation or traversing
the mesh. Ordinary primitive meshing never invokes it. It is a read-only
foundation for later feature-driven local sizing and surface patch planning.

```cpp
#include <paacs/algo/SurfaceFeatures.hpp>
paacs::algo::SurfaceFeatureOptions options;
options.enabled = true;
options.angleDegrees = 30;
options.angularToleranceDegrees = 1e-8;
options.minimumLength = .01;
options.includeTagBoundaries = false;
auto features = paacs::algo::detectSurfaceFeatures(mesh, options);
if (features.completed) {
    for (const auto& edge : features.edges) {
        // edge.first/second: mesh node IDs; firstFace/secondFace: cell IDs.
        // flags can combine Sharp, OpenBoundary and TagBoundary.
    }
}
```

The runnable installed-SDK example is `examples/uniform_features.cpp`.

```sh
paacs-primitive-mesh elbow .2 bend.paacsmesh --elements tet10 \
  --detect-features on --feature-angle 30 --feature-min-length .01
```

The CLI reports counts, array bytes, and analysis time in its JSON output.
Feature options alone do not enable analysis; use `--detect-features on`.
`--feature-tag-boundaries on` also selects seams between different face tags.
`--feature-max-edges` bounds the retained list. The library additionally exposes
edge-use and separate working/output-byte limits. CLI analysis uses the mesher's
working-byte limit and caps output at the smaller of 64 MiB and its output limit.
Analysis occurs before writing the mesh archive. A failed analysis request
preserves an existing output file. Feature lists are returned by the C++ API;
the CLI JSON contains counts rather than the complete edge list.

## Exact meaning of controls

- **Angle:** angle between the two oriented adjacent face normals, in degrees.
  An edge is sharp only if its measured angle is greater than `angleDegrees +
  angularToleranceDegrees`. Flat edges have angle zero. At the threshold,
  the edge is excluded. The tolerance is an angular decision band, not CAD
  healing, welding, or a surface-distance tolerance.
- **Minimum length:** omits selected edge segments shorter than this value from
  the returned list. It is not a minimum element size, a whole-feature-chain
  length, or permission to delete geometry. Increasing surface resolution can
  split a long feature into short segments, so use zero to retain every segment.
- **Open boundaries:** included by default when enabled. A surface edge with
  one incident face is an open boundary; it does not have a dihedral angle.
- **Tag seams:** optional, independent of geometric sharpness. Coplanar faces
  with different valid face tags can be selected. Candidate tag counts are
  reported even when selection is disabled.

The angle criterion is an established mesh-processing technique, rather than
a novelty claim; see CGAL's [feature-detection documentation](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__detect__features__grp.html)
and VTK's [feature-edge definitions](https://vtk.org/doc/nightly/release/9.5/html/classvtkFeatureEdges.html).
PAACS adds its own explicit budgets, deterministic output, default-off contract,
and failure handling around its mesh representation.

## Scope, robustness and cost

Tri3, Tri6 and Quad4 surface cells are supported; line and volume cells are
ignored. The routine does not extract a surface from a volume-only mesh.
Tri6 analysis uses corner geometry, and Quad4 uses a normalized Newell normal.
It therefore does not detect curvature variations inside high-order faces.
Coarse tessellation can also create large facet-normal angles on a smooth CAD
surface. These are mesh features, not authoritative recovered CAD edges; choose
the threshold accordingly and retain CAD ownership when it is available.
Normals are computed after translation and scaling to avoid naive area-product
overflow/underflow. Degenerate or numerically unresolved faces are refused.

One flat edge-use array is sorted by node pair and face identity. Adjacent
faces must traverse an edge in opposite directions. Edges with more than two
incident faces, or inconsistent winding, refuse the result. These topology
checks are independent of angle and length filters. Failure returns no partial
edge list; diagnostic counters may describe the completed portion of the pass.
Open surfaces are allowed. This does not certify that a closed surface is
self-intersection-free, that it encloses the requested CAD solid, or that
high-order midsides conform. Keep the existing volume/boundary audits.

For F faces and U edge uses, complexity is O(F + U log U). On the reference
64-bit build, working array payload is `40F + 16U` bytes; retained edges require
56 bytes each. The count pass allocates exact arrays; a second edge scan sizes
the output before allocation. Small sort stack, allocator metadata, the source
mesh, and process runtime overhead are not part of these array-payload counts.
There is no per-edge hash allocation. Cancellation is polled in all linear
passes and before/after sort; the bounded in-place sort is not interruptible.

Tests cover triangulated and quad cubes, scales from 1e-150 to 1e150, planar
tag seams, folded/open surfaces, Tri6 primitive rims, smooth tori, exact thresholds,
resource-boundary refusals, malformed topology, allocation faults and disabled
zero-allocation behavior. The CLI tests compare archive bytes before/after
enabled analysis to prove it does not change density or geometry.

## Next step toward feature-based variation

Do not feed every short tessellation edge directly into a local size field.
First group selected segments into persistent feature chains/patches and link
them to CAD ownership, then derive sizes from physical curvature/proximity and
through-thickness requirements. Integrate with the existing `MeshControls` and
compiled spatial-field/growth components; avoid an independent competing size
contract. Keep that variation opt-in, retain explicit minimum/maximum sizes and
growth limits, and report constraints that geometry or memory budgets cannot
satisfy. No local adaptive remeshing is enabled by this release.
