# Primitive elements and sharp-feature quality — PAACS 1.0.0

[Adaptive box Tet4](LOCAL_SIZE_CONTROLS.md) and [planar Prism6/Hex8 layers](PRISM_LAYER_CONTROLS.md)
are separate APIs. Primitive cylinder Tet4/Tet10 defaults to
[uniform tetrahedra](UNIFORM_CYLINDERS.md); `cylinderStrategy=Swept` selects the
older O-grid conversion. Solid torus/elbow Tet4/Tet10 also use balanced sections;
`--circular-strategy swept` selects the older route.

Every explicit primitive supports `HexDominant`, `Tet4`, and `Tet10`: cylinder,
cone, frustum, tube, torus, and elbow, including the existing hollow variants.
These choices apply to `meshPrimitive`; they do not turn the planar STEP/BREP
importer into a general curved-CAD mesher.

| Choice | Volume cells | Boundary cells | Geometry |
|---|---|---|---|
| `HexDominant` (default) | Hex8; cone tip also Pyramid5 | Quad4; cone tip also Tri3 | Piecewise linear facets |
| `Tet4` | Tet4 only | Tri3 | Piecewise linear facets |
| `Tet10` | Tet10 only | Tri6 | Straight-sided quadratic elements, shared midsides |

Tet10 increases the finite-element approximation order. This route deliberately
retains the linear geometric envelope: midside nodes are arithmetic edge
midpoints, with their actual binary64 rounding included in the Jacobian checks.
It does not project these nodes onto the exact cylinder or torus. The separate
existing `Tet10CurvingMesher` provider route remains available; integrating
analytic curving into these primitive choices is future work. Boundary midsides
may lie off the exact CAD surface, within the reported chord-error bound.

## C++ API

```cpp
paacs::algo::PrimitiveShape shape; // cylinder, radius 1, height 2
paacs::algo::PrimitiveMeshingOptions options;
options.elements = paacs::algo::PrimitiveElements::Tet10;
options.targetSize = 0.2;
options.threads = 4;
options.warningQuality = 0.1;
options.qualityPolicy = paacs::algo::PrimitiveQualityPolicy::Warn;
paacs::mesh::Mesh mesh;
auto result = paacs::algo::meshPrimitive(shape, options, mesh);
if (result.published && result.qualityWarning) {
    // Counts cover every volume cell. The first 32 affected cell IDs are
    // available in result.quality.examples; worstCell is reported separately.
}
```

See `examples/primitive_elements.cpp` for a complete installed-SDK example.
The default element family remains unchanged for existing callers. Rebuild
consumers: the options/result structures gained fields; binary ABI compatibility
with earlier release candidates is not claimed.

## Command line

```sh
paacs-primitive-mesh cylinder .2 cylinder-tet4.paacsmesh --elements tet4
paacs-primitive-mesh cylinder .2 cylinder-tet10.paacsmesh --elements tet10
paacs-primitive-mesh cone .15 cone-tet4.paacsmesh --elements tet4
paacs-primitive-mesh cone .15 swept.paacsmesh --elements tet4 --cone-strategy swept
paacs-primitive-mesh cone .15 warn.paacsmesh --elements tet4 --warn-quality .3 --quality-policy warn
paacs-primitive-mesh cone .15 strict.paacsmesh --elements tet4 --warn-quality .3 --quality-policy reject
```

`--min-quality Q` remains the older **hard** floor, independently of the warning
policy. Setting that floor and `--quality-policy warn` does not disable the floor.
All output-affecting checks happen before mesh publication. A failed meshing
request preserves the caller's mesh and the CLI's existing output file. Archive
I/O failure has the existing writer's separate contract.

## Quality diagnostics

Quality is the signed mean ratio of a tetrahedron, normalized to one for an
ideal regular tetrahedron. A hex or pyramid reports its worst checked subtet.
Tet10 reports the quality of its four corner nodes; its full quadratic Jacobian
is a separate validity gate. This scalar alone is not a solver error estimate.

The result contains the complete poor-cell count, the count of poor cells
incident to a tagged rim or apex, a ten-bin quality histogram, the worst cell
ordinal and quality, and the first 32 affected cell ordinals with their quality
and sharp-feature flag. `examplesTruncated` explicitly marks omitted examples.
This bounds diagnostic memory; it does not imply only 32 cells were checked.
Cell ordinals refer to the published archive's cell order. The CLI emits the
same information as JSON and a concise warning on stderr. Save stdout beside
the archive when retaining diagnostics: the archive format has not been extended
with warning metadata. `warningQuality=0` disables quality warnings, not validity
checks. The default warning threshold is 0.1 and may be changed by the caller.

`Warn` publishes a fully validated mesh with diagnostics. `Reject` refuses the
mesh if any cell falls below `warningQuality`. Neither policy permits zero or
negative volume, invalid boundary incidence, nonfinite coordinates, or an
unresolved Jacobian. An intrinsically sharp angle, thin wall, or requested
resolution may make good isotropic elements impractical. A published warning is
explicitly useful in that case; it is not a claim of suitability for every solver.

## Cone strategy

`Automatic` uses the existing O-grid for HexDominant, and graded tetrahedral
refinement for Tet4/Tet10 cones. `Swept` requests the old disk-sweep construction,
followed by conforming tetrahedral conversion when requested. `Graded` is valid
only for a cone with Tet4 or Tet10; incompatible explicit requests are refused.

The graded route jointly chooses a small base fan (8, 12, 16, 24, or 32 sectors)
and a bounded refinement depth. Each tetrahedron splits into eight children.
Boundary midpoint projection preserves the cone, base rim, and single apex;
interior nodes stay inside. The central octahedron uses the best of three local
diagonals without changing any shared boundary face. This leaves fewer angular
sectors next to the tip and progressively more farther away. It avoids the
O-grid's shrinking, highly stretched terminal layer. Final actual cell signs,
volume, chord bounds, and full boundary incidence still have to pass.

Uniform refinement of the macro tetrahedra is deterministic and simple, but can
overshoot a requested size when another refinement level is needed. It is not
local adaptive remeshing, and `targetSize` remains nominal. Very slender cones
can still have poor elements. Resource limits may refuse a fine Tet10 request;
no silent coarsening or partial mesh is published. For a very slender cone,
`Swept` can be cheaper when its quality warnings are acceptable; choosing that
trade-off stays explicit. Hex-dominant cone-tip quality
is unchanged and is now reported through the same warning mechanism.

## Conformity, validity, and cost

* Balanced cylinders and solid circular sweeps build tetrahedra directly from
  triangular sections. They do not allocate or audit an intermediate hex mesh.
  Their section ordering preserves shared diagonals, including the torus seam.
* Hybrid-to-tetrahedral conversion uses a global vertex ordering on each face
  (pulling triangulation): six tets per hex and two per pyramid. Shared quads
  choose identical diagonals, including periodic seams. No center nodes or
  per-face hash table are needed. The emitted Tet4 mesh receives its own audit.
* Tet10 edges are interned using flat, sorted integer records. Every edge has
  exactly one midside node, shared with adjacent Tet10 and Tri6 cells. Actual
  emitted edge bindings are verified before publication. Corner incidence is
  inherited from the already audited, unchanged Tet4 connectivity. The public
  `audit` field describes that corner audit; `quadraticConnectivityVerified`
  describes the additional high-order checks.
* `positiveTet10JacobianFast` bounds the affine Jacobian matrix over the complete
  reference tetrahedron using outward-rounded interval arithmetic. A strictly
  positive determinant interval is sufficient for acceptance. Unresolved cases
  use the existing interval-Bernstein certifier. A sampled Jacobian never accepts
  an element. Actual Tet10 volume uses degree-three tetrahedral quadrature, so
  rounded midpoint positions are included in the final volume-error gate.
* Generation, quality, and certification use disjoint fixed chunks. IDs,
  connectivity, reports, reductions, and archive hashes are independent of the
  requested thread count. Edge sorting and exact face incidence remain serial.
* Memory is preflighted before each allocation phase, using checked arithmetic.
  Tet10's unique-edge count is resolved before allocating its final storage.
  The auxiliary working budget includes retained source meshes and edge tables.
  A candidate output payload is charged to the output budget; stack, thread-pool,
  allocator overhead and the caller's pre-existing destination are not RSS
  estimates. There is no global all-cell-pair search. Conversion is O(cells),
  edge interning is O(edge uses log edge uses), and records are contiguous.
* Cancellation is checked between chunks and allocation phases. A single in-place
  edge sort and individual allocation are not interruptible; cancellation is
  checked again before publication. All exception and cancellation paths retain
  the strong destination guarantee.

The geometric basis for compatible face triangulations is the pulling
triangulation of a polytope under one vertex ordering; see the authors' treatment
in [Beck and Sanyal, *Combinatorial Reciprocity Theorems*, §5.7](https://matthbeck.github.io/papers/crt.pdf).
This is an established method, not a novelty claim. The implementation additionally
checks emitted cells because swept finite-element cells and floating-point
coordinates require their own geometric validation.
