# Curved primitive meshing — PAACS 1.0.0

Cylinder Tet4/Tet10 defaults to the [uniform route](UNIFORM_CYLINDERS.md).
Balanced solid torus/elbow sections are documented in [circular sweeps](UNIFORM_CIRCULAR_SWEEPS.md).
Explicit `Swept` retains the earlier construction; optional feature analysis is
[separate](SURFACE_FEATURES.md).

All six primitive families offer Tet4 and straight-sided Tet10 as well as
Hex-dominant output. Tet cones select a graded tip strategy by default. Poor
valid elements have bounded diagnostics; strict rejection is optional. See
[ELEMENT_CHOICES.md](ELEMENT_CHOICES.md) for choices and limitations.

## Hex-dominant / swept construction

`paacs::algo::meshPrimitive` is a dependency-free, explicitly parameterized
volume mesher for cylinders, cones, circular frusta, hollow tubes, ring tori
and torus segments (elbows). It emits Hex8 and, at a cone tip, Pyramid5.
The ordinary `Mesh`, archive format, boundary audit, CMake target and writers
are shared with the other engines. The core requires C++26 by default.

```cpp
#include <paacs/algo/PrimitiveMesher.hpp>
paacs::algo::PrimitiveShape shape;
shape.kind = paacs::algo::PrimitiveKind::Torus;
shape.majorRadius = 2;
shape.radius = .6;
paacs::algo::PrimitiveMeshingOptions options;
options.targetSize = .08;
options.maximumChordError = .001; // absolute model units
options.maximumRelativeVolumeError = .005;
options.threads = 4;
paacs::mesh::Mesh mesh;
auto result = paacs::algo::meshPrimitive(shape, options, mesh);
if (!result.published) { /* inspect primitiveReasonName(result.reason) */ }
```

`origin`, `axis` and `xDirection` define a right-handed local frame. Axis
length is normalized; the x hint is projected onto its perpendicular plane.
Near-parallel/zero/nonfinite axes are refused. An axial primitive extends in
the positive axis direction. A torus is centered at the origin, with its main
circle in the local xy plane. An elbow starts on +x and turns towards +y.

| Kind | Active dimensions |
|---|---|
| Cylinder | radius, height |
| Cone | radius, height; one zero-radius top apex |
| Frustum | radius, endRadius > 0, height |
| Tube | radius, 0 < innerRadius < radius, height |
| Torus | majorRadius > radius, radius; innerRadius = 0 or hollow tube radius |
| Elbow | same radii as torus; 0 < sweepAngle < 2*pi |

Unused dimension fields have no effect. Intersecting spindle tori, horn tori
and zero-thickness shells are refused. Reorient the axis and origin to make
a reversed cone. Torus/Elbow `innerRadius` hollows the tube cross-section;
the central donut hole already exists when `innerRadius` is zero.

## Geometry accuracy and validity

A square core and four surrounding patches form an O-grid disk. The outer
perimeter has `4*segmentsPerQuarter` equally spaced points on the true circle.
Hollow shapes use an annular section. There are no duplicate nodes on periodic
seams and no collapsed axis Hex8. A cone has one apex shared by pyramids, not
a plane of coincident hex corners. These construction choices give bounded
local connectivity as the resolution grows, except the intentional apex fan.

Curves are approximated with **linear** facets. `targetSize` is a nominal
resolution, not a guarantee that every diagonal or edge has that length.
Angular counts also satisfy a conservative absolute sagitta bound and a
relative analytic-volume error target. Smaller accuracy targets can refine
the mesh even at an unchanged nominal size. A zero `maximumChordError`
selects 0.0025 times the largest end radius. The default volume error is 1%.
For hollow shapes the combined angular error is additionally limited to a
quarter of one radial layer's thickness. A very thin wall can therefore force
finer angular resolution even with a coarse target size. This fixes a retained
regression where a faceted thin tube had cell centers inside its true CAD bore.
The reported chord bound includes a coordinate-roundoff allowance; it is not
a general Hausdorff certificate for arbitrary CAD.

Validation checks the emitted coordinates, not just ideal construction
formulas. Every hex gets a sufficient Bernstein-coefficient positivity test
for its entire trilinear reference cell. The determinant has degree at most
two in each parameter; 27 evaluations are converted to tensor coefficients.
Strict floating-point arithmetic, normalization and a conservative error
guard are used. An unresolved test refuses publication; there is no exact
arithmetic fallback or recursive subdivision. Positive corner samples alone
are insufficient, and a retained adversarial fixture exercises that case.
This check is also available as `mesh::positiveHex8Jacobian`.

Signed sub-tetrahedron volumes and mean ratios are checked on all volume
cells. Pyramids check both base diagonals and four exact orientation signs.
Optional `minimumQuality` rejects
poor elements. Exact full face-incidence auditing is mandatory. These tests
do not prove the absence of overlap in an arbitrary unrelated input mesh.
The allowed primitive parameters exclude self-intersecting analytic shapes.

**Cone-tip quality decreases with section refinement.** Very thin-wall cells
can also have low isotropic mean-ratio quality because they are anisotropic.
The current fixed
section topology is scaled towards one apex; it does not coarsen locally or
remesh that fan. The results and scaling figure expose this limitation.
Pyramid-cap coarsening or an independently qualified tetrahedral tip region
is needed before promising uniform solver quality under refinement.

## Resource and publication contracts

Exact output counts, connectivity, payload bytes and a conservative sum of
auxiliary arrays are checked before mesh allocation. `maximumCells` includes
boundary cells. Working bytes include section arrays, trigonometric tables,
validation reductions and exact audit scratch; they exclude allocator
overhead, thread-pool storage, caller-owned meshes and archive serialization.
`outputBytes` and `plannedWorkingBytes` are not process RSS. Benchmarks record
RSS separately. The present fast audit preflight restricts node count to
INT32_MAX, in addition to caller budgets; public mesh indices stay 64-bit.

Generation uses exact-size arrays and disjoint writes. Validation folds
fixed-size chunks in a stable order. Shared seams are addressed arithmetically,
not deduplicated with spatial searches. Section preparation, generation and
cell validation take O(nodes + cells) work; the inherited exact audit has its
documented per-bucket sorting cost. Tested hashes match across 1/4 threads.
No new GPU route is enabled, and the existing pool still permits only one
external dispatcher at a time.

Cancellation is checked before allocation, between bounded chunks and
before publication. Allocation itself is not interruptible. A failed call
leaves the destination's logical contents unchanged. Budget, cancellation,
allocation, quality, geometry and numerical failures have explicit reasons.
The output file writer has the existing archive API's non-atomic write
contract; failed meshing requests never start writing it.

Boundary tags are semantic primitive tags: outer surface 1, start cap 2,
end cap 3, inner surface 4; volume 1. Rim edges use 101–104 and the cone apex
uses vertex 201. They are not original STEP/BREP entity identifiers.

## CLI, examples and reproducibility

```sh
build/gcc16-release/paacs-primitive-mesh torus .12 donut.paacsmesh
build/gcc16-release/paacs-primitive-mesh cone .15 cone.paacsmesh --threads 4
build/gcc16-release/paacs-primitive-mesh cylinder .15 cylinder.paacsmesh
build/gcc16-release/paacs-primitive-mesh elbow .15 elbow.paacsmesh --inner-radius .8
python3 benchmarks/run_primitive_scaling.py build/gcc16-release/tests/paacs-primitive-corpus
python3 benchmarks/render_primitive_gallery.py build/gcc16-release/paacs-primitive-mesh gallery.png
```

Rendering needs matplotlib and numpy; generation and measurement do not.
`examples/primitives.cpp` is also compiled through the installed SDK consumer.
The corpus covers thin/slender/shallow shapes, expanding and small-cap
frusta, hollow and nearly horn tori, rotations, translations and scale changes.
Independent tests check analytic polygonal volumes, implicit surface distances,
boundary orientation, seam incidence, Euler characteristic and repeatability.

Optional Open CASCADE tests create matching analytic solids, round-trip them
through STEP/BREP, then independently compare CAD volume, sampled boundary
distances and interior classification. **These use explicit primitive
parameters; they are not curved-CAD import support.** `paacs-cad-mesh` retains
its planar-only import contract. Automatic recognition of imported curved
solids, arbitrary trims, spline faces, blends and industrial assemblies
remain separate work.

Method background: Johnen, Weill and Remacle, *Robust and efficient validation
of the linear hexahedral element* (2017), https://arxiv.org/abs/1706.01613.
The present implementation uses a simple 27-evaluation sufficient test;
it does not implement the paper's optimized 20-evaluation/subdivision method.
