# CAD and figure provenance

## CAD in the source package

| Files | Count | Source |
| --- | ---: | --- |
| `tests/fixtures/cad/` | 33 | [OCCT corpus generator](../tests/occt_corpus.cpp) builds 16 shapes and writes each as STEP and BRep. `malformed.brep` is an intentional parser-failure input. |
| `tests/fixtures/curved-reference/` | 30 | [Primitive CAD reference generator](../tests/primitive_cad_reference.cpp) builds and exports 15 parameterized shapes in both formats. |
| `examples/elasticity_cantilever/cantilever.step` | 1 | [Cantilever generator](../examples/elasticity_cantilever/make_step.cpp) builds the example's rectangular beam. |
| `examples/elasticity_plate_hole/plate-hole.step` | 1 | [Plate generator](../examples/elasticity_plate_hole/make_step.cpp) cuts a cylinder from a rectangular plate. |
| `examples/flow_straight_duct/straight-duct.step` | 1 | [Duct generator](../examples/flow_straight_duct/make_step.cpp) builds the fluid region used by the straight and hybrid duct examples. |

All 66 were reproduced with the project generators. The STEP `DATA` sections
matched after normalizing OCCT's generated product labels. The BRep files
matched after normalizing the OCCT checked-state flag and line endings.

## Figures in the source package

| Figures | Origin |
| --- | --- |
| `docs/images/construction-budget-study.png` (1) | Actual archived Tet4 meshes of the project-owned cylinder fixture. The [study](CONSTRUCTION_BUDGET_STUDY.md), [adjacent figure receipt](images/construction-budget-study-provenance.json) and [renderer](../benchmarks/render_construction_budgets.py) retain source/mesh/report hashes, native cell IDs, common camera and whole-element cut. Saved native reports are hash-paired, without a new native audit. |
| `docs/images/refinement-rc3-*` (5) | Actual published-RC3 meshes of the project-owned cylinder fixture, rendered from exact saved archive coordinates. The [study](REFINEMENT_FINDINGS.md) and [numerical/figure receipt](data/refinement-cylinder-rc3.json) record settings, IDs, hashes and limitations. These are scientific renderings, not generated interpretations of geometry. |
| Cantilever, plate with a hole, straight duct and hybrid duct (7) | PAACS meshes and solver results from the generated example fixtures above. Their READMEs link the source generators and give plotting commands. The README uses the plate example's `achieved-sizing.png`. |
| `examples/python_step_sizing/sizing-sections.png` (1) | Actual native STEP meshes generated through the Python interface from the project-owned cylinder fixture. The [example](../examples/python_step_sizing/README.md) gives commands, control settings, achieved-size measurements and archive hashes. |
| `docs/mesh_gallery*.png` and the RC4–RC8 figures in `docs/` (13) | PAACS-generated primitive meshes, size and layer controls, and benchmark measurements. The [primitive gallery](../benchmarks/render_primitive_gallery.py), [element gallery](../benchmarks/render_element_gallery.py), [circular comparison](../benchmarks/render_circular_comparison.py) and [control plots](../benchmarks/render_local_controls.py) record the underlying geometry and render paths. |
| `native-cad-gallery.png` and `native-cad-sizing.png` (2) | Locally generated OCCT primitives and a disconnected assembly, then meshed by PAACS. The figure hashes match the [gallery receipt](images/native-cad-gallery-provenance.json) and [sizing receipt](images/native-cad-sizing-provenance.json). |
| RC9–RC11 and 1.0.0 controls figures in `docs/images/` (8) | Generated PAACS meshes and sampled size fields. See the [profiled](../benchmarks/render_profiled_controls.py), [complex](../benchmarks/render_complex_controls.py) and [spatial](../benchmarks/render_spatial_controls.py) renderers and the adjacent archive manifests. |
| `paacs-1.1.0-ui-bridge.png` (1) | Scientific render of a PAACS cylinder from the HTTP bridge. The figure records its node, cell and quality measurements. |
| `paacs-cad-coverage-25-distributions.png` (1) | Plot of PAACS cell volume and quality measurements on the ABC benchmark. The plot contains measured distributions, and its case labels refer to the [archived benchmark receipt](data/cad-coverage-25-figure-provenance.json). |

## ABC benchmark material

The [release corpus results](step-corpus-1.2.0-rc.1.csv) record case IDs, input
paths, hashes and mesh results. Obtain the STEP inputs from [ABC Dataset chunk
0000](https://archive.nyu.edu/handle/2451/44309?mode=full). The release corpus
files and the 66 fixtures above use different input sets. Cite the [ABC
dataset and paper](https://deep-geometry.github.io/abc-dataset/)
when using the benchmark results.

The NYU archive describes the dataset as MIT-licensed. The [dataset license
page](https://deep-geometry.github.io/abc-dataset/#license) identifies the
original CAD creators as model copyright holders and refers to [Onshape's
public document terms](https://www.onshape.com/en/legal/terms-of-use). The
following 14 figures that render model geometry are held in private review
evidence. The numerical results and provenance receipts are included in the
source package.

| Archived figure | ABC cases and evidence |
| --- | --- |
| `step-sizing-cad06.png` | `00000006` in the [release corpus](step-corpus-1.2.0-rc.1.csv) |
| `native-cad-abc.png` | `00000007`, `00008126` in the [render receipt](images/native-cad-abc-provenance.json) |
| `native-cad-unsupported.png` | `00000011` alongside generated cases in the [geometry receipt](images/native-cad-unsupported-provenance.json) |
| `native-cad-recovery.png` | `00005454` in the [recovery receipt](images/native-cad-recovery-provenance.json) |
| `paacs-cad-coverage-25-gallery.png` and the nine `paacs-cad-coverage-25-abc*.png` figures | `00000004`, `00000005`, `00000008`, `00000009`, `00000017`, `00000020`, `00000029`, `00005454`, `00007272` in the [coverage figure receipt](data/cad-coverage-25-figure-provenance.json) |
