# Native STEP meshing

Build the optional `PAACS::occt` target with `PAACS_WITH_OCCT=ON`.
See the [automatic STEP build](BUILDING.md#automatic-step-build) and the
[installed CMake package](BUILDING.md). Native STEP meshing is a preview.

## Import and mesh

```sh
paacs-native-cad-preview part.step part-preview.json 0.05
paacs-cad-mesh part.step 2 part.paacsmesh --native-cad --report-json
```

The preview JSON contains display triangles and original CAD face tags.
The mesher uses the STEP geometry and topology directly. STEP coordinates,
sizes and geometric error limits are in millimetres.

The CLI chooses a mesh route based on the solid. It can use a recognized
primitive, a supported planar or swept strategy, or native Tet4 recovery.
Check `archive_published`, `audit_valid` and each body's `route` in the JSON
receipt. A refusal leaves an existing output archive intact.

## Native size controls

A single Tet4 solid can use selected CAD faces and volume boxes together.
For the checked-in cylinder, face 1 is the curved wall.

```sh
paacs-cad-mesh tests/fixtures/cad/cylinder.step .7 cylinder.paacsmesh \
  --native-cad --face-size 1 .38 --size-box -.3 -.3 .7 .3 .3 1.3 .31 \
  --min-size .34 --growth 1.55 --report-json
```

- `--native-sizing` selects the native size-field route for a face-only
  request. A box selects it automatically. Use `--face-size TAG SIZE` and
  `--size-box X0 Y0 Z0 X1 Y1 Z1 SIZE` to add controls.
- `--min-size` clamps the requested field. CAD curvature and boundary recovery
  can require smaller edges. The field uses the smallest active size.
- `--growth G` controls how quickly the requested size increases with distance
  from a selected face or box. Inspect achieved lengths and adjacent Tet4 size
  ratios when the transition matters.
- Recognized cylinder and tube cap profiles remain available through
  `--face-size` on the analytic route. Native size-field controls currently
  apply to one Tet4 solid.

Use [native_step_sizing.py](../scripts/native_step_sizing.py) to save a face
tag and box with the STEP input hash. Its `mesh` command reopens the project,
runs the native CLI and verifies face and body groups in a C3D4 solver include.
The `plot` command uses exact intersections of the stored Tet4 cells.

## Recovery, assemblies and limits

| Option | Effect |
|---|---|
| `--recovery-mode advanced` | Try bounded native edge and face recovery on routes that reach a certified shell. The default is `standard`. |
| `--local-partition-recovery`, `--adaptive-tet-recovery` | Enable optional bounded recovery routes. |
| `--curvature-points-per-turn N` | Enable curvature-based surface sizing on eligible native routes. |
| `--improve-tet-quality` | Try bounded quality improvement on native pieces and after spatial sizing, keeping the prescribed surface fixed. |
| `--independent-parts` | Mesh supported bodies separately and report material interference. The result certifies each part individually. |
| `--max-work`, `--max-working-bytes` | Set work and temporary storage limits. |
| `--warning-quality`, `--min-quality` | Report a weak mesh or require a strict quality floor. |
| `--trace-native-stages` | Record active stages and diagnostics on standard error. |

The optional spatial Tet4 quality pass uses edge swaps and partial interior-node
movements. It preserves CAD boundary coordinates and tags, and retains the native
geometry and validity checks. Its separate work and query limits can leave weak
cells unresolved; inspect the quality distribution rather than assuming the
switch meets a requested quality floor. `volume_quality` appears in the stage
trace when this pass is requested.

The default material checks reject exterior tetrahedron centroids. On eligible
Tet4 recovery routes, `--allow-centroids-within-chord-tolerance` applies a
bounded CAD distance check and additional interior samples. Check the receipt
for chord and volume error, body ownership, and whether all native samples
were checked. Sampled CAD error describes the checked locations.

The [release validation](RELEASE_READINESS.md) gives the matched corpus results.
The [native CAD development notes](NATIVE_CAD_HISTORY.md) describe the recovery
experiments, evidence and detailed option behavior.
