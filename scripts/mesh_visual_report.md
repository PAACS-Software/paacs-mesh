# Actual mesh visual reports

`mesh_visual_report.py` produces cyan/gold cutaways, exact cell sections,
volume and equivalent-size distributions, tetrahedron mean-ratio histograms,
and a provenance JSON with input, script and figure SHA256 hashes. It reads
PAACS `PAACSM8` version 1 archives directly, preserving double coordinates,
or the legacy ASCII `UNSTRUCTURED_GRID` files produced by `io::writeVtk`.
The reader accepts native archives and supported legacy ASCII VTK cell types.

Install the pinned plotting requirements in a local environment, then run:

```powershell
python -m pip install -r scripts/mesh_visual_requirements.txt
python scripts/mesh_visual_report.py --self-test
python scripts/mesh_visual_report.py --manifest out/visual-validation-2026-10-03/abc-manifest.json --output out/visual-validation-2026-10-03/abc
```

The reproduction driver runs the existing native sizing test and installed
primitive CLI, captures their receipts and produces all three initial groups:

```powershell
python scripts/run_mesh_visual_validation.py --archive-verification evidence/local-2026-10-02/cad-native-evidence/archive-verification.json --spatial-executable build/windows-msvc-release/tests/Release/paacs-spatial-controls.exe --primitive-executable out/windows-install/bin/paacs-primitive-mesh.exe --output out/visual-validation-2026-10-03
```

This command requires the five archives referenced by the existing native
verification receipt. Add `--skip-generation` to redraw existing example inputs
without rerunning their generators. A local `pip --target` installation also
works: set `PYTHONPATH` to that directory before running the command. The local
session uses `out/visual-validation-2026-10-03/python-deps` and the bundled Python
runtime; the renderer records the actual Python/numpy/matplotlib versions.

Paths in the manifest resolve relative to the manifest file. A minimal input:

```json
{
  "title": "PAACS / Supported ABC planar meshes",
  "gallery_mode": "cutaways",
  "columns": 3,
  "distributions": false,
  "cases": [
    {
      "name": "ABC 22",
      "mesh": "../../evidence/local-2026-10-02/abc-screening/runs/0022-00-a9624d63/mesh.paacsmesh",
      "units": "mm",
      "source": "../abc_0000_step_v00/00000022/00000022_ad34a3f60c4a4caa99646600_step_002.step",
      "section": {"axis": "y", "fraction": 0.5},
      "verification": {"verified": true, "divergence_volume": 15625, "minimum_quality": 0.7559526299369238}
    }
  ]
}
```

Optional case fields:

- `report: false` creates only a gallery for that case.
- `section.value` specifies an absolute coordinate instead of a bounding-box
  fraction; `section.keep` is `+1` or `-1` for the retained halfspace.
- `section.expected_section_area` supplies an independent analytic area oracle
  (the generated unit-cube sizing gallery uses one).
- `camera: {"elevation": 27, "azimuth": -58}` controls the view.
- `show_midsides: true` marks visible straight Tet10 exterior midside nodes.
- `field_csv` hashes native field samples for traceability. The measured section
  color remains cube root of actual volume, not that requested field.
- `metadata` carries generation controls and the source native receipt.
- `verification` may contain the full native archive verification object. The
  renderer checks its volume/minimum tet quality independently when present.

Set `gallery_mode` to `sections` for a shared scale sizing gallery. Only combine
commensurate units and comparable definitions in a shared scale. `distributions`
adds an overlaid volume, equivalent-size and quality comparison figure.

`geometry_cases` optionally describes **actual** OCCT display tessellations:

```json
{"name":"ABC 11", "geometry":"abc11-tessellation.json", "source":"model.step", "status":"Refused", "reason":"Unsupported curved faces", "metadata":{"linear_deflection_mm":0.1}}
```

Each tessellation is legacy ASCII VTK surface cells or JSON with `points`,
`triangles` and optional `face_tags` arrays. Apply OCCT face locations and face
orientations in the native exporter. The overview explicitly labels these as
geometry approximations, with no claim of a supported PAACS volume mesh.

Sections intersect the original convex cell edges. They do not discard whole
cells by centroid or invent replacement tetrahedra. A cut through Tet4 may be a
quadrilateral; a shared coplanar face is rendered once. Mixed Prism6/Hex8/
Pyramid5 cells must be convex with planar corner faces; the renderer refuses
nonplanar/nonconvex cells rather than claiming an exact cut. Tet10 midpoint
positions are checked and curved Tet10 is refused. Tet quality is the standard
corner mean ratio. Linear mixed quality reproduces the native fixed minimum
`SubTetMeanRatio` decomposition from `MixedMeshQuality.cpp`; it is a diagnostic
whose value depends on that prescribed decomposition, not a family-independent
ideal-shape score. These diagnostic subtets are never used to replace the mesh
in a figure. Native full boundary, orientation, geometry and resource checks
remain separate evidence.

Volume CV is population standard deviation / mean volume. Cube root of volume
is an equivalent size, not target edge size. Uniformity, shape quality and
geometric accuracy measure different things; none replaces the other. The JSON
records volumes, quantiles, minima, section area, native comparisons and scope.

An optional case `quality_warning_target` records the number and fraction of
volume cells below the preferred quality and marks that threshold on the
individual histogram. This fraction counts cells; it is not volume weighted.
The plotting threshold never changes the mesh or its acceptance decision.
