# Python requests and failure diagnostics

This guide describes additions to the development interface. Use a wheel built
from the corresponding source commit; a development artifact's version number
alone does not identify its Python features.

`CadModel.mesh()` exposes its supported keyword arguments and defaults through
Python inspection and `help()`:

```python
import inspect
import paacs_mesh as pm

print(inspect.signature(pm.CadModel.mesh))
```

The native defaults are preserved: standard recovery, independent parts, one
meshing worker and one validation worker. Sizes and coordinates use millimetres.
Unspecified native budgets retain their native defaults.

## Use a guarded selection directly

Both direct meshing and saved recipes accept `selected_face_sizes` entries as
`(FaceSelection, size)` or `(FaceSelection, size, growth)`. The selection must
belong to the imported model; its body ownership and geometric property guards
are checked before conversion to native face controls.

For an input containing exactly one cylindrical face of radius 1 mm:

```python
model = pm.read_step("cylinder.step")
wall = model.select_faces(surface="cylinder", radius=1, expect=1)
controls = [(wall, .3, 1.2)]

mesh = model.mesh(
    .7, selected_face_sizes=controls, assembly="conforming",
    min_size=.18, growth=1.4,
)

recipe = pm.create_recipe(
    model, .7, selected_face_sizes=controls, assembly="conforming",
    min_size=.18, growth=1.4,
)
pm.save_recipe(model, recipe, "cylinder.recipe.json")
```

The direct request and recipe use the same guarded conversion. Existing
`face_sizes=[pm.FaceSize(tag, size)]` requests continue through native validation.
Controls for the same face are not silently merged. `None` means no optional
controls for `face_sizes`, `selected_face_sizes` and `volume_controls`.

Spatial controls still require one solid and an explicit
`assembly="conforming"` request. Adding controls does not switch the default
`preserve_parts` policy, fuse bodies or relax native validity checks. Selection
expansion is bounded; an oversized request fails before native meshing.

## Interpret the failure that occurred

`pm.diagnose()` accepts structured PAACS exceptions as well as native reports
and saved CLI failure envelopes:

```python
try:
    model.select_faces(surface="plane", expect=1)
except pm.SelectionError as error:
    diagnostic = pm.diagnose(error)
    print(pm.format_diagnostics(diagnostic))
```

A selection failure retains its actual stage, reason, expected count and
candidate information when those fields exist. Candidate face tags remain
selection evidence; they are not labelled as confirmed mesh classifications.
An error before meshing does not imply that volume construction was attempted.
The source hash remains unknown if the supplied error report does not carry it.

For native meshing refusals, the summary's `native_hint` preserves the native
report's existing hint, and the text formatter displays it. No hint is invented
when the native report has none. Recipe validation identifies the affected face
tag or volume-control index where that information is available.

Diagnostics copy their inputs. They do not modify the exception's report, retry
the request, change assembly policy, certify incomplete checks or publish a
partial mesh. A sizing change is still a request to compare and verify, rather
than a guarantee of improved cell quality or simulation accuracy.
