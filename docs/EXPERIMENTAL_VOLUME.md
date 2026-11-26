# Frozen-boundary volume experiments

This opt-in interface belongs to the unpublished 1.3 development line. Published
PAACS Mesh 1.2.0 remains an immutable reference and does not provide these calls.
An experiment must pin its development source revision and wheel hash; matching
the package version alone does not identify an experimental build. No learned
interior model is included in this interface.

The first Python scope is one imported solid, uniform sizing, conforming Tet4, Standard
recovery and the Normal advancing interior constructor. Assemblies, local face
or box sizing, boundary-first ordering and other unsupported controls refuse.
The normal `CadModel.mesh()` entry point and its routing defaults are unchanged.

```python
import paacs_mesh as pm
from paacs_mesh import experimental

model = pm.read_step("part.step", retain_source=True)
boundary = experimental.prepare_boundary(model, size=0.7, max_seconds=300)
boundary.save("saved-boundary")  # Save before attempting any volume fill.

loaded = experimental.load_boundary("saved-boundary")
classical = loaded.mesh()
points = experimental.initial_points(classical, max_bytes=16 * 1024**2)
supplied = loaded.mesh(points)

classical.write_archive("classical.paacsmesh")
supplied.write_archive("supplied.paacsmesh")
print(supplied.report["experimental_boundary_volume"])
```

Use fresh output directories. Bundles contain the original STEP geometry and
are CAD artifacts subject to that geometry's permissions. They contain fixed
local filenames, SHA-256 identities, the original requested controls, preparation
receipt, effective chord bound, units and CAD incidence fingerprint. Boundary
coordinates, connectivity and classifications use the existing mesh archive.
Process-local model tokens and tessellation ordinals are not portable authority.
This first bundle version requires the same native-extension binary hash when
reloaded. Its portable archive encoding does not imply cross-build replay.

Preparation stops after native CAD surface construction, source error/volume
checks and exact shell certification, before interior point placement or volume
construction. Its successful receipt qualifies that boundary under the recorded
checks; it does not promise a successful volume fill or solver adequacy. Saving
first preserves usable surfaces even when the subsequent volume attempt fails.

Loading snapshots bounded regular files, verifies exact source/archive hashes
and build/kernel identities, reimports the STEP, checks the incidence fingerprint,
and revalidates the untrusted boundary against that CAD. Node/face incidence,
native geometric error, CAD face-use winding, volume and embedded-shell checks
remain necessary even when every saved hash is correct. Loading refuses a
mismatch rather than retessellating, remapping tags or repairing the geometry.

`loaded.mesh(points=None)` uses ordinary classical initial placement.
`loaded.mesh(points=[])` deliberately supplies no interior points. Other supplied
point sequences must contain finite millimetre triples within the admitted
count/storage limits. Their order and exact IEEE64 input hash are recorded.
Duplicate or inadmissible points refuse the attempt. Points are not clipped,
projected, silently dropped, or replaced by an automatic cloud after failure.
Every supplied point must remain fixed and be referenced by a published Tet4.
Recovery may add other points, which must be included when comparing mesh size.

All methods use the same downstream construction and native publication gates.
A published result must preserve the tagged, oriented boundary exactly and pass
the usual CAD checks, linear-volume audit and quality checks. An outside node,
changed boundary or missing supplied point cannot become a valid publication.
Failure retains diagnostics without publishing a partial valid volume mesh.

```python
try:
    result = loaded.mesh(points)
except pm.MeshingError as error:
    print(error.report)
    # Available only if the initial cloud was captured before the refusal.
    # An unavailable cloud is distinct from an intentionally empty cloud.
    if getattr(error, "initial_points_bytes", None) is not None:
        try:
            failed_points = experimental.initial_points(error, max_bytes=16 * 1024**2)
        except ValueError as copy_limit:
            print("Initial cloud copy refused:", copy_limit)
```

Initial clouds on failures use a separately bounded immutable binary attachment,
not an unbounded JSON list or a published failed mesh. Requested input and actual
staged cloud hashes are distinct. Tuple copies require an explicit decoded-memory
allowance. A refusal before cloud capture reports that the cloud is unavailable.

Reports retain named source admission, initial placement, constraint recovery,
material fill/postprocessing and validation spans, including failed work where
measured. These are not an exhaustive additive partition: placement includes its
auxiliary numerical topology, and other setup/finalization remains in total time.
Global audit/quality clocks remain in the native report. The Python replay total
also includes point packing and post-call hashes. Preparation, bundle I/O/reload,
sampling, inference, exports and every failed attempt must be charged separately
when reporting a complete CAD-to-mesh workflow.

Resource allowances are explicit and frozen with the prepared boundary.
`max_seconds` is a cooperative deadline for native boundary preparation,
rebinding and volume replay after the operation gate is acquired. It does not
cover fresh STEP import, fingerprint inspection or Python bundle work, and is
not a hard process deadline or a global retry budget. A study must supervise
processes and enforce its aggregate allowance across all those stages and
retries. Modeled memory/output limits do not bound all OCCT/Python allocations
or process RSS.

The 1,024-site limits in the released implementation concern its legacy proposal
cloud and post-recovery interior-front stage. They are not a universal limit on
the main advancing initial constructor. Supplied initial clouds use this
interface's explicit count/output/work/memory admission. Raising a stage limit
or changing the construction policy is a separate treatment, applied equally to
the classical and learned comparisons.
