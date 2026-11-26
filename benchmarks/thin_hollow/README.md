# Thin hollow STEP qualification fixtures

This standalone generator creates owned CAD for testing the main STEP mesher.
The round outside diameter is 40 mm, the square is 40 by 40 mm, and the
rectangle is 60 by 30 mm. Each has an open through-bore. No meshing, repair or
defeaturing is performed by this tool.

Build against an installed Open CASCADE 8 SDK:

~~~sh
cmake -S benchmarks/thin_hollow -B build/thin-hollow-generator \
  -DOpenCASCADE_DIR=/path/to/occt/lib/cmake/opencascade
cmake --build build/thin-hollow-generator --config Release
build/thin-hollow-generator/paacs-make-thin-hollow \
  build/thin-hollow-fixtures --group full
~~~

For a multi-configuration generator, the executable is in the Release
subdirectory. On Windows, add the SDK bin directory to PATH. On other systems,
make the SDK libraries available through the normal runtime loader. Use the
same SDK resource configuration as PAACS's STEP tests.

Visual Studio projects disable global vcpkg integration and its app-local
dependency copying for this explicitly selected SDK. Use a fresh build
directory when changing SDKs: existing OCCT libraries beside an executable
can take precedence over PATH. The runtime must match the SDK used to build
the generator.

| Group | Fixtures |
|---|---|
| basic | Three unrotated 1 mm wall profiles, length 40 mm. |
| core | Three profiles at lengths 40 and 160 mm, walls 1 and 4 mm: 12 cases. |
| full | Core cases, six rotated/translated thin cases and four rounded rectangular thin cases: 22 cases. |

Rotated cases use a 0.73 radian rotation about (1, 2, 3), followed by a
translation of (11, -7, 5) mm. Rounded rectangles use 4 mm outside and 3 mm
inside corner radii, preserving the 1 mm wall. The manifest records dimensions,
analytic volume, complete placement matrices, topology counts and world bounds.

Construction and STEP reopening check valid BReps, one solid, analytic volume,
millimetre units, matching topology counts and bounds, five void probes
including both bore ends, an interior wall probe and an exterior probe.
Subsequent meshes still require PAACS's native geometry and connectivity checks.
The probes and analytic volume qualify these fixtures, not simulation adequacy.

The output directory must not already exist. A completed manifest is renamed
into place only after all fixtures pass. Failed runs can leave diagnostic
STEP/BRep files and an incomplete manifest staging file; use a fresh directory
for a rerun. STEP timestamps and metadata can change on regeneration, so record
file hashes per trial instead of assuming byte-identical STEP files.

The fixtures support qualification of general initial construction, resource
use, achieved sizing and quality without control volumes or geometry changes.
Future explicit CAD preparation needs separate changed-model evidence. These
valid baselines do not require repair.

## Run through the Python STEP API

The benchmark scripts require Python 3.11 or newer. This is separate from the
Python versions supported by the library. Rendering additionally uses
Matplotlib; it is not a dependency of the meshing engine.

The serial qualification runner records every request, refusal and timeout
separately. It pins the input files, staged Python package, native extension
and supplied SDK libraries. Each attempt uses one meshing worker and one
validation worker. Its output includes the native report and, after success,
the native archive and VTK export.

~~~sh
python benchmarks/thin_hollow/qualify.py \
  --manifest build/thin-hollow-fixtures/manifest.json \
  --stage /path/to/staged-python-package \
  --sdk /path/to/occt \
  --output build/thin-hollow-trial \
  --cases round-t1-l40 square-t1-l40 rectangle-t1-l40 \
  --assembly conforming \
  --profiles default front front-boundary-first \
  --sizes 4 2 --chord 0.1 --memory-gib 2 --seconds 180
~~~

The stage directory contains the paacs_mesh package, not just the extension.
The SDK option is needed for an unbundled development build. Sizes, wall
thicknesses and chord limits use millimetres. The combinations of cases,
profiles and sizes form separate trials; all costs are retained. Add
--plan-only to inspect these requests without importing or running the native
extension.

Explicit front construction requires the conforming assembly policy, even
for one solid. Keep that policy the same across a paired constructor
comparison. The runner defaults to preserve_parts for ordinary import and
meshing trials; it does not silently change assembly policy for an explicit
constructor.

| Profile | Options changed from the ordinary request |
|---|---|
| default | None; automatic constructor selection remains unset. |
| auto-order | Optional source-clearance phase ordering; constructor selection stays automatic. |
| boundary-first | Recover the boundary before initial interior sites on a compatible front route. |
| front | Explicit advancing interior front. |
| front-boundary-first | Explicit front and boundary-first ordering. |
| front-auto-order | Explicit front and optional source-clearance phase ordering. |
| front-source-fallback | Explicit front and the optional source-provider fallback. |
| legacy | Explicitly disable the advancing interior front. |
| advanced | Advanced recovery, boundary-first ordering and source-provider fallback together. |

These names describe requests, not guaranteed constructor outcomes. Use the
returned route, source, sizing and resource diagnostics to establish what
actually ran. An advanced-profile gain does not identify which of its three
options caused the change.

The geometry and generator are original PAACS Software Ltd work under MPL 2.0.
No downloaded CAD assets or modelling histories are used. Generated files and
trial evidence should initially stay under the ignored build directory until
qualification and intended redistribution are reviewed.

## Explicit preparation prototype

The separate paacs-prepare-thin-hollow executable exercises preparation guards
on the owned, unrotated rounded-rectangle-t1-l40 fixture only. It is a benchmark
prototype, not a public PAACS Python repair or defeaturing API. Ordinary STEP
import and meshing do not call it. The operation defaults to none.

~~~sh
build/thin-hollow-generator/paacs-prepare-thin-hollow \
  --input build/thin-hollow-fixtures/rounded-rectangle-t1-l40.step \
  --fixture-id rounded-rectangle-t1-l40 \
  --output build/preparation-unchanged
~~~

The none result supplies an original face catalog with cylindrical radii.
Defeaturing requires an explicit operation value of defeature and a
--remove-faces comma-separated list of original face IDs from that catalog.
It supports only selected cylindrical corner faces in this prototype.
--protect-faces takes another comma-separated selection. Selecting a protected
face for removal refuses immediately. If the operation changes, deletes,
splits or makes another protected face ambiguous, publication is also refused.
There is no nearest-face rebinding.

Repair is a separate explicit operation:

~~~sh
build/thin-hollow-generator/paacs-prepare-thin-hollow \
  --input build/thin-hollow-fixtures/rounded-rectangle-t1-l40.step \
  --output build/preparation-valid-repair \
  --operation repair --max-tolerance 0.000001
~~~

This initial repair probe requires a valid input and checks a valid owned
copy; it does not qualify repair of malformed or invalid raw CAD. The measured
maximum output tolerance must respect the explicit cap, limited here to at
most 0.01 mm. ShapeFix history is not treated as a proven entity mapping:
repair reports history unavailable and refuses protected-face rebinding.
The Perform return value is recorded, not interpreted as proof that an input
needed repair.

Every operation uses a deep copy. Source STEP bytes and serialized original
BRep are compared before publication. Copy lineage is checked explicitly;
defeaturing records provider modified/generated/deleted history and retained
result-face references. Original and prepared integer face IDs are local to
their respective models. Unattributed prepared faces remain explicit.
History references the prepared BRep before export; STEP-reimport face-ID
mapping is not qualified by this prototype.

After valid closed-solid, open-bore, wall-material, tolerance and STEP
roundtrip checks, the prepared STEP and report are published together by
renaming a staging directory to OUTPUT/result. A refusal writes
OUTPUT/report.json and does not create a published result directory.
Incomplete staging files may remain as diagnostics. Output directories must
be fresh.

The prototype bounds input/output to 8 MiB, topology to 128 faces and 1,024
edges/vertices, and checks a cooperative --max-seconds deadline (default 30,
maximum 300) before and after kernel calls. A kernel call can overrun that
deadline; no hard preemption is claimed. Reports include volume, bounds,
maximum topology tolerance and unavailable geometric-distance measurements.
Changed models need separate geometric-deviation and meshing evidence and do
not count as successful meshing of the unchanged baseline.
