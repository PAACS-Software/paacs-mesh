# UI integration in PAACS 1.1.0

The mesher can now run behind a frontend without making that frontend part of
the geometry or meshing implementation. C++26 remains the default. Existing
algorithms, defaults, numerical checks and archive v1 are unchanged.

## Choose the application boundary

| Interface | Coverage | Intended use |
|---|---|---|
| `sdk::JobService` | Five owning adapters, complete native options and reports | Native desktop controllers, application servers, future CAD hosts |
| `MeshTask` extension | Calls to additional native operations; host owns additional non-flat outputs | Prepared Tet, hybrid, feature analysis and custom pipelines |
| Optional `PAACS::ui` shared library | Primitive and adaptive-box adapters, versioned C ABI | Foreign-language bindings without exposing C++ objects |
| `examples/ui/paacs_bridge.py` | Same two routes, JSON Lines protocol and Python interface | Process-based frontends or Python applications |
| `examples/ui/demo_server.py` | Local reference web UI using that Python interface | Experiment with controls and inspect tagged surface geometry |

The C/JSON bridge is **not yet an exhaustive serialization of every native API**.
Native applications retain all C++ functionality. Geometry models, per-body
controls, layered patches and mapped blocks already have owning C++ job adapters;
portable CAD/model-registration endpoints are a later addition. Native detailed
reports remain available; C/JSON exposes reason codes, a diagnostic, summary,
quality warning, counts and mesh pages, not every engine histogram/body report.

## Build and try it

```sh
cmake --preset gcc16-release -DPAACS_BUILD_UI_BRIDGE=ON
cmake --build --preset gcc16-release --parallel 2
ctest --test-dir build/gcc16-release -R 'paacs.sdk' --output-on-failure

# Linux example; use the actual shared-library path from your build.
python3 examples/ui/demo_server.py --library build/gcc16-release/libpaacs_ui.so
```

Open the printed loopback URL. Choose the shape, element type, target size and
strategy, then **Generate mesh**. Expand **All backend controls** for the rest.
Use **Load JSON** with one of `examples/ui/recipes/*.json`; **Save request JSON**
exports the recipe. The field/profile editor accepts the optional top-level
`field`, `profile_controls` and `axial_positions` sections.

For a primitive profile, enable `profile_enabled` and set `profile_length` to
the primitive height. `profile_global_size`, local size and growth specify its
physical-distance envelope. An explicitly disabled profile stays disabled.
The Python protocol enables supplied profile controls by default only when the
caller did not specify `profile_enabled`. The C setter only stores the controls;
enable the profile separately through its named option.

The demo is a replaceable reference, not a committed frontend choice. Python
uses only its standard library. The browser uses HTML, JavaScript and Canvas,
with no bundler or internet assets. Stop the server with Ctrl+C. It binds only
to `127.0.0.1`, checks host/origin and a session token, and never exposes the
stdio bridge's filesystem export command over HTTP. It is not a multiuser or
production network service.

Windows/macOS require building with a suitable compiler and matching C++ runtime;
the shared-library suffix/location changes (`paacs_ui.dll`, `libpaacs_ui.dylib`,
or configuration subdirectories). The new implementation has Linux qualification,
not Windows/macOS qualification. Do not bundle a compiler runtime arbitrarily;
use the chosen toolchain's deployment requirements.

## Native jobs preserve controls and input ownership

```cpp
#include <paacs/sdk/JobService.hpp>
paacs::sdk::JobService jobs({8, 128ULL << 20, 256ULL << 20});
paacs::algo::PrimitiveMeshingOptions options;
options.elements = paacs::algo::PrimitiveElements::Tet10;
options.targetSize = .3;
auto submitted = jobs.submit(paacs::sdk::primitiveTask({}, options));
if (!submitted.job) { /* inspect submitted.reason */ }
// On a GUI timer: submitted.job.snapshot(); submitted.job.cancel();
// After success: auto result = submitted.job.result();
```

`examples/background_jobs.cpp` is a complete installed-package consumer.

| Factory | Ownership and options |
|---|---|
| `primitiveTask(shape, options)` | Copies profiles, axial positions and audit-volume spans; preserves every scalar option |
| `adaptiveBoxTask(request, fieldOwner)` | Retains a shared immutable field; copies audit-volume spans |
| `geometryTask(GeometryJobInput)` | Retains the immutable model/global/body fields; copies per-body profiles; preserves all body controls |
| `layeredPatchTask(sourceOwner, request)` | Retains immutable source patch; copies profile/audit spans; preserves schedule and limits |
| `mappedHexTask(request)` | Copies all three axis fraction spans; keeps classifications, grading and validation mode |

For fields, a raw pointer must be null or match the supplied shared owner;
mismatches throw before submission. Geometry body-field owners are empty or have
one entry per body control. Null body fields inherit the global field. Shared
models, fields and source meshes must remain immutable, including through other
aliases. Factories allocate/copy and can throw; the service's `submit` returns an
admission reason. Raw input spans/profiles must be valid only until their factory
returns. The adapters replace all external cancellation pointers, including
nested audit/profile pointers, with the job's owned atomic flag.

The native extension is a `MeshTask(JobContext&, mesh::Mesh&) -> TaskOutcome`.
Capture ownership, honor `maximumResultBytes()`, propagate the cancellation flag
to the chosen operation and return `published=true` only for validated output.
Additional non-flat results remain application-owned; the service's immutable
result is the ordinary flat mesh plus its supported engine report. Do not wait
for another job, call shutdown, or destroy this service from its own task.

## Threading, cancellation and results

The current core pool supports one external dispatcher. The service serializes
generation across all `JobService` instances linked to the same core instance;
each job can still use its requested internal thread count. Route **all** meshing
through that service while it is active. Raw calls and legacy sessions do not
participate in its gate. Avoid embedding independent static/shared copies of the
engine and assuming they share a scheduler.

States are queued, waiting for the engine, running, succeeded, failed and
cancelled. Poll snapshots on the UI thread. No UI callbacks run from workers.
Stage progress is real; `total == 0` means indeterminate. Current owning adapters
report preparing/generating/finalizing, not detailed per-element percentages.
Custom tasks may report measured work through `JobContext::progress`.

Queued cancellation completes immediately and releases captured input storage.
Cancellation of a running operation is cooperative. A successful `cancel()` is
linearized against publication: the job cannot publish afterward. Cancelling an
already terminal job returns false. Mapped Hex8 has no internal cancellation;
it completes its native call and discards the result when cancelled. Other routes
also have phases, such as a bucket sort or final summary/hash scan, which cannot
be interrupted instantly. Shutdown cancels and joins; initiate it away from the
GUI event thread when a long job may still be running.

Results are `shared_ptr<const JobResult>`. Fetching a result does not copy mesh
arrays, and it survives destruction of its job/service. Failed/cancelled jobs
never expose a partial mesh. Native `report()` includes structured engine refusal
details. Keep the last successful result as the application's displayed mesh
until the new result is accepted. The example also checks a settings revision
before displaying a completed result.

## Memory and large-mesh handling

- `maximumJobs` includes retained terminal job handles. Release handles when done.
- `maximumResultBytes` caps one flat mesh's logical node/cell/connectivity payload.
  The adapters lower the engine output budget accordingly; mapped blocks receive
  an overflow-checked exact payload preflight.
- `maximumRetainedResultBytes` charges published payload until the **last result
  reference** is released, even when the job handle or service is gone. Admission
  of output fails when the remaining budget is exhausted.
- These are logical payload limits, **not a process RSS cap**. Vector capacity,
  allocator overhead, queued model/field snapshots, reports, generation workspace,
  Python/JSON and renderer buffers have additional costs. Keep the engine's own
  work/temporary/audit limits and an application's input/preview budgets.
- C pages have a maximum of 65,536 entries. JSON pages have 4,096. Offset and node
  indices are checked before copying. Coordinates, connectivity, original cell
  ordinals and `(dimension, tag)` are preserved. Empty end pages are valid.
- The reference preview caps 50,000 source nodes, 200,000 scanned cells and 4,000
  surface faces. It reports truncation. It draws Tri6 by corners; it does not
  change Tet10 connectivity or the solver mesh. A production viewport should
  add level of detail and GPU buffers on top of paged/native immutable access.

## Portable controls and protocol v1

The bridge discovers **49 primitive and 19 adaptive-box scalar options** from
the same table used by its setters. This includes shape/frame dimensions, element
and strategy choices, quality policy, thread count, output/working/audit limits,
profile distribution controls and adaptive depth/growth/work limits. Choice
names and representation bounds are included. Bounds alone do not guarantee a
valid geometric combination; the unchanged engine remains authoritative.

Optional arrays expose all 12 `ControlField` node operations, capped profiles
and explicit axial positions. Caller-supplied per-cell audit-volume spans,
custom field subclasses and external cancellation pointers are native-only.
No unsupported key is silently ignored. Field graphs are compiled before a job
is admitted; an invalid graph produces its named reason. Enum/integer/finiteness
errors are rejected transactionally. ABI version/layout checks precede Python
calls; all C functions contain their exceptions.

Run the worker without a server:

```sh
python3 examples/ui/paacs_bridge.py --library /absolute/path/to/libpaacs_ui.so
```

Send one JSON object per line; receive one response per line. The optional `id`
is a string (up to 128 characters), echoed for request correlation.

```json
{"version":1,"id":"controls","method":"schema","params":{"route":"primitive"}}
{"version":1,"id":"start","method":"submit","params":{"spec":{"route":"primitive","options":{"kind":"cylinder","elements":"tet10","target_size":0.4,"maximum_chord_error":0.03}}}}
{"version":1,"id":"poll","method":"status","params":{"job":"1"}}
{"version":1,"id":"mesh","method":"summary","params":{"job":"1"}}
{"version":1,"id":"page","method":"cells","params":{"job":"1","offset":"0","count":1024}}
```

Only ask for a result after `state == "succeeded"`. Other methods are `nodes`,
`cancel`, `release` and `export` (`job`, `path`; local stdio/Python only). Release
after termination; the bridge rejects releasing an active job so cancellation
and completion can be observed. EOF shuts down the service cooperatively.
Success is `{version,id,ok:true,result}`; failure is `{version,id,ok:false,error}`.
Malformed/duplicate-key/nonfinite JSON and lines above 1 MiB are rejected; the
worker drains an oversized line and continues with the next request.

Job IDs, node/cell ordinals, counts and connectivity use decimal strings so
JavaScript does not truncate uint64 values. Hashes are `null` by default. To request the existing strict mesh hash, set
`JobServiceLimits::computeResultHash`, pass `PAACS_UI_COMPUTE_HASH` to
`paacs_ui_service_create_with_flags`, use `Bridge(..., compute_hash=True)`, or
start the JSON worker with `--hash`. Enabled hashes use 16-digit hexadecimal
strings. This optional extra full-mesh pass is separate from mandatory mesh validation. Scalar integer controls are capped at `2^53-1` or their smaller native
range. Treat a job ID as service-local; it is not a persistent project identity.
The summary's ten cell-type counts follow `mesh::CellType` order.

Field nodes carry `operation`, optional `first`/`second` indices, `a`/`b`/`c`
three-vectors and `parameters` with four numbers. Unused dependency indices
default to `UINT32_MAX`. See `ControlField.hpp`, [the field guide](SPATIAL_CONTROL_FIELDS.md)
and `recipes/local-ball.json` for parameter meaning and composition semantics.

## Before expanding OCCT integration

The new application boundary is ready for a CAD adapter. These remain the next
specific additions, not capabilities claimed by this release:

1. An immutable geometry document with a document ID, revision, unit metadata and
   entity catalog. Scope every saved selection to `(document, revision, dimension,
   tag)`; do not assume a tag survives healing, booleans or re-import.
2. Portable model registration and geometry/body/layer endpoints using the same
   controls as the native adapters, plus serialization of detailed body reports.
3. CAD preflight reporting for units, tolerances, open/nonmanifold shells and
   unsupported surfaces, followed by explicit healing choices.
4. Measured stage counters and finer cancellation points in long engine phases.
5. Windows/macOS UI-bridge qualification, packaging of the chosen frontend and
   renderer level of detail; then curved CAD meshing and matched comparisons.

The test record, including inherited qualifications and current limitations, is
in [VALIDATION.md](../VALIDATION.md). This integration work does not establish
arbitrary-CAD support, solver convergence or superiority over another mesher.
