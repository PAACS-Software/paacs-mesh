# PAACS Mesh 1.2.0

Version 1.2.0 retains the RC4 meshing algorithms, defaults and Python interfaces.
It provides a versioned baseline for external applications and reproducible
research. Stable package numbering describes the documented release scope;
native STEP meshing remains a preview with explicit unsupported cases and
quality limitations. Packaging maturity is classified as Beta.

The final versioned source must pass the qualification gates below before
publication. The matching GitHub release carries its source and publication
receipts. Earlier local studies retain their own source, platform and settings;
they do not substitute for the stable release checks.

## Python installation and compatibility

Install the exact release with:

```sh
python -m pip install paacs-mesh==1.2.0
```

The source and Python package version are both `1.2.0`. Existing release-candidate
tags and published files remain immutable; stable wheels require their own build
and qualification.

The required matrix remains regular 64-bit CPython 3.12–3.15 on Windows x64 and
Windows 11 ARM64, Linux x64/ARM64 with glibc 2.28 or newer, and macOS ARM64 with
minimum macOS 15. The same macOS wheels must pass installation on macOS 15, 26
and 27. The release requires 20 repaired wheels and 28 independent consumers.
The macOS 27 public-preview runner remains a required check. Intel macOS, PyPy,
32-bit, debug and free-threaded Python builds remain outside this wheel matrix.
CAD runtime dependencies are bundled and NumPy remains optional.

## Resource admission and construction

New Python calls, new saved recipes and direct STEP CLI requests default to
`resource_policy="normal"`. With `max_work` omitted, Normal has no practical
logical-work ceiling; an explicit positive limit remains authoritative. Omitted
modeled working memory retains 512 MiB. Memory, node, cell, output, validation,
cancellation and cooperative time limits still apply. Normal admits actual live
construction storage and grows it on demand, including overlapping temporary
allocations. Modeled memory is not an operating-system hard RSS limit.

Use `resource_policy="bounded"` for historical resource admission and its
50-million omitted-work allowance. Old recipes missing the policy and the C++
`NativeCadMeshOptions{}` aggregate retain Bounded semantics. C++ callers can
explicitly use `normalNativeCadMeshOptions()`.

Construction effort remains `"legacy"`. The separate `"auto"` and `"thorough"`
policies and explicit surface/interior stage caps are opt-in. Reports distinguish
construction limits, achieved sizing, measured quality and native validity.
Increasing a budget does not guarantee completion, better quality or finer cells.
See [construction controls](docs/CONSTRUCTION_BUDGETS.md).

Normal AUTO considers verified primitive, affine-block and linear-sweep methods
for compatible constant-size requests, and resolves eligible independent assembly
parts separately. Guarded fallback retains established successful routes and
charges failed attempts to the same request. Face and volume controls retain
their compatible spatial constructor. The default assembly policy remains
`preserve_parts`; independent part meshes do not certify conforming interfaces
or global non-overlap.

Native changes also preserve checked interior topology through boundary recovery,
regularize curve sizing, retain knot-sided and pcurve correspondence, recover
crossing trim chords and local source facets, and reuse bounded immutable trim
and rejection queries. Bounded local topology edits address numerical slivers;
the final quality and original-CAD checks remain independent. Cancellation is
forwarded into Native Boolean and Section operations, with committed child and
source-recovery work retained in failure receipts.

## Inspection, diagnostics and guarded identity

Python and the CLI add bounded per-element Tet4 quality inspection, native-ID
pagination, spatial distributions, trimmed-CAD centroid distances and
quality-field VTK export. Deterministic SVG views use whole-element cuts and
retain real mesh IDs; comparisons can select a common reported-valid region.
Independent Tet4 shape metrics and guarded mixed-cell volume calculations improve
quality reporting. These measurements do not certify solver accuracy.

Python request signatures and guarded face selections are exposed consistently.
Reports distinguish a bare archive's body tags from source-bound CAD identities,
admit large body summaries within explicit budgets, and preserve source hashes
and query provenance. Capability discovery and readable/JSON diagnostics describe
the available controls, missing measurements, requested policies, actual body
routes, failed attempts and recovery counters. Explicit bounded CAD body-relation
inspection can report contact or overlap for selected pairs; it does not fuse
parts or change assembly policy.

See [mesh inspection](docs/MESH_INSPECTION.md),
[API usability](docs/API_USABILITY.md) and
[diagnostic contracts](docs/DIAGNOSTIC_SCHEMAS.md).

## Hollow parts, physical cap sizing and experimental controls

Eligible Normal single-body native cap sweeps without local controls first try a
physically spaced cap source. A permitted refusal can try the historical cap
source once within the remaining request budget. Cancellation, allocation,
source mutation and internal or kernel failures do not trigger that retry.
`bodies[].native_cap_sizing` records requested size, application, cap counts and
all attempted cap work and time. A successful coarse fallback is reported as
unapplied sizing. Bounded admission retains the historical cap dispatch.

`thin_region_ordering=False` remains the default. The opt-in boolean makes a
bounded source-clearance query before choosing the order of compatible Normal
Tet4 front stages. Explicit `boundary_first=True` takes precedence; selected-face
controls abstain from the unqualified advisory query. The option neither forces
the front constructor nor changes CAD, sizing requests or validity limits.

The separate `boundary_first`, `native_source_fallback`, `trimmed_ray_witnesses`
and advanced recovery controls retain their documented opt-in scope. No importer
or standard meshing call silently repairs or defeatures CAD.

The [owned hollow study](docs/THIN_HOLLOW_QUALIFICATION.md) completed 22 of 22
models locally on Windows with the final experimental setting, preserving all
20 baseline successes and their CAD tags; 14 archives stayed byte-identical.
Its rounded profiles retain poor quality tails. This is a separate owned-model
study, not a rerun or a new score for the original 51 downloaded CAD models.
The [historical 51-file qualification](docs/RELEASE_READINESS.md) records its
32-published-mesh baseline. Stable-release corpus results must identify their own
package, options, validation and timing scope. Paired construction, source-recovery
and resource studies retain their complete attempts and limitations.

A separate benchmark preparation prototype qualifies guarded copying of valid
owned CAD and explicit selected-fillets defeaturing. Invalid-CAD repair and
STEP-reimport face history remain unqualified. Modified models are separate
inputs and do not count as unchanged-source coverage. No new solver or simulation
qualification or arbitrary-CAD completion claim follows from these results.

## Release gates and compatibility

Before tagging, the final versioned tree must pass all seven native qualification
jobs, including GCC address/undefined and thread sanitizers, and all 55 Python
qualification jobs. PR virtual-merge trees must match the merged main/tag tree.
Publication repeats the entire wheel and consumer matrix from one immutable tagged
source and requires resource-policy native/Python parity, exact metadata and
hashes, compiler configuration and dependency receipts, and strict Twine checks.

The GitHub release must carry matching source and qualification evidence before
OIDC Trusted Publishing. The final index job checks all published wheel hashes
and performs a fresh hash-pinned pip install, import, mesh, export and CLI check.
See [compatibility](docs/PYTHON_COMPATIBILITY.md) and
[publication](docs/PYPI_RELEASE.md) for the required process.

The CMake package, SDK and CLI report `1.2.0`.
Public C++ request/report layouts changed: rebuild C++ consumers with matching
headers, libraries and toolchain. Mesh archive format v1 and C bridge ABI version
1 remain unchanged. The C++ core builds independently of OCCT. The project remains
MPL 2.0, copyright PAACS Software Ltd.

The [RC4 notes](https://github.com/PAACS-Software/paacs-mesh/blob/v1.2.0-rc.4/RELEASE_NOTES.md),
[RC3 notes](docs/RELEASE_NOTES_1.2.0-rc.3.md) and [changelog](CHANGELOG.md)
retain the earlier release history.
