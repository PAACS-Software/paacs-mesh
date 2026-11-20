# Contributing to PAACS Mesh

Start with the [build instructions](docs/BUILDING.md). Use the GitHub issue
tracker for build questions, bug reports and feature requests. Open a pull
request to propose a code or documentation change.

For a bug report or fix, include the smallest input that reproduces the issue
and the command that exercises it. For CAD issues, include the target size,
recovery settings and JSON receipt. Use a STEP file you can share or construct
a fixture in `tests/occt_corpus.cpp`.

## Where things live

| Path | Purpose |
|---|---|
| `include/paacs/` | Installed C++ headers and SDK interfaces |
| `src/algo/`, `src/bnd/` | Mesh generation, boundary recovery and algorithms |
| `src/geometry/` | Geometry adapters, native STEP meshing and CAD validation |
| `src/sizing/`, `src/mesh/` | Size fields, mesh storage, quality and validity audits |
| `apps/` | Command-line entry points |
| `tests/fixtures/`, `tests/` | Reproducible inputs and regression checks |
| `examples/`, `scripts/` | SDK consumers, solver examples and report tools |

## Checking a change

```sh
cmake --preset gcc16-release
cmake --build --preset gcc16-release --parallel 4
ctest --preset gcc16-release --output-on-failure
```

Changes to native STEP should also pass the [STEP fixture check](docs/BUILDING.md#automatic-step-build).
For changes to recovery or sizing, report whether the mesh passed the audit,
its element count and the achieved size. Include end-to-end time and work limits
for performance changes.

## Preparing a pull request

Keep `main` buildable and use a short-lived branch for each fix or feature.
Group commits by purpose, with a description of the problem and resulting
behaviour. Separate unrelated cleanup so each change can be reviewed on its own.

Include the commands and results used to check the change. Preserve existing
validity checks, and put experimental meshing routes behind an explicit option.
Document any supported-input or output-format change before merging. Add only
fixtures and figures whose origin and redistribution terms are known.
