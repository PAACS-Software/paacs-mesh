# Stable 1.2.0 qualification protocol

This document defines the required local STEP qualification and release gates.
It records no new success count or completed stable-release result. The measured
outcomes belong to the matching release assets, with their artifact identities,
options, validation scope and costs.

The qualification uses the installed Python distribution through its public
API. Workers run in isolated processes without source-package injection,
manually configured CAD SDK paths or loader overrides. Local qualification
wheels and final tagged-source publication wheels retain their own source and
binary receipts; development artifacts are not relabelled as published files.
The [installed-STEP runner](../benchmarks/qualify_installed_step.py) freezes a
plan before running the serial attempts.

## Fixed 51-file corpus

The case membership is the exact set in the committed
[RC1 corpus table](step-corpus-1.2.0-rc.1.csv), identified by its relative input
paths and SHA-256 fingerprints. All 51 local source files have been checked
against those fingerprints. Selecting a new first-51 sample from the larger
download directory would not reproduce this set.

For each imported file, let `d` be the diagonal length of the union of the
original bodies' axis-aligned bounding boxes, in millimetres. The requested
global size is `h = 0.15 * d`. This is a requested size, not a certificate of
achieved element spacing.

| Setting | Fixed corpus profile |
|---|---|
| Element family | Tet4 |
| Resource policy | `normal` |
| Assembly policy | `preserve_parts` |
| Recovery | `standard` |
| Construction effort | `legacy` |
| Meshing / CAD validation workers | 1 / 1 |
| Maximum nodes / cells | 200,000 / 1,000,000 |
| Modeled working memory | 512 MiB |
| Mesh output / native audit allowance | 128 MiB / 128 MiB |
| Native cooperative time allowance | 300 seconds |
| Hard worker-process deadline | 600 seconds |
| Independent inspection allowance | 180 seconds |

Normal admission inherits its omitted logical-work policy; explicit count,
memory, output, audit and time limits still apply. The frozen plan records all
effective request and inspection limits, including any process-monitor limits.
Cases run serially so another qualification case does not compete for resources.

The [historical RC1 qualification](RELEASE_READINESS.md) used a different
assembly policy, a 90-second deadline and concurrent processes. Its 32-mesh
result remains historical evidence. The new protocol does not produce a direct
like-for-like coverage or speed comparison with that score. Report changes in
policy and allowances alongside any displayed counts.

## Engineering STEP inputs

Two additional local inputs exercise an IPE beam and a swivel-beam assembly.
Their fingerprints and planned options are fixed before the attempts. The
source CAD is retained locally.

| Input | Global size | Resource policy | Modeled working memory | Native time | Hard process deadline |
|---|---:|---|---:|---:|---:|
| IPE beam | 5 mm | `normal` | 2 GiB | 300 s | 600 s |
| Swivel-beam assembly | 10 mm | `normal` | 4 GiB | 600 s | 900 s |

Both inputs have a default profile. The swivel-beam assembly also has a
separately predeclared advanced profile, giving three planned engineering
attempts. The default preserves parts, uses standard recovery and retains the
documented automatic construction selection. The advanced profile selects
advanced recovery with `native_source_fallback=True` and `boundary_first=True`.
Every profile uses one meshing worker and one CAD validation worker, limits of
500,000 nodes and 3,000,000 cells, 256 MiB output and 512 MiB native audit
allowances, and 180 seconds for independent inspection. The complete request
and process-monitor limits are recorded in the frozen plan.

The advanced result does not erase a failed or slower default attempt. Retain
each planned attempt, its actual options, native report, process outcome and
elapsed time. A result obtained with additional recovery must be identified by
that profile.

## Validity, identity and measurement

No attempt automatically changes geometric tolerances, repairs CAD or removes
features to obtain a success. Requested and effective settings remain distinct
from achieved sizes and measured shape quality. Low-quality valid cells are
reported using their measured distributions; a finer request is not evidence
of improved solver accuracy.

A successful native result must publish a mesh and have complete, valid native
audit and quality reports, with zero invalid connectivity, nonfinite nodes,
inverted cells and reported boundary/orientation defects. Source identity,
original body ownership and expected face/body tags must remain consistent.
Preserving separately valid parts is distinct from establishing a globally
conforming, non-overlapping material domain. Record interference and unresolved
body relations explicitly, including unavailable checks.

The encoded archive is reopened for an independent structural, positive-Tet4
orientation, quality and tag screen. Keep this screen separate from the saved
native CAD validation: reopening an archive does not re-import the CAD, rerun
its geometry audit or certify global element intersections. A VTK output is
hashed as a native export artifact; this qualification does not independently
parse VTK or certify its use by a downstream solver.

Record import, meshing, export, inspection and end-to-end process costs where
available. Keep every failed attempt and internal retry's reported cost; mark
unavailable timing components rather than inventing a decomposition. Report
timeouts, resource stops, invalid input and unsupported construction separately.

Windows process-memory measurements are sampled working-set/RSS observations.
A sampled watchdog is not an operating-system hard memory cap and can miss
short peaks. Native modeled working memory is a separate admission measure.
Native deadlines are cooperative; the recorded hard process deadline bounds a
worker that fails to return. Neither measure is presented as a scaling proof.

## Published evidence and release gates

Measured aggregate results and provenance will be attached to the
[1.2.0 GitHub release](https://github.com/PAACS-Software/paacs-mesh/releases/tag/v1.2.0)
as:

- `step-qualification-summary.json`
- `step-qualification-summary.json.sha256`

These assets contain public aggregate measurements and artifact/protocol
provenance. Private CAD, mesh coordinates, archives, figures, local paths and raw
private-input diagnostics are not included. Detailed local receipts retain the
evidence needed to investigate each attempt. Results are release assets so the
predeclared source protocol stays fixed while measurements are collected.

Local qualification alone does not publish or qualify the release. Before
tagging, the exact final versioned tree must pass all seven native jobs,
including both sanitizer configurations, and all 55 Python jobs. The reviewed,
merged and tagged trees must agree. Publication must then pass all 60 workflow
jobs from the immutable tag, including the fresh 20-wheel/28-consumer matrix,
source receipts, strict package checks, Trusted Publishing and the fresh PyPI
index installation check. See the [publication guide](PYPI_RELEASE.md).

Stable release qualification is complete only after the required local checks
and all release gates have completed, with failures and limitations retained in
the published evidence. This protocol does not claim that those gates have
already passed.
