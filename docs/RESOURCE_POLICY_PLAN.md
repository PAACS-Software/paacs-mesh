# Resource policy implementation and qualification

Status: native implementation, default-policy qualification, paired corpus
checks and visual evidence are complete. Final installed-wheel qualification
is in progress. No release is scheduled by this plan.

The immediate defect is capacity-dependent admission: a larger permitted
auxiliary complex can require a larger bootstrap reservation even before those
cells exist. A valid finer request can consequently start with no interior
points. The implementation must separate future ceilings, current storage,
temporary peak memory and charged computational work.

## Scope

Preserve recovered boundary constraints, body identity, native geometry checks
and element validity gates. Keep historical Bounded construction unchanged;
allow Normal curve integration to adapt to finer spacing. Normal
resource admission must not automatically enable the experimental source/bulk
scheduling policies. Shape quality, sizing completion, geometric validity and
solver adequacy remain separate measurements.

The existing single live numerical complex remains construction machinery;
guard and exterior cells must still pass through the established material,
boundary-recovery and publication checks. It is not a finished body mesh.

## Implementation sequence

1. **Separate bootstrap admission from future capacity.** Build the initial
   candidate once within real temporary-memory and explicit-work allowances.
   Account for its actual cells, mapping work and retained payload while
   preserving peak overlap with the imported complex. Do not allocate or
   charge a full future capacity in advance. Increasing an otherwise inactive
   future ceiling must not make the same bootstrap fail.
2. **Grow the live complex safely.** Between candidate batches, with no
   outstanding insertion plan, compact dead cells when useful and otherwise
   grow storage under a conservative peak-memory allowance. Retry the same
   proposal after admitted growth. Capacity-only growth preserves cell IDs;
   compaction invalidates plans and triggers the existing heap rebuild.
   Account for cavity, relocation, reallocation and compaction scratch, rather
   than replacing the existing peak model with raw cell size alone.
3. **Grow the bulk queue by demand.** Reserve space for actual entries plus
   bounded headroom. Charge actual traversal, comparisons and insertion work.
   A configured queue ceiling is a limit, not completed work.
4. **Expose explicit resource policies.** Add `resource_policy="normal"` and
   `resource_policy="bounded"` across native, Python, CLI and saved recipes.
   Normal with an omitted `max_work` has no practical logical-work ceiling;
   an explicit cap always wins. Memory, node, cell, output, cancellation and
   native validity guards stay active. Bounded retains historical work
   admission. Old recipes missing the policy retain Bounded, independent of
   the eventual default for new calls.
5. **Qualify the default.** Keep the new policy opt-in until evidence supports
   changing new-call defaults. Preserve an explicit Bounded mode and report
   requested and effective resource settings, stage stop reasons, achieved
   sizing and quality distributions. An unfinished construction queue does
   not certify sizing completion; native validity alone does not demonstrate
   simulation convergence.

## Required evidence

| Check | Completion criterion |
| --- | --- |
| Future-capacity admission | A larger unused ceiling cannot block the same bootstrap. Work is charged for actual completed operations; peak memory still covers retained and temporary overlap. |
| Growth and compaction | Stable IDs on capacity-only growth; stale plans rejected after compaction; no repeated full numerical rebuild merely to raise capacity. |
| Resource exhaustion | Small explicit work, memory and output limits stop distinctly. Cancellation remains effective. No invalid or partially published result. |
| Cylinder wall sizing | Global size 0.5 mm with wall sizes 0.3, 0.2, 0.15, 0.1 and 0.08 mm. Report requested/effective/achieved sizes, all three quality metrics, native stop reasons, stage costs and peak process observations. |
| Independent resource effects | Hold stage caps and geometry fixed. Cross 50M, 100M and 500M work caps with 512 MiB and 1 GiB modeled memory; also exercise inherited Normal defaults. Retain every refusal and timeout. |
| Broader and scaling cases | Box, sphere, through-hole and oblique fixtures at multiple global sizes. Retain unchanged meshes where policy is inactive; identify geometry routes unaffected by the new admission. |
| Original CAD regression set | Rerun all 51 cases using fixed settings, one meshing worker and one validation worker. Account for all retries and compare gained, regressed, failed and timed-out cases against the recorded baseline. |
| API and packaging | C++ API, Python signatures, recipes, direct CLI, archive/export tags, optional NumPy, source protection and clean installed-wheel checks pass. Historical recipes preserve allocation meaning. |
| Visual inspection | Fixed-camera cutaways expose interior cells; highlights link to measured archive cell IDs. Include failures and adverse results alongside improvements. |

Native regression tests must use active runtime checks in Release builds.
Timing medians require at least three complete comparable repetitions; both
failed attempts and process time remain visible. Do not select an output only
because its element count increased or its picture looks finer.

## Reproducible experiment plan

The benchmark writes a bounded plan before any worker starts. `--plan-only`
performs no CAD import or meshing. The independent grid keeps all other
settings fixed:

```sh
python -I benchmarks/construction_budget_trials.py --output build/resource-grid --profiles quadruple normal-quadruple --wall-sizes .1 .08 --work-values 50000000 100000000 500000000 --working-byte-values 536870912 1073741824 --repetitions 3 --plan-only
```

Remove `--plan-only` to execute once the rebuilt resource-policy bindings are
available. The grid contains 72 isolated attempts. Every directory name records
policy profile, wall request, work cap, memory cap and repetition. The default
hard deadline remains separate from the native cooperative deadline.

For inherited Normal resource defaults:

```sh
python -I benchmarks/construction_budget_trials.py --output build/resource-normal --profiles normal --wall-sizes .3 .2 .15 .1 .08 --max-work none --max-working-bytes none --repetitions 3
```

`none` means inherit the selected native policy. It is not itself a claim of
unlimited work or memory. Native `resources` records effective settings.
Repeated summaries retain all attempt times. A native mesh-time median is
available only after three or more successful repetitions with the same
source, native binary and archive hash. Raw native reports and refusal reports
remain authoritative; absent counters are reported as unavailable, not zero.

Use the existing `render_construction_budgets.py` for retained mesh cutaways.
Rendering does not rerun native validity checks. Keep source, binary, archive
and saved-report identities paired throughout analysis.

## Default decision

Promote Normal only after the admission defect is reproduced and removed,
resource/refusal tests pass, the baseline is preserved, and refinement/scaling
results support its behavior. Any remaining quality, recovery, topology or
sizing gap must be reported explicitly. Changing policy does not establish
that arbitrary CAD or every requested mesh resolution can be completed.

The completed qualification supports Normal for new Python meshing calls, new
recipes and direct STEP CLI requests. This removes hidden historical admission
ceilings while retaining explicit limits and the 512 MiB modeled memory
default. It does not promise faster execution or better quality for every
request. Historical recipes missing the policy and C++ aggregate options retain
Bounded defaults; construction effort remains Legacy.

Release tags and PyPI publication remain outside this work.

## Measured qualification checkpoints

- The frozen parent and corrected Normal profiles each publish the same 32 of
  51 original inputs within the 90-second process deadline, with one meshing
  worker and one validation worker. Both explicit-work and inherited-work
  Normal profiles preserve that set. This is preservation, not a coverage gain.
- Independent archive screens cover every published mesh through verified
  source, archive and body-tag identities. They supplement the saved native
  CAD checks; they do not repeat CAD classification or global intersections.
- All 36 owned-fixture scaling attempts and 72 fixed work/memory-grid attempts
  publish. Three-repeat summaries require stable identities. The results show
  both improvements and adverse quality/runtime changes.
- The five-size inherited-policy cylinder series completes all 30 attempts,
  with ten stable three-repeat groups. Normal removes the fine-size count
  reversal on this fixture, but very flat individual cells and modeled-memory
  stops remain. The 0.08 mm native median is 111.146 seconds, outside the
  separate original-corpus 90-second criterion.
- Verified whole-element cutaways retain original cell IDs and expose every
  selected poor cell through an x-ray overlay. The
  [resource study](RESOURCE_POLICY_STUDY.md) records measured limitations and
  the distinct final-wheel qualification, rather than inferring that more
  resources complete every requested resolution.

Remaining mesh-quality, fallback-sizing and recovery limitations are separate
algorithm work. They must not be hidden by relaxing validity checks or by
reclassifying an unfinished placement queue as completed sizing.
