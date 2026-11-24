# Tet4 quality and sizing qualification

This development pass follows the resource-policy study. It does not schedule
a merge, tag or package publication. Preserve native geometry, connectivity,
orientation, material and resource checks throughout.

## Questions to resolve

1. Which thin and curved cases are qualified by the general STEP mesher, as
   opposed to a prepared sweep or geometry-specific example?
2. At which volume-construction stage do the near-flat cylinder cells first
   appear? Hold the surface fixed when isolating volume changes.
3. Why does the sphere fallback return the same mesh for different global
   sizes, and how should incomplete sizing be reported?
4. Can a person or agent distinguish successful geometric validation, shape
   quality, achieved sizing and solver adequacy from one bounded report?

## Implementation and evidence

- Freeze the passing resource-policy build (`c154ab8`) and retain a compact
  thin/curved baseline with source and native fingerprints. Run cases serially
  with one meshing worker and one validation worker. Retain every failure and
  retry. Compare actual counts, size distributions, CAD errors and all three
  Tet4 shape metrics.
- Trace the initial numerical complex, point insertion, constraint recovery,
  material selection and final publication. Determine the first stage with
  severe cells before choosing a correction. Prefer general placement or
  local topology corrections using the existing robust predicates. Do not add
  a fixture-specific control volume or change CAD tolerances to force success.
- Repair sizing dispatch or fallback construction from the observed cause.
  Resource-limited or unattempted construction must not certify sizing
  completion. Requested, effective and achieved measures remain distinct.
- Make quality and incomplete-sizing warnings useful in the default report
  and machine-readable diagnostics. Reuse existing strict quality options;
  retain historical recipe behavior where practical. Poor shape alone must
  not automatically recommend finer elements.
- Test each correction in isolation, then compare the combined build on the
  same fixtures and relevant regression cases. Use actual whole-element
  cutaways, consistent IDs and measured quality values. Qualify the installed
  package before presenting the work as a release candidate.

Thin geometry needs special interpretation: a coarse isotropic request across
a very small thickness can produce valid, highly anisotropic tetrahedra. A
successful volume audit does not establish suitable through-thickness
resolution, acceptable conditioning or stress/flow convergence. Report the
geometry dimensions and achieved thickness resolution rather than treating
every small mean ratio as the same defect.

The earlier 90-second comparison deadline is a benchmark condition, not a
native engine limit. Keep that fixed-deadline score reproducible, but run
larger qualification requests with explicit longer deadlines and retained node,
cell and stage-work measurements. Record all cancelled attempts as well as
completed runs. A timeout alone establishes neither geometry incompatibility
nor linear scaling; test growth in actual output and work across sizes before
making a scaling claim. Never silently override a caller's explicit deadline.
Direct Python `CadModel.mesh()` calls default to `max_seconds=None`; the
bounded recipe workflow has its own adjustable 90-second default. Qualification
scripts must record their explicit native and hard process deadlines.

The first qualification set will be selected from owned fixtures and existing
examples. The earlier plate-with-hole convergence study uses a prepared planar
triangulation and conforming sweep; its evidence must remain separate from
general STEP-mesher results.
