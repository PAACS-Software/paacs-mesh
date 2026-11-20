# Known limitations

1. **No OCCT/STEP/BREP integration.** V0.10 finishes the provider boundary but
   does not import or heal real CAD.
2. **Prepared-surface generality differs by API.** Multiple loops/holes are
   available through `SurfaceMesherAdapter`; the simple legacy `GFace` model
   still exposes one boundary loop.
3. **No zero-copy compact topology through the complete Tet pipeline.** Initial
   Delaunay uses the compact private arena; qualified constraint/refinement
   stages use their accepted authoritative workspaces.
4. **Refinement is bounded.** It can return `Incomplete` with retained valid
   private evidence when a target cannot be completed under budgets.
5. **Quality improvement is not universal sliver elimination.** The existing
   bounded T4-A operation and smart smoothing may leave poor but positive Tets.
6. **SPD metrics are a foundation.** Current volume refinement consumes them
   conservatively through a scalar adapter; general anisotropic Delaunay is not
   claimed.
7. **Model-derived orientation is suitable for built-in/prepared geometry.**
   The OCCT adapter should provide authoritative oriented face-use data.
8. **Curved layers remain prepared.** Automatic CAD corner/junction detection
   and general local termination are not complete.
9. **Parallel Delaunay remains experimental and slower.** `SerialExact` is the
   V0.10 production default.
10. **Pre-release C ABI symbol prefix.** The opaque C session retains the
    `paacs_v8_` prefix for compatibility; new Tet request types use `v010`.
11. **Qualification is focused, not exhaustive.** No complete MSVC/macOS/TSan
    or full allocation-failure campaign is claimed for the new orchestration
    layer.
