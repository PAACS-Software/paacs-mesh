# V0.10 Tet10 report

Tet10 remains an elevation of the accepted Tet4 topology.

```text
Tet4
  ↓
unique global source edges
  ↓
one shared midside node per edge
  ↓
geometry-provider placement
  ↓
interval-Bernstein Jacobian certification
  ↓
classified Tet10 output
```

V0.10 adds a default `ModelTet10PlacementProvider` over the abstract geometry
model. Curve edges use the model edge projection capability and surface edges
use face parameterization/evaluation. No OCCT type is present.

The existing T4-B implementation retains:

- deterministic unique-edge identity;
- shared midside nodes;
- curve/surface/interior classification;
- provider stamp/revision validation;
- bounded placement and certificate work;
- conservative whole-element quadratic Jacobian certification;
- atomic publication.

Focused result:

| Model | Output |
|---|---|
| Unit box | 118 nodes, 43 Tet10, zero invalid Tet10, strict hash `16508864188371921865` |

The V0.10 CLI and archive audit now treat a Tet10-only volume as valid when its
Tet10 audit/certificate is valid; it no longer prints the old uninitialized
linear-Tet radius-ratio sentinel.

The existing T4-B algorithm can reject or report incomplete when a geometry
placement cannot be certified under its limits. General real-CAD curving remains
a provider/qualification task for the OCCT phase, not a topology redesign.
