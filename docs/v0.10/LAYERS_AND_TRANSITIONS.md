# V0.10 layers, sweeps and transitions

V0.10 retains the established prepared generation paths:

- `StraightPrismLayerMesher`;
- `CurvedPrismLayerMesher`;
- `PreparedSweepMesher`;
- mapped Hex8;
- Pyramid5/Prism6/Hex8/Tet4 transition catalogue;
- structured polyhedral generation.

The focused V0.10 probe verifies:

```text
prepared Tri3 source
        ↓
2 geometric layers
        ↓
2 positive Prism6 cells
        ↓
transition catalogue
        ↓
positive Tet4 mesh
```

The retained V8 release probe covers prepared mixed Tri3/Quad4 sweep and hybrid
execution. The retained volume-engine probe also passed.

V0.10 does not add another layer framework. Future OCCT work must identify
wall patches, feature corners, junctions and collision surfaces and then feed
the existing prepared layer requests.

Remaining difficult real-CAD cases include local curved-layer termination,
concave multi-wall junctions and full 2:1 hybrid transition catalogues.
