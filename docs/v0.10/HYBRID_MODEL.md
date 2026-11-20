# Hybrid model integration

The V0.10 Tet production output is an ordinary classified `mesh::Mesh`. It can
participate in the retained hybrid executor through:

```cpp
HybridBodyInput body;
body.solidTag = solidTag;
body.sourceKind = HybridBodySourceKind::PreparedFixedMesh;
body.preparedFixedMesh = &tetOutput;
```

Other retained body sources are mapped Hex8, prepared mixed sweep, structured
polyhedral block and a prepared polyhedral mesh.

Shared fixed-cell interfaces use explicit `HybridPortalNodeBinding` records.
The executor verifies local node identity, coordinates and classifications and
publishes one global node. It does not merge by coordinate tolerance alone.

The existing `MeshingPlan` continues to choose among requested/proven engine
capabilities, while the actual generation algorithms remain separate. V0.10's
new Tet orchestrator supplies the previously missing normal way to prepare a
classified Tet4 or Tet10 body before hybrid assembly.

The retained V8 release probe exercises planning and hybrid execution. The
V0.10 focused probe additionally covers a prepared prism layer and its Tet4
transition.
