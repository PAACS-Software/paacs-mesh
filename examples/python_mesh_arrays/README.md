# Inspect a STEP mesh with NumPy

Install a development wheel containing the [NumPy interface](../../docs/PYTHON.md)
and its optional NumPy dependency. From the repository root:

```powershell
python examples/python_mesh_arrays/inspect_mesh.py tests/fixtures/cad/cylinder.step out/cylinder-arrays.npz --size .3 --face-tag 1
```

The example extracts owned arrays from a validated native mesh, computes Tet4
signed volumes and selects boundary triangles by CAD face tag. It prints the
array payload size and body tags, then saves every array to a compressed NumPy
archive. It refuses to overwrite an existing destination.

The archive contains `points`, the node classification arrays, and per-type keys
such as `tet4_connectivity`, `tet4_entity_tags` and `tet4_cell_indices`. The
`metadata` field records units and native validation scope as a JSON string.
Load it using `numpy.load(path, allow_pickle=False)`. It is an array interchange
example; it does not replace the native `.paacsmesh` archive.

Connectivity is zero-based and retains native orientation. Grouped cell blocks
retain their original cell positions and face/body tags. A face tag selects the
model-local CAD identity; it is not a universal semantic label.
