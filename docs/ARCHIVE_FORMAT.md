# PAACS mesh archive version 1

All integers are little-endian; coordinates are IEEE-754 binary64. The layout
is unchanged from v0.10. No ABI-dependent structures are dumped to disk.

| Offset | Size | Field |
|---|---:|---|
| 0 | 8 | Magic bytes: `50 41 41 43 53 4d 38 00` |
| 8 | 4 | Format version, uint32 = 1 |
| 12 | 4 | Reserved uint32 = 0 |
| 16 | 8 | Node count, uint64 |
| 24 | 8 | Cell count, uint64 |
| 32 | 8 | Total connectivity count, uint64 |
| 40 | 32 per node | x/y/z doubles, classification dim/tag int32 |
| After nodes | Variable per cell | Cell record described below |

A cell record is type uint8, three zero reserved bytes, classification dim/tag
int32, followed by signed int64 node indices. Node count is derived from type:
Line2=0 (2), Line3=1 (3), Tri3=2 (3), Tri6=3 (6), Quad4=4 (4), Tet4=5 (4),
Tet10=6 (10), Prism6=7 (6), Hex8=8 (8), Pyramid5=9 (5).

Exact length is `40 + 32*nodes + 12*cells + 8*connectivity`. Extra bytes,
nonzero reserved fields, non-finite coordinates, classifications outside
`dim in [0,3], tag >= 0`, and invalid indices are rejected. All size arithmetic
is checked before allocation. Counts must fit the signed 64-bit Index type
and the host size_t. Structure is validated before reserving decoded storage.

`MeshArchiveLimits` defaults to 1 GiB encoded input, 50 million nodes,
50 million cells, and 500 million connectivity indices. Applications can pass
other budgets to the span deserializer or three-argument file reader. Encoded
size is not a bound on peak memory: decoded cells, input bytes, and destination
mesh coexist until publication. Serializer limits are based on representable
sizes; it can write files larger than the default reader budget.

Every failed serialization/deserialization leaves the caller's destination
container logically unchanged. File writes are not atomic replacement; errors
at close are detected. Archives have no checksum/authentication field; valid
numeric mutations are not necessarily detectable. Use the package's SHA256
manifest for source-integrity checks, and an application-level digest for mesh
files when needed.
