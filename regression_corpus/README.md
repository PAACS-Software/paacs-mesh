# V0.10 focused pre-OCCT regression corpus

This package retains the existing source-tree tests and adds one focused
integration probe. The entries below identify the representative prepared
models exercised during the V0.10 build. They are geometry-independent and
are intended to be rerun after the OCCT adapter begins supplying real
geometry and topology.

| Case | Purpose | Expected result |
|---|---|---|
| Random point cloud, 50k | Exact compact Delaunay regression | 335,606 Tet4; hash `0xe763362f526a11b1` |
| Random point cloud, 200k | Exact compact Delaunay memory/performance regression | 1,346,539 Tet4; hash `0xfcd63f974312b432` |
| Built-in box, Tet4 | End-to-end surface, recovery, classified Tet4 and audit | 26 nodes; 43 Tet4; strict hash `16723556976053534829` |
| Built-in box, Tet10 | Shared midsides, model placement provider and Jacobian certificate | 118 nodes; 43 Tet10; zero invalid Tet10; strict hash `16508864188371921865` |
| Built-in cylinder, Tet4 | Curved prepared geometry, boundary recovery and audit | 182 nodes; 428 Tet4; strict hash `8085090672163299017` |
| Shared canonical box curves | Build each topological edge once and reuse reversed | 12 canonical edge records; deterministic hash checked by installed consumer |
| Prepared Tri3 sweep | Prism6 layer construction | 2 positive Prism6 cells |
| Prism6 to Tet4 | Transition-catalog integration | Positive Tet4 result |
| T3-A bounded refinement | Local transaction/invalidation behavior | One committed insertion; structured incomplete outcome with boundary preserved |
| Prepared surface adapter | Holes/multiple-loop capable adapter and fault atomicity | Probe passed; zero escaped allocation faults |
| Mapped Hex8 | Retained structured-engine regression | 64^3 cells; hash `10940713674713804929` |

Raw logs are in `../evidence/`.
