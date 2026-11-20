# V0.10 performance notes

V0.10 is an integration release. Delaunay and structured-engine performance
were treated as regression contracts rather than new headline targets.

## Compact exact Delaunay regression

Clang 17, C++20, `-O3 -DNDEBUG -march=native`, warmed median:

| Cloud | Final Tet4 | Certified time | Throughput | Canonical hash |
|---|---:|---:|---:|---|
| 50k | 335,606 | 0.2781 s | 1.207 M Tet4/s | `0xe763362f526a11b1` |
| 200k | 1,346,539 | 1.2205 s | 1.103 M Tet4/s | `0xfcd63f974312b432` |

The 200k compact memory counters remain:

```text
arena bytes:     61,605,504
scratch bytes:   25,070,584
workspace bytes: 86,676,088
```

The source Delaunay implementation was not changed by V0.10. Timing variation
relative to earlier V0.9.3 evidence reflects the current runner/build, while
the canonical topology hashes remain identical.

## Structured regression

The retained reusable mapped-Hex path produced 262,144 Hex8 cells in a median
0.002254 s on the 64^3 regression, or 116.3 million Hex8/s, with hash
`10940713674713804929`.

## Integrated workflow timings

These are small prepared-model engineering timings, not competitor claims:

| Case | Median |
|---|---:|
| Box quality Tet4 | 0.569 s |
| Box certified Tet10 | 0.579 s |
| Cylinder quality Tet4 | 0.964 s |

The production routes include surface construction, constraint recovery,
classification, audit and publication; they are not raw point-set Delaunay
benchmarks.
