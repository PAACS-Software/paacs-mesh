# Native CAD failure taxonomy — complete normal CLI screen

This retained screen has **25/51 completed whole-source meshes** and **26 nonsuccess cases**. Every original STEP hash was rechecked. Its normal CLI SHA256 is `161d0dc0a4eb228df1fa58bd341f46b9d22e885bac4d1decad151f9a8abae686`. The profile uses 240s, 20M logical work, 200k nodes/1M cells and explicit contact policies; quality 0.25 is a warning, not a required floor.

## Recorded outcomes

| Runner outcome | Cases |
| --- | ---: |
| unsupported | 7 |
| timeout | 5 |
| succeeded | 25 |
| resource-limit | 11 |
| invalid-geometry | 3 |

All nonsuccess attempts leave no published archive. A timeout does not supply a final native diagnosis. A top-level `budget_exceeded` can preserve an earlier local admission refusal while later attempts fail geometry or kernel search; it does not establish global-work or OS-memory exhaustion. Native process peak memory is observational, not a hard process-memory cap.

## Mechanism groups

| Observed group | Cases | Source IDs |
| --- | ---: | --- |
| Other unresolved original-body interfaces | 3 | 0, 1, 28 |
| Timeout; no final native phase/mechanism receipt | 5 | 2, 3, 11, 12, 21 |
| Proven original-body material overlap | 1 | 10 |
| Derived native split fails volume conservation | 4 | 13, 26, 27, 9090 |
| Mixed curved-body contacts; exact interface type unqualified | 3 | 14, 15, 16 |
| Surface stitching and derived split guards | 1 | 23 |
| No strict split-star kernel found by bounded search | 6 | 24, 2727, 4545, 6363, 8181, 9999 |
| Original BRep fails OCCT validity | 1 | 25 |
| Raw-valid original BRep includes faces outside solids | 2 | 909, 3636 |

These groups select the most useful recorded guard for each case; they do not erase other failed attempts. `NoKernel` means the bounded proposal search found no accepted strict kernel, not that a valid general volume mesh is impossible. The JSON retains every native body and partition-attempt receipt.

## Per-case evidence

| Source ID | Runner outcome / native reason | Observed mechanism and remaining uncertainty |
| --- | --- | --- |
| 0 | unsupported / `unresolved_body_interface` | Current interface gate refuses the assembly; independent topology/distance survey has 2 near-zero-distance body pairs. Empty solid Common/zero distance alone does not establish point, edge or face contact. No completed permitted full interface was certified. |
| 1 | unsupported / `unresolved_body_interface` | Current interface gate refuses the assembly; independent topology/distance survey has 1 near-zero-distance body pairs. Empty solid Common/zero distance alone does not establish point, edge or face contact. No completed permitted full interface was certified. |
| 2 | timeout / `no final JSON` | Runner stopped the process after 240s; no final native JSON diagnostic returned. Failure phase and cause are unresolved; earlier private evidence does not identify the phase of this timed-out invocation. |
| 3 | timeout / `no final JSON` | Runner stopped the process after 240s; no final native JSON diagnostic returned. Failure phase and cause are unresolved; earlier private evidence does not identify the phase of this timed-out invocation. |
| 10 | unsupported / `unresolved_body_interface` | Independent native Common proves positive-volume original-body overlap: 1/3=646.235695713 mm³, 2/3=33.3979625265 mm³. Separate contact meshes cannot admit overlapping material without a different explicit material/Boolean model. |
| 11 | timeout / `no final JSON` | Runner stopped the process after 240s; no final native JSON diagnostic returned. Failure phase and cause are unresolved; earlier private evidence does not identify the phase of this timed-out invocation. |
| 12 | timeout / `no final JSON` | Runner stopped the process after 240s; no final native JSON diagnostic returned. Failure phase and cause are unresolved; earlier private evidence does not identify the phase of this timed-out invocation. |
| 13 | resource-limit / `budget_exceeded` | Final cartesian_fifths native split reports volume_mismatch; measured relative native CAD volume change 4.15758652459e-06 exceeds the native split guard 1e-6. Original-body stitched surface diagnostic: nonmanifold. |
| 14 | unsupported / `unresolved_body_interface` | Mixed native assembly: 7 bodies, 214 faces, 131 curved faces; current interface gate refuses it. Neither complete native skin/Common contact dimension nor a conforming permitted interface is established by the retained lightweight topology/distance receipts. The earlier 120s pair-Boolean diagnostic timed out; it cannot establish contact dimension or exclude material overlap. |
| 15 | unsupported / `unresolved_body_interface` | Mixed native assembly: 3 bodies, 42 faces, 25 curved faces; current interface gate refuses it. Neither complete native skin/Common contact dimension nor a conforming permitted interface is established by the retained lightweight topology/distance receipts. |
| 16 | unsupported / `unresolved_body_interface` | Mixed native assembly: 3 bodies, 29 faces, 8 curved faces; current interface gate refuses it. The six-plane box is body 1; its distance to body 3 is 5.423872 mm. Actual touching pair 2/3 involves curved bodies, so the matched full-box-face extension does not address it. Neither complete native skin/Common contact dimension nor a conforming permitted interface is established by the retained lightweight topology/distance receipts. |
| 21 | timeout / `no final JSON` | Runner stopped the process after 240s; no final native JSON diagnostic returned. Failure phase and cause are unresolved; earlier private evidence does not identify the phase of this timed-out invocation. |
| 23 | resource-limit / `budget_exceeded` | Original surface guard: nonconforming_edge; final native split guard: invalid_result. Original-body stitched surface diagnostic: nonconforming_edge. |
| 24 | resource-limit / `budget_exceeded` | Final cartesian_fifths split has 85 pieces; bounded strict-star search refuses piece 1. Original-body star-volume diagnostic: boundary-genus-not-zero. |
| 25 | invalid-geometry / `invalid-shape` | Native import returned invalid-shape. Independent original BRep validity issues: SelfIntersectingWire, UnorientableShape. |
| 26 | resource-limit / `budget_exceeded` | Final cartesian_fifths native split reports volume_mismatch; measured relative native CAD volume change 2.39255843375e-06 exceeds the native split guard 1e-6. Original-body star-volume diagnostic: boundary-genus-not-zero. |
| 27 | resource-limit / `budget_exceeded` | Final cartesian_fifths native split reports volume_mismatch; measured relative native CAD volume change 4.53375871547e-06 exceeds the native split guard 1e-6. Original-body star-volume diagnostic: boundary-genus-not-zero. |
| 28 | unsupported / `unresolved_body_interface` | Current interface gate refuses the assembly; independent topology/distance survey has 2 near-zero-distance body pairs. Empty solid Common/zero distance alone does not establish point, edge or face contact. No completed permitted full interface was certified. |
| 909 | invalid-geometry / `invalid-shape` | Native import returned invalid-shape. Independent raw BRep validity passes, but 3 original faces are outside every solid; solid-only import refuses them. This is a source-representation/import contract limitation, not evidence that the raw BRep fails OCCT validity. |
| 2727 | resource-limit / `budget_exceeded` | Final cartesian_fifths split has 136 pieces; bounded strict-star search refuses piece 35. Original-body star-volume diagnostic: boundary-genus-not-zero. |
| 3636 | invalid-geometry / `invalid-shape` | Native import returned invalid-shape. Independent raw BRep validity passes, but 4 original faces are outside every solid; solid-only import refuses them. This is a source-representation/import contract limitation, not evidence that the raw BRep fails OCCT validity. |
| 4545 | resource-limit / `budget_exceeded` | Final cartesian_fifths split has 150 pieces; bounded strict-star search refuses piece 119. Original-body star-volume diagnostic: no-strict-star-kernel-found. |
| 6363 | resource-limit / `budget_exceeded` | Final cartesian_fifths split has 105 pieces; bounded strict-star search refuses piece 1. Original-body star-volume diagnostic: boundary-genus-not-zero. |
| 8181 | resource-limit / `budget_exceeded` | Final cartesian_fifths split has 116 pieces; bounded strict-star search refuses piece 2. Original-body star-volume diagnostic: boundary-genus-not-zero. |
| 9090 | resource-limit / `budget_exceeded` | Final cartesian_fifths native split reports volume_mismatch; measured relative native CAD volume change 4.37915832124e-05 exceeds the native split guard 1e-6. Original-body star-volume diagnostic: boundary-genus-not-zero. |
| 9999 | resource-limit / `budget_exceeded` | Final cartesian_fifths split has 102 pieces; bounded strict-star search refuses piece 56. Original-body star-volume diagnostic: no-strict-star-kernel-found. |

## Source geometry and research boundaries

- ABC10 original Common volumes prove overlap; contact-only extensions cannot qualify that assembly.
- ABC14/15/16 lightweight topology establishes mixed curved bodies, not definitive interface dimension; the old ABC14 Boolean timeout remains unresolved.
- ABC25 has native SelfIntersectingWire/UnorientableShape. ABC909/3636 raw BReps are valid but have 3/4 free faces outside solids; distinguish representation refusal from invalid BRep.
- Native split volume_mismatch/invalid_result refer to derived Boolean partitions, not a claim that the unchanged source CAD is invalid.
- Retained private ABC24/9999 adaptive whole trials and ABC4545 Grid9 profiles still contain rejected pieces and no whole archive; individual piece/fan gains do not qualify whole-file coverage.
- The private ABC9999 core-phase-work-budget receipt is an alternative constrained-construction admission refusal, not a successful mesh or proof of geometric impossibility.
- The 12 original successful archives remain byte-identical per the separately retained baseline comparison. Quality acceptance uses warning .25, not a strict quality floor.

## Reproduction and retained evidence

Generated from the repository root with `build-failure-taxonomy.py --screen evidence/local-2026-10-03/coverage-25/final-complete-screen --output docs/LOCAL_CAD_COVERAGE_FAILURES_2026-10-03`. Choose a fresh output stem for another screen because the generator refuses to replace existing evidence.

No new CAD execution, source changes, tolerance relaxation or inferred timeout phase was used.

- `final-complete-screen/manifest.json` — SHA256 `2c8d4c37a1589e350ef7da789227e53bf6e82487f072513ca36bb8e076014419`
- `final-complete-screen/results.jsonl` — SHA256 `8252018beef210e064895bad724f79a45973caf2bfad35b97469364751693a0a`
- `taxonomy-complete.jsonl` — SHA256 `2d22746d5a1ae853d57aeeb3decb1e38604eaae944419a213993f3a1dc83d301`
- `separate-box-contact-private/extra-candidate-topology/summary.json` — SHA256 `77ee2774e79f965575a863f6000e4195c9adc2247ff408da26e928ba7edc7959`
- `taxonomy-summary.md` — SHA256 `4121c88610b730ce6311dbaf77a4297aa2ed79ff7e42e6e7b1461e9c712dbd5e`
- `adaptive-24-9999-whole-summary.json` — SHA256 `fea46596413c04b84a546c2ec5e355c6912325d458964f9100dd7e6d173a6c48`
- `private-grid7-9-candidates/00004545-grid9-1024-v2-sweeps256/receipt.json` — SHA256 `8c48ec8adf946468343e4dbdec89a068523e449d9be6873e504795b62462ad79`
- `fixed-core-whole9999-torus2-ymid1-20m/receipt.json` — SHA256 `68eab53f2d594d235994472820540fb4cbe53e178a9f18c125c1f762d0084819`
- `final-baseline-comparison.json` — SHA256 `a331c9f5f72214bd5d0ac544129aa678d6bc2340a76e9675296fe803f08da79e`
