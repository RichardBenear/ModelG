# MG12-H2 higher-arity fixed-substrate development — MG110

We have extended the explicit rewrite engine so tokens can move through relations connecting several sites at once, while preserving the identity of the relation responsible for each proposed move. The initial checks passed: the implemented transport matches an independently constructed operator on small irregular substrates and known two- and three-dimensional grids, with conserved token counts and the agreed one-use-per-tick behavior. This removes the restriction to binary transport links in the reference engine. The chosen linear transport remains equivalent to a weighted pairwise process, so it does not yet establish irreducible higher-order dynamics, emergent spatial dimension or particle behavior. Those distinctions are now explicit rather than hidden in the implementation.

**Status: H2 INITIAL BUILD/DEVELOPMENT PASS; formal H2 statistical qualification remains pending.** H1 was signed off and closed under MG109. Its sources, data and registered scope remain unchanged.

## What changed

ActiveCell, ACell and BCell are typed unordered relations with two or more distinct incident sites. A token can propose a move to any other incident site. For a cell of arity k, the per-destination rate is the species coefficient times cell weight divided by k-1. Each incident cell therefore contributes the same total exit rate regardless of arity. Several incident cells contribute additively; a high-degree node can have a higher total exit rate. This is equal-volume undirected transport, not an automatic degree/volume normalization.

The engine matches the original multi-site incidence, rather than silently replacing the substrate with an ordinary graph. Overlapping cells that offer the same source/destination move retain separate channel keys. Reservations lock participating tokens, not shared read-only cells or sites. Independent tokens can therefore use one cell during a tick; a single token cannot use overlapping proposals twice.

Reaction matching is carried through from frozen H1, but this development batch turned chemistry off to isolate transport. Mixed higher-arity transport with full chemistry is not newly qualified by these tests. The substrate does not change topology.

## Independent operator and endpoint diagnostics

For this particular single-token linear law, the summed rate between two sites can be written as a weighted pairwise generator. We independently construct that matrix from cell incidence blocks, then compare every explicit event/count-state rate to its corresponding matrix row. All **63 reachable fixture states** passed the frozen 1e-12 relative/absolute comparison. Constant fields, row/column conservation and symmetry passed. Pairwise compilation also matched the initial total exit rates of the frozen H1 matcher. This establishes the intended representation equivalence, not a new higher-order force.

| Fixture | Reachable states | Initial total exit rate | Exact endpoint CDF distance | Strict tick endpoint CDF distance |
|---|---:|---:|---:|---:|
| Binary ring | 4 | 2.000 | 0.033977 | 0.018352 |
| Ternary cell | 3 | 1.000 | 0.036045 | 0.027788 |
| Overlapping ternary/quaternary cells | 4 | 1.600 | 0.018803 | 0.018803 |
| Known 2D periodic grid | 9 | 4.000 | 0.064566 | 0.041129 |
| Known 3D periodic grid | 27 | 6.000 | 0.039225 | 0.087777 |
| Distinct A/B resource cells | 16 | 1.800 | 0.030044 | 0.021002 |

Six fixtures × two executors × 256 independent seeded paths gave **3,072 development paths**. Times were 0, 0.05 and 0.1, with strict maximum tick 0.005. All six exact-executor endpoint state-index CDF distances passed the registered development smoke bound **0.117676703**, nominal 99% family coverage under fixed-sampling/reference assumptions. The strict-tick distances are descriptive only. These are small development diagnostics, not the stronger H1 formal joint-law bounds and not a formal H2 accuracy qualification.

The known 2D and 3D periodic grids also reproduced analytic Fourier-mode decay within 1e-11. The mixed-cell test preserved all five overlapping outgoing channels and its total rate 1.6; competing proposals accepted only one move for the single token. Same-cell moves by two distinct tokens were both allowed. Relabeling and explicitly corresponding proposal marks preserved effects. Four malformed-cell cases were rejected.

All paths passed snapshot, reservation, cell-provenance, transport and closed-token checks. Separate rereading of the saved raw counts and eligible event histories reconstructed every observation and **866 accepted events**, and recomputed endpoint distances through sorted/searchsorted arithmetic. It shares fixture/reference data and author, but does not execute the primary dynamics; it is not external replication.

## Dimension calibration

We counted Manhattan balls at fixed radii 4 and 8 on infinite integer-lattice fixtures, independently by coordinate enumeration and the exact combinatorial formula. These are separate calibration structures, not radius-8 balls on the tiny periodic simulation grids.

| Known asymptotic dimension | Sites within radius 4 | Sites within radius 8 | Finite-scale growth slope |
|---|---:|---:|---:|
| 1 | 9 | 17 | 0.917538 |
| 2 | 41 | 145 | 1.822357 |
| 3 | 129 | 833 | 2.690945 |

The slope is log(B8/B4)/log(2). At these finite scales it is below the known asymptotic dimension. This is a calibration result, not an error to eliminate by choosing favorable radii. A future substrate's dimension must be reported with its scale dependence; cell arity and species count cannot substitute for it. No dimension is assigned to the irregular tiny-cell fixtures.

## Execution, retention and limits

User instruction after H1 sign-off authorized the next phase. Contract, fixtures, all seeds, source hashes and the 300-second single-worker development cap were published at **afef059** before random draws or matrix integrations. Development took **9.310173 seconds**, with one numerical-library thread. There was no failure, retry, parameter selection or criterion change. Existing closed kernels remain frozen.

The full count/reference arrays, seeds, fixture definitions, eligible-event traces and diagnostics are retained in development/. The main result is development/results.json; the separate reread is saved_evidence_check.json. See retention_verified.json for post-publication binary and metadata verification. Raw trace/data comparisons distinguish JSON values from Git line endings.

The next step is a prospectively specified formal higher-arity transport qualification, including joint-law bounds, normalization/label controls and a measured sample/cost plan. That should precede particle tests on new substrates. Cooperative multi-token transport, variable cell volumes, topology evolution and shared-substrate resource regeneration remain separate choices. H1 qualification does not automatically cover these new fixtures or the extension's mixed chemistry/transport dynamics. No particle or emergent-space result is claimed here.

Publication verified at c040b49: all retained binaries, metadata and report match published bytes/values after stated text normalization; source hashes unchanged. See retention_verified.json.
