# MG025 independent review of MG7 results

**Recommendation to the PI: accept NOT SCORED (feedback), with CHANGED at b = 27.5 and SURVIVES at b = 28.5.** The registered mass-scaling comparison is unavailable because the low endpoint loses its particle. The recorded feedback outcome is scientifically informative. I recommend correcting the survival-boundary interpretation and completing the descriptive reporting before closing the record. No replacement production runs are needed to adjudicate MG7.

Reviewed on 2026-11-04 at repository revision `0b86bf7`, against registered MG7 Draft 2, frozen at `5b44cde`, and the pre-run specification in GATES.md.

## Replication and outcome

I ran the repository's independent replay in strict mode on Windows with NumPy 2.2.4. All 16 candidate hashes matched. It verified endpoint-file hashes against the manifest, the recorded b values and absence of a steady-state shift, regenerated node sets, reproduced disc and annular counts where a particle was present, and redetected all 336 scheduled checks. It reproduced both feedback outcomes and the gate, exiting successfully.

My additional script independently regenerated all candidate hashes and reconstructed total populations at all 31,984 edited boundaries, including periods after particle loss. It checked acceptance totals, absence of accepted ephemeral candidates, population accounting, reported bounds, quadrature comparison tolerances, run completion, count-based window means and the production ledger. It also used the independent replay detector to recompute the alternative-threshold counts. The 16 streams collapse to the expected four hashes, one per start shared across arms and endpoints.

| First scheduled check with no particle, b = 27.5 | D0_p0 | D0_p1 | D1_p0 | D1_p1 |
|---|---:|---:|---:|---:|
| CP2 | 20 | Single throughout | 20 | 140 |
| SH2 | 160 | 90 | 60 | 110 |

These are first observed zero-count times at ten-unit spacing. They are not exact continuous-time dissolution times. There was no return to a detected particle at later checks. All eight b = 28.5 runs stayed single. CP2 therefore scores CHANGED at 27.5 (1/4 single) and SURVIVES at 28.5 (4/4 single). Seven low-endpoint runs violate the single-throughout-window requirement. The primary gate correctly stops at NOT SCORED (feedback), before any test of the endpoint change.

The alternative MG4c thresholds yield the same single-in-window verdict for every run. In the seven dissolved runs, the absolute final smoothed-X minimum is at most 1.34e-8, compared with negative thresholds near −4. The final loss is robust to this threshold choice. The data do not support recovering a mass-scaling result by averaging surviving runs, replacing missing particle-centred measurements with zero, or changing the window.

The independent Part 0 count-mean differences reproduce the registered checks. Step refinement gives 0.167 nodes at each endpoint, below 2. Tick refinement gives 0.5 and 1.0 nodes, below 6. Depth and count checks pass. The production ledger sums to 73.00584 core-hours for exactly 16 completed runs. Part 0 is separately recorded at about 8.63 core-hours.

## Corrections to the interpretation

### 1. A precise new survival edge has not been located

GATES.md says the particle's lower survival edge lies above 27.5. The supported statement is that **b = 27.5 failed this finite-panel survival criterion through t = 200 under both tested arms**. A single CP2 run survived there. There was no survival-edge search, and monotonic survival probability in b was not tested. The earlier MG6 b = 28 panel used different streams. Report the loss of robustness at the low endpoint without claiming an estimated boundary.

Likewise, dissolution in SH2 establishes that local alignment of the birth bias with the particle is unnecessary for the observed loss. SH2 still has field-dependent births elsewhere on the torus. MG7 omitted DEC, so these data do not isolate uniform birth–death rewriting as the cause. The result is compatible with a narrower robust-survival region under rewriting, but does not separately identify the contribution of shifted coupling, field transfer, and graph evolution.

### 2. Retain the descriptive output that remains defined

The evaluator's early return correctly prevents the primary score, but also suppresses aggregate descriptive sections requested in §5. The high-endpoint contrast and threshold sensitivity were already requested there, even though the separate summarizing script was written after the data. Distinguish that post-data code addition from a newly chosen scientific analysis.

In particular, the assertion that threshold sensitivity was not computed is inaccurate. The evaluator computes it in run_row and stores it under each run's alt_threshold field before returning. I independently reproduced those counts. Bring the resulting unchanged survival verdict into the closing summary.

Show the available early/late means rather than implying a settled response. For CP2 at b = 28.5, late minus early D₀ means are +19.4, −6.0, +13.8 and +16.4 nodes. For SH2 they are +19.2, −15.0, −20.6 and +5.4. The corresponding changes in the paired contrasts are +0.2, +9.0, +34.4 and +11.0 nodes. These are descriptive temporal changes over the fixed window, not a new gate or proof of secular drift. They matter when interpreting a small difference between b = 28 and 28.5 as an equilibrium scaling effect.

The low-endpoint contrast, endpoint residuals, endpoint ratio and complete twelve-run regression have undefined inputs and should remain unavailable. Report that explicitly. Any retained high-endpoint profile can be shown separately.

## High-endpoint result

From integer disc counts, the four CP2-minus-SH2 contrasts at b = 28.5 are

−75.0909, −71.0000, −74.8182 and −67.1818 nodes.

The mean is −72.0227 nodes and the conditional paired t is −38.6320, agreeing with descriptive.json. This is a replicated local deficit among four surviving paired starts. The high-endpoint minus MG6-midpoint differences are negative on all four starts, as reported, but the experiments use different streams and the planned low endpoint is missing. Keep that cross-study comparison descriptive and retain the registered NOT SCORED outcome.

## Scope and files

This review is a records audit. I did not rerun the ODE evolution or independently rebuild the full-field quadratures from PC-only snapshots. Agreement between the recorded quadrature implementations was checked from their saved outputs. Successful numerical validation on the unregistered Part 0 graph does not establish indefinite survival on every registered graph.

The replication script, output and strict replay log are under `review/MG025_replication/`. My recommendation is to sign off MG7 with the recorded feedback outcome and the interpretation corrections above. The mass-scaling question remains open.
