# MG9b — result

*Claude, 2026-11-05. Scored by the registered rules (draft 3, frozen at `185aceb`, with amendment 1) from the
committed records in `results/MG9b/` (run on the PI's PC, commit `5a54f7a`). Awaiting ChatGPT result review and PI
sign-off.*

## Result

**NOT RESOLVED** (confirmatory). Branch B3 (all three b), calibrated. No qualifiers: no LINEAR MODEL MISFIT, no
NON-MONOTONE, no LATTICE-FLAGGED, CAPPED, BUDGET-SHORT or BUDGET-TRUNCATED.

The trend in the activation slope is S₁ = **−0.023** per unit Ω per unit b, 95% profile CI **[−0.094, +0.048]**. The
interval contains 0, so no change of S with b is shown over b = 29.00–29.35. The design had limited power for modest trends: the ordinary calibration
scenario gives about 0.73 at S₁ = −0.10 and 0.25 at S₁ = −0.05 (scenario-specific design properties; no direct estimate
at −0.023). *(Wording corrected after ChatGPT MG035.)*

| b | Ω cells (runs / events) | S (per unit Ω) | 95% profile CI | Lack of fit p |
|---|---|---:|---|---:|
| 29.00 (MG3b, as recorded) | 130, 120, 110, 100 | 0.0642 | 0.0465–0.0817 | 0.89 |
| 29.20 | 170 (32/28), 155 (32/22), 140 (48/25) | 0.0479 | 0.0284–0.0672 | 0.08 |
| 29.35 | 250 (64/26), 235 (32/21), 220 (32/29) | 0.0570 | 0.0396–0.0745 | 0.68 |

Trend fit S_b = S₀ + S₁(b − 29.175): S₀ = 0.0571. Linear-model check (three separate slopes vs the linear trend): LR =
1.09, 1 df, p = 0.30.

## Controls (all PASS)

| b | V-b1 reference | V-b2 late window | V-b3 tight tolerance |
|---|---|---|---|
| 29.20 | PASS (drift 1.1e-15) | PASS, 10/10, all complete | PASS, 31 confirmations, all complete, none changed |
| 29.35 | PASS (drift 2.8e-14) | PASS, 10/10, all complete | PASS, 33 confirmations, all complete, none changed |

## Descriptive (not scored)

- **Placement.** At 29.35, Part A stopped at the first grid point (Ω = 250, 4/4 events), so Ω_T sits at the top of the
  registered grid. The Part B cells (250, 235, 220) are all inside the stop20 rule and none is capped. At 29.20, Part A
  gave 0/4 at 250 and 200, then 3/4 at 170.
- **Rates.** At comparable rates, 29.35 needs Ω about 65 higher than 29.20 (k ≈ 8e-5 at Ω = 220 vs Ω = 155), as expected near the upper edge of the
  window (b_r ≈ 29.396 from MG9).
- **Startup.** Runs with a first event before t = 50 add no exposure (registered rule). Counts: 29.20 cells 170/155/140
  have 1/10/23 startup runs; 29.35 cells 250/235/220 have 7/6/2. Risk sets and landmarks at t = 100 and 150 are in
  `evaluation.json`.
- **Event shape.** 40 of 75 first events at 29.20 and 30 of 76 at 29.35 had a single deep minimum at t = 300 in their
  confirmation. The rest had 2–5 persistent minima at confirmation (histograms in `evaluation.json`).

## Deviations and cost

- **Ledger 26.28 core-h** (cutoff 50, no budget stop; projection about 25–35). Rate collection (Parts A and B with
  pilots and screens) used about 15.7 core-h.
- **Two very slow tight-tolerance confirmations.** V-b3 at 29.20 used 10.24 core-h. Two confirmations in cell 170 took
  19,312 s (rep 5) and 15,697 s (rep 22); the median was 18 s. Both completed and agreed with the standard solver. The
  run's wall time was about 8 h against the 4–6 h estimate. V-b3 has no per-job time limit, and the registration set
  none. This did not affect any scored quantity. A time limit (with an over-limit job counted as incomplete) should be
  set in any future study that uses tight-tolerance confirmations.
- No other deviations. All 280 run records have b_det = b. All 280 runs reached their planned endpoints; all 16
  placement and 240 production records are complete for first-passage scoring. The 8 T = 20 clip pilots (4 in Part A, 4 in the Part B screens) and 15 of the 16
  T = 50 crowding-screen runs carry `complete=false` under the full first-passage definition, as intended (the other
  screen run had a startup event); none entered a rate fit. *(Corrected after ChatGPT MG035.)*

## Independent recomputation

`src/check_mg9b.py` recomputes the result from the raw run JSONs without going through `mg9b.py` or the evaluator. It
recounts events and exposure from each run's t_first and startup flag, checks the counts against `partB_*.json` and the
stop20 batch boundaries, and refits with scipy (BFGS, Brent for the profile limits). It matches the evaluator:

- the per-b slopes to 1e-8;
- S₁ (−0.023208 vs −0.023207);
- the CI limits to 1e-10;
- the LR statistic (1.0895).
