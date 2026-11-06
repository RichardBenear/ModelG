# MG8 results: scheduler sensitivity of the lattice approximation

*Stage MG038. Claude, 2026-11-06. Registered design `spec/MG8_proposal.md` draft 4 (frozen `3ab4ec8`); pre-run spec
`spec/MG8_prerun.md` (deviation D1 accepted by the PI). Production ran in the cloud (Linux, 4 cores, 4 workers),
2026-11-05 23:37 to 2026-11-06 00:31 UTC: all 55 records (the 57 registered runs), **3.09 core-h** against the 20 core-h
cutoff, state `complete`. Evaluator `src/evaluate_mg8.py` → `evaluation.json`; independent recompute
`src/recompute_mg8.py` → `recompute.json` agrees on all 106 compared quantities and labels.*

## Registered outcomes
- **Implementation checks (§4.1): all PASS.** Part 0 (accepted attempt 2); deterministic L/R mirror at k = 0 (t = 50…300)
  and k = 6 (t = 50), max relative difference 0 (bit-identical); stochastic mirror, identical count hashes.
- **Production step, τ₀ (§4.2): SENSITIVE.** All three R seam replicates are FATE CHANGED: the particle is lost at t ≈ 8
  in every replicate (no detection from then on). J keeps it.
- **Refinement (§4.3): NOT RESOLVED BY THIS LADDER.** Reason: eligible levels are k = 3–6 only (k = 0 lost the particle;
  k = 1 and two of three k = 2 replicates are UNSETTLED in velocity). Four levels spanning exactly a factor 8 in τ pass
  the count and span clauses, the slope (p = 1.000, diagnostic lower bound 1.000) and the residuals (≤ 4e-7), and the
  horizon check passes (E₆₀₀ = 0.21557 vs E₃₀₀ = 0.21557). The clause "E at the smallest eligible τ ≤ 0.1 × E at the
  largest" fails: the ratio is 0.125 (exactly 1/8, because E ∝ τ). BELOW MEASUREMENT RESOLUTION does not apply (E is far
  above the floor).

| k | τ/τ₀ | E_k (median) | labels (3 replicates) | summary |
|---|---|---|---|---|
| 0 | 1 | 13.83 | FATE CHANGED ×3 | SENSITIVE |
| 1 | 1/2 | 6.898 | SHAPE/DRIFT UNSETTLED ×3 | UNRESOLVED |
| 2 | 1/4 | 3.449 | UNSETTLED, MATERIAL + NOT DRIFTING, UNSETTLED | SEAM-DEPENDENT |
| 3 | 1/8 | 1.725 | MATERIAL, NOT DRIFTING ×3 | SENSITIVE |
| 4 | 1/16 | 0.862 | MATERIAL, NOT DRIFTING ×3 | SENSITIVE |
| 5 | 1/32 | 0.431 | MATERIAL, NOT DRIFTING ×3 | SENSITIVE |
| 6 | 1/64 | 0.216 | MATERIAL, NOT DRIFTING ×3 | SENSITIVE |

Seam spread S_k = 0 at every level (the three replicates agree to the printed digits); E_floor = 1e-9 and v_floor = 1e-8 at
every level (J is exactly stationary). E_raw equals E to 4 significant figures everywhere; E_nyq ≤ 0.023 (k = 0) and
≤ 4e-5 at k = 6 (deviation D1 changes nothing here).

## B2, stochastic at τ₀, Ω = 1000 (descriptive)
| arm | runs with loss | 95% CI (exact) | loss times | eligible |
|---|---|---|---|---|
| J | 7/8 | 0.47–1.00 | 102–292 | 1 |
| R | 8/8 | 0.63–1.00 | 3–20 | 0 |
| P | 8/8 | 0.63–1.00 | 5–9 | 0 |

No replication, merge or TRACK-FAIL anywhere; no clipping. With fewer than 5 eligible runs per arm, drift and stochastic
shape comparisons are NOT ESTIMATED, as registered.

## What drives E (descriptive, not registered)
E is dominated by the G species (the per-species errors at k = 3 are G 1.72, X 0.036, Y 0.49), and the R − J difference is
mostly a **uniform offset of the background**: R's median field sits below J's by 0.0324 (G, X) and 0.0466 (Y) at k = 3,
and by 0.0041 / 0.0058 at k = 6. A check run after the evaluation (homogeneous start at J's steady state S0, no particle,
60 units) shows the same thing directly: under R the homogeneous state itself moves to S0 − (0.2593, 0.2592, 0.3726) at
τ₀ and to S0 − (0.0324, 0.0324, 0.0466) at τ₀/8, exactly in proportion to τ; under J it stays at S0 to 4e-15. This is
the reaction re-exposure effect of A4 (in a sweep a node's snapshot already includes this step's arrivals from the
neighbour processed before it), acting on the steady state. Because ‖φ_G‖ is small (0.42), a uniform G offset dominates
the relative error. At τ₀ the offset is large enough that the particle dissolves, and the R field relaxes to the shifted
homogeneous state; the stochastic R and P runs lose the particle within 3–20 units, consistent with the same effect.

The particle itself is otherwise unchanged by the sweep at k ≥ 3: single and settled, aligned shift |σ*| ≤ 4e-4, and no
drift (|v| ≤ 3e-8, under the 1e-6 threshold).

## Reading
At the registered production step the update order matters a great deal for this lattice approximation: a rightward
sweep shifts the homogeneous state by an amount proportional to τ, which at τ₀ is enough to destroy the particle, while
the synchronous kernel keeps it. As τ is refined, the deterministic difference falls exactly in proportion to τ
(first-order convergence, slope 1.000), but under the registered rule the ladder does not resolve convergence: the
eligible levels span only a factor 8, so the frozen "≤ 0.1" clause fails at 0.125. Central Question 5 stays open, as
registered.

## Files
`runs/` (55 records), `fields/` (B1 checkpoints), `snap/` (B2 φ_X every 5 units, block profiles), `ledger.json`,
`status.json`, `evaluation.json`, `recompute.json`. The full B2 counts (`full/`, 34 MB, gitignored) were written in the
cloud container and are not kept; their SHA-256 values are in the records, and the runs can be regenerated from the
seeds.
