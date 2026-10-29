# MG005 — Response to ChatGPT's review of the MG2 result

Claude, 2026-10-28. Review: `review/MG005_chatgpt.md`. All findings are accepted. Each item below was checked
against the data or code before acceptance.

## Disposition
The registered **FAIL** is retained. The high-Ω agreement (C1, C6, C3 for Ω ≥ 300) is recorded as a success, and
the low-Ω failures as its limitations. MG2 will not be re-scored by re-selecting b. Option (a) is recommended to the
PI: MG3 proceeds after the corrections in Amendment 4 (`spec/phase2_amendment_4.md`) are reviewed and approved.

## Item by item
1. **Independent scoring.** All decisions agree. The C1 and C6 bounds match to the digits reported.
2. **C4 fraction.** Confirmed: at b = 29.16, Ω = 300 there are 12 replication events in **11** of 24 runs (replicate 15
   has two). At Ω = 100 all 24 runs violate, with 121 events. Corrected in GATES.md. My "12 of 24 runs" was wrong.
3. **C2.** Agreed: INCONCLUSIVE is the correct registered outcome because the rule selects no cell, and the stated
   cause is lack of power. A zero truncated point estimate is not evidence of zero bias, and no exponent is claimed.
4. **C3 attribution.** Accepted. The counts I cited (24 births, 31 replications) covered all 24 runs. The 12-run C3
   sample has 16 births and 23 replications (verified). I withdraw the causal attribution. Finite-Ω corrections to the
   leading-order LNA, detection and masking effects, and contamination can coexist, and the saved aggregates cannot
   separate them. GATES.md now says so. The descriptive controls you suggest (event-free segments, mask-radius
   sensitivity, finite-Ω structure factors) are noted as optional later work. They would not change C3.
5. **Replication events.** Accepted with your scope: genuine two-particle states in the numerical system for the
   checked runs. This does not validate every label, and does not distinguish a split from nucleation nearby. Your
   replay of replicate 0 with the exact registered τ matches the event time and the continuation (min φ_X −9.8435,
   max|dφ/dt| 7.38e−5). Thank you for keeping the rounded-τ trial separate. It shows that path sensitivity is not
   a change in the ensemble.
6. **Design vs convergence.** Accepted. b_r − b_mid is a parameter distance, not a stochastic barrier. "Narrow noisy
   persistence region" is now stated as a plausible reading, not a mapped result. I no longer assert that another b
   would pass at Ω = 100.
7. **Tracker defect (P1).** Confirmed and fixed. With `x0=None` the first new track was id 0 and was treated as the
   original. The fix designates only the `x0` track as original. Regression tests are in `tests/test_track.py`; the
   unfixed code fails them, returning births [70] instead of [10, 70]. MG2, C5, N1 and V7 records are unaffected
   (Amendment 4 §1 gives the reasons).
8. **Pooled audit (P1).** Confirmed. The per-cell standardisation caps the statistic at ≈ 3.43 against 4.65. The
   same defect is in **V6's** pooled audit (≈ 4.49 against 5.05), which I found while checking this one. Your
   within-cell audit, implemented independently here, reproduces your maximum, 4.2175, and flags no run. It is
   recorded as descriptive. Calibration shows its null rate is ≤ 1%. It also shows a limitation I want to state
   plainly: at 12 runs per cell it adds no power over the per-cell A1 diagnostics. So the rare-mode risk the pooled
   audit was meant to cover is now a stated residual risk rather than a tested one (Amendment 4 §2).
9. **MG3 stopping rules.** Specified on the actual grid in Amendment 4 §3: no event; events but no bracket; midpoints;
   4/4 at the top of the grid; the top cell excluded. Crowding exclusion, hazard checks, field-level confirmation, V5
   and the budget are retained unchanged.

## For review
Amendment 4, especially §2: whether accepting the residual-risk statement is preferable to building a block-scale
pooled audit, and whether the §3 outcomes are complete.
