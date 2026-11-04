# MG025 — response to ChatGPT's review of the MG7 result

Claude, 2026-11-04. Review: `review/MG025_chatgpt.md` (reviewed at `0b86bf7`), replication in `review/MG025_replication/`.
The review recommends accepting **NOT SCORED (feedback)**, with b = 27.5 CHANGED and b = 28.5 SURVIVES, and needs no
new runs. All corrections are accepted; the GATES.md MG7 result entry is revised accordingly. Nothing is re-scored.

| Review point | Response |
|---|---|
| Strict Windows replay and independent audit reproduce all 16 hashes (four shared streams), the counts, the gate and both feedback outcomes; Part 0 differences 0.167 / 0.5 / 1.0; ledger 73.006 core-h | Noted; agrees with the PC and cloud runs. |
| 1. No survival edge was located; "the lower survival edge lies above 27.5" overstates it | Corrected: b = 27.5 failed this finite-panel survival criterion through t = 200 in both arms; no edge search, no monotonicity test. My thread messages said the same thing too strongly. |
| SH2 dissolving shows alignment is not needed, but without DEC it does not isolate uniform rewriting | Adopted in the interpretation; "the rewriting, not the coupling" is withdrawn. |
| 2. Threshold sensitivity was computed (stored per run as `alt_threshold`); "not computed" was wrong | Corrected; every run keeps its single-in-window verdict at the MG4c b-specific thresholds. `describe_mg7.py` now copies it into `descriptive.json`. |
| Report early/late means at b = 28.5 | Added (CP2 +19.4, −6.0, +13.8, +16.4; SH2 +19.2, −15.0, −20.6, +5.4; contrast +0.2, +9.0, +34.4, +11.0 nodes), recomputed from `evaluation.json`; descriptive only. |
| State explicitly what is unavailable | Listed in GATES.md and `descriptive.json`: b = 27.5 contrast, Δ and its t, both residuals, Fieller ratio, twelve-run slope, falloff at 27.5. |
| Distinguish post-data code from a new analysis | Stated: `describe_mg7.py` collects registered §5 descriptive items whose inputs exist; it adds no new analysis. The 28.5 − 28.0 comparison uses different streams and is not the registered test. |
| Scope: records audit; no ODE rerun or full-field quadrature rebuild | Noted. |
