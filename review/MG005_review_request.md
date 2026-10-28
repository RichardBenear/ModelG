# MG005 — Review request: MG2 gate result (protocol v1.3.3)

For ChatGPT (independent replication and adversarial review). Prepared by Claude, 2026-10-28.

## What to review
- Registered criteria: `GATES.md`, section MG2 (unchanged since pre-registration except the PI budget note, which was recorded before any C2–C5 evaluation).
- Result record: `GATES.md` → "MG2 result"; full output `results/MG2/MG2_evaluation.md` and `results/MG2/evaluation.json`.
- Raw data: `results/MG2/cells/*.pkl` (every attempted run; stripped of fields except C5 confirmation fields). Run log: `results/MG2/run_log.txt`.
- Code: `src/run_mg2.py` (driver), `src/evaluate_mg2.py` (scoring), `src/stats_mg2.py`, `src/measure.py`, `src/confirm.py`, `src/runner.py`.
- Post-gate diagnostic: `src/diag_mg2_replication.py` → `results/MG2/diag_replication.txt`.

## Summary
Registered verdict **FAIL**: C3 FAIL (one cell: b = 29, Ω = 100, X variance +17.1 %), C4 FAIL (replication at b_mid in 12/24 runs at Ω = 300; widespread events at Ω = 100). C1, C5, C6, N1, N2 PASS; C2 INCONCLUSIVE; pooled A1 audit ok.

## Specific questions
1. **Scoring.** Re-score C1–C6, N1, N2 independently from the cell files. Does any decision change?
2. **C2.** No ladder cell has a bias lower bound above 0, so the selection rule selects nothing. Point estimates (X, b = 29): 1.80 % (Ω = 1000), 0.46 % (3000), 0 (10⁴) with upper bounds 4.8 %, 2.3 %, 1.1 %. Is INCONCLUSIVE the correct registered outcome, and is the cause correctly stated as insufficient power?
3. **C3 at Ω = 100.** Is attributing the failure to far-field contamination by births and replications (24 births, 31 replications in that cell) supported, or should the LNA's leading-order error at Ω = 100 be considered first? The X excess decreases as +17.1, +6.2, +4.2, −0.9 % over Ω = 100, 300, 1000, 3000.
4. **C4 replications.** The diagnostic re-runs three replicating runs bit-identically and continues the field deterministically to t = 300: two persistent particles each. Is this sufficient to call the events physical rather than detector artifacts?
5. **Design error.** b_mid = 29.16 is 0.24 below the h = 0.2 replication edge (29.40). Do you agree the C4 failure reflects a narrow noisy window rather than a defect of the token rules' mean-field convergence?
6. **Next step.** Options for the PI: (a) accept the FAIL record and proceed to MG3 unchanged (MG3 is at b = 29 and uses field-level confirmation, so it does not rely on b_mid); (b) an amendment re-testing persistence at a b re-chosen for the noisy window — necessarily post hoc and labelled so. Which do you recommend, and does anything in this result require an MG3 change before it runs?
