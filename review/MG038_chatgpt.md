*ChatGPT result review of MG8 (commit `42bebf0`), pasted by the PI 2026-11-06 01:10 UTC. Saved verbatim.*

**I support accepting MG8’s registered results: SENSITIVE at the production step, and NOT RESOLVED BY THIS LADDER under refinement.** I found no result-blocking scoring error. The results narrative needs a few corrections before closure.

I reviewed commit [`42bebf0`](https://github.com/RichardBenear/ModelG/tree/42bebf0db3a3464195db0f403980d04ab1cebeca), including the implementation fixes, production records, field checkpoints, evaluator, and recomputation code. No repository files were changed.

Verification supports the reported outcomes:

- Re-running the repository’s separate recomputation reproduced all **106 compared quantities and labels**.
- My additional checks reproduced **187 field-checkpoint detection counts**, both deterministic mirrors exactly, all stochastic loss classifications from saved detection lists, the binomial intervals, and the **3.0881 core-hour** ledger total.
- The endpoint and implementation-check gating fixes are present. Nyquist filtering is recorded as PI-accepted deviation D1.
- The stochastic mirror’s stored hashes agree; discarded full counts prevent independently rehashing those arrays.

| Registered result | Review |
|---|---|
| Production-step sensitivity | Supported: all three deterministic R replicates lose the particle; J retains it. |
| Refinement unresolved | Correctly scored: eligible levels are k=3–6; their error ratio is **0.12500025**, exceeding the required ≤0.1. |
| Stochastic drift/shape comparisons | Correctly NOT ESTIMATED: eligible counts are J=1, R=0, P=0. |

The eligible refinement levels have slope **0.99999907**, consistent with first-order decay over the measured range. That descriptive evidence must remain separate from the registered unresolved label. [Evaluation](https://github.com/RichardBenear/ModelG/blob/42bebf0db3a3464195db0f403980d04ab1cebeca/results/MG8/evaluation.json)

1. **Correct the interpretation of the J-control losses.**

   MG8 observes stochastic runs to **T=500**, whereas MG2 used **T=100**. All seven J losses occur after time 100, at approximately 102–292. The results therefore do not contradict the earlier short-horizon observation.

   The planning assumption that Ω=1000 would keep J single was insufficient for this longer horizon. Report that explicitly. R/P losses occur much earlier—3–20 and 5–9 respectively—but the final loss counts alone, 7/8 versus 8/8, provide little separation. Earlier loss is descriptive evidence; MG8 did not register a survival-time comparison. [MG2 setup](https://github.com/RichardBenear/ModelG/blob/42bebf0db3a3464195db0f403980d04ab1cebeca/GATES.md), [MG8 results](https://github.com/RichardBenear/ModelG/blob/42bebf0db3a3464195db0f403980d04ab1cebeca/results/MG8/MG8_results.md)

2. **Replace “the particle itself is otherwise unchanged.”**

   Single, settled, and below the drift threshold do not establish unchanged shape.

   I confirm that uniform offsets dominate E. However, after subtracting the spatial mean difference at the scored alignment, replicate 0 still has normalized Y discrepancies of approximately **2.56% at k=3** and **0.31% at k=6**. These are additional descriptive calculations, not replacement scores.

   A supported statement is: “The particle remains single and settled at k≥3, with no drift resolved by the registered criterion; background offsets dominate the measured discrepancy.” At k=6, the registered E remains **0.216**, so refinement has not made that observable negligible. [Results interpretation](https://github.com/RichardBenear/ModelG/blob/42bebf0db3a3464195db0f403980d04ab1cebeca/results/MG8/MG8_results.md)

3. **Qualify the causal explanation and preserve its evidence.**

   The committed fields support the reported background offsets: for replicate 0, approximately −0.03241 in G/X and −0.04668 in Y at k=3, halving with refinement.

   They do not by themselves establish that the background shift is sufficient to destroy the particle, or uniquely attribute the effect to reaction re-exposure. Those are plausible interpretations of coupled scheduler effects.

   The post-hoc homogeneous-start diagnostic described in the report has no corresponding committed diagnostic script or output among the reviewed MG8 files. Archive those inputs and outputs if retaining its numerical claims; otherwise identify the explanation as provisional. This does not require rescoring production.

4. **Correct seam-spread wording and record the retention deviation.**

   Seam spread is **small, not zero**: the saved values range from approximately \(6.47\times10^{-8}\) to \(2.84\times10^{-6}\). Replace the exact-zero statement.

   Also record explicitly that the full stochastic count arrays were discarded, departing from the planned retention arrangement. Saved detections and deterministic fields support substantial independent checking, but hashes cannot replace the missing arrays. Regeneration is a separate reproducibility step. [Results and retention statement](https://github.com/RichardBenear/ModelG/blob/42bebf0db3a3464195db0f403980d04ab1cebeca/results/MG8/MG8_results.md)

**Disposition:** accept the numerical results and preserve the registered labels. Make these reporting corrections before closure; no additional production runs are needed. Central Question 5 remains open.
