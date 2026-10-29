# MG007 — ChatGPT review of the MG3b Part 3 calibration

Reviewer: ChatGPT. Relayed verbatim by the PI in the project thread, 2026-10-29. Reviewed: `results/MG3b/calibration/`
and `src/calib_mg3b.py` on `main` after the merge of PR #4 (commit `18d438f`). The reviewer's supporting scripts are on
the PI's machine (`review/MG007_replication/` in ChatGPT's workspace) and are not yet in the repository.

---

**The calibration is useful, and the registered fallback is handled correctly. But I found two remaining gaps in the post-selection analysis.** The standalone results should not yet be described as validating every interval used after timestep selection.

I checked the implementation and ran its ten tests; all passed. I also performed independent synthetic checks.

**1. The reused Ω = 120 rate needs its own coverage qualification.**

Garwood was calibrated on standalone cells. The full decision-tree calibration checks the *slope*, but not the rate interval for the arm selected for reuse.

In 100,000 independent simulated studies with equal true rates and 35% startup loss in each arm, the selected rate’s Garwood coverage was:

| Production branch | Coverage |
|---|---:|
| τ after EQUIVALENT | 98.6% |
| τ/2 after INCONCLUSIVE | **90.8%** |
| τ/2 after Part 1b EQUIVALENT | **84.8%** |

These are **branch-conditional** coverages—the standard already applied to the slope. Averaging over branches answers a different question.

Therefore, the standalone minimum coverage of 95.0% cannot be transferred to the selected Ω = 120 rate within each branch. Add this diagnostic and, where necessary, label that interval approximate. That option requires no additional lattice runs.

**2. Part 1b’s equivalence error is not covered by the standalone comparison.**

Part 1b reuses τ/2 data after those data helped trigger DIFFERENT. Its comparison is consequently selected, rather than an independent repetition of Part 1.

For an illustrative case with first ratio 0.71 and second ratio 0.5, my simulation found:

- Part 1b declared EQUIVALENT in **9.06% of triggered comparisons**.
- The whole-study probability of triggering Part 1b and then declaring EQUIVALENT was **2.285%**.

This does **not** demonstrate a whole-study error above 5%. It shows why the claimed error rate must distinguish conditional from whole-study control. Check the triggered comparison at both equivalence boundaries before assigning it the standalone calibration’s interpretation.

The rest of the report is appropriately cautious:

* Keep the approximate profile-likelihood slope interval; do not use diagnostic reruns to reverse the registered failure.
* Describe 0.937 as an **observed minimum over qualifying calibration cells**, not a guaranteed coverage floor.
* The median-split test’s low power means non-rejection will usually provide little evidence for stationarity.
* The disclosed implementation deviations should remain explicit.

One smaller implementation mismatch: hypothetical 96-run tree simulations also use 96 τ/4 runs, whereas Part 1b specifies 64. The selected actual design is 64, so that mismatch does not affect its current branch simulation.

**My recommendation:** retain the existing method decisions, then add the selected-rate and triggered-comparison checks or qualifications before calling the production analysis fully calibrated. This does not require reopening MG3 or redesigning the entire study.

The supporting review and scripts are saved in [MG007 replication review](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG007_replication/README.md).
