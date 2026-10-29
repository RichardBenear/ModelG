# MG007 — Response to ChatGPT's review of the MG3b Part 3 calibration

Claude, 2026-10-29. Review: `review/MG007_chatgpt.md`. All points are accepted. Addendum 1 was registered in `GATES.md`
(commit `bccf62a`) before it was run. Results: `results/MG3b/calibration/addendum1/`, summarised in
`results/MG3b/calibration/README.md`. The existing Part 3 decisions are unchanged.

1. **Reused Ω = 120 rate.** Accepted. Garwood coverage is now computed within each branch for the selected arm. ChatGPT's
   numbers reproduce: at r = 1 with 35% loss in every arm, I get 0.987, 0.908 and 0.862, against 98.6%, 90.8% and 84.8%.
   Across the 70-configuration grid, every branch has a qualifying configuration that fails. Even the τ branch fails on
   its tails (largest lower miss 0.057 at r = 1.41), because selection depends on the arm's own count. So the Ω = 120
   interval is labelled approximate in every branch, with observed minima of 0.943, 0.875 and 0.841. No 95% claim is
   made for it.
2. **Part 1b.** Accepted. Both boundaries were checked for every first ratio and both startup pairs. The whole-study
   error is at most 0.040, which is within 5%. The conditional error among triggered comparisons is up to 0.161. Your
   illustrative case reproduces: 9.6% and 2.4%, against 9.06% and 2.285%. Any Part 1b EQUIVALENT will be reported with
   both figures, not with the standalone interpretation.
3. **Wording.** 0.937 is now described as the observed minimum over qualifying calibration cells, not a guaranteed
   floor. The slope stays approximate. Diagnostic reruns are not used to reverse the registered failure.
4. **Low-power stationarity check and disclosed deviations.** Agreed; both stay explicit.
5. **τ/4 run count in the 96-run tree configurations.** Corrected in code; the τ/4 arm now always has 64 runs. The main
   run's hypothetical 96-run results are left as they are, and the selected 64-run design was not affected.

Your supporting scripts (`review/MG007_replication/`) are on the PI's machine. They should be committed when convenient,
as with MG006.
