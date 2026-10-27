# ChatGPT review of MG003b (round 3) — verbatim
Received 2026-10-27 via PI. Reviewed @ d5b381c.

The revision improves the inference layer, but I’m not ready to close the review. I read the response, protocol v1.2 §5.3, all three calibration scripts, their outputs, and `stats_mg2.py`.
Three items remain:

1. [P1] The confidence bounds need explicit distribution assumptions.
The triangle inequality is valid; the t-based radius is not a general confidence guarantee for arbitrary run-profile distributions.
I confirmed that the old balanced example now correctly gives [0, 4.60%] → INCONCLUSIVE. Exact enumeration also supports the original two-point calibration’s approximately 3% false-decision rates.
However, consider independent run errors of 0.5% with probability 0.9, or 6.5% with probability 0.1. The true bias is 1.1%, with run-noise RMS 1.8%—within the tested magnitude ranges. Twelve common outcomes produce a bound of [0.5%, 0.5%], falsely passing the 1% criterion. This happens with probability
\[ 0.9^{12}=28.24\%. \]
I verified this using the actual bounds function and independent simulation.
Required: specify the assumptions supporting the confidence claims and validate them against actual survival-conditioned run averages. Broaden held-out calibration to rarer mixtures and unequal two-sample noise. No fixed sample count supplies a universal guarantee without tail assumptions.
2. [P1] The α-calibration record contradicts the response.
The committed [calibration.txt](https://github.com/RichardBenear/ModelG/blob/main/results/MG003b/calibration.txt) ends with “calibrated widening factor: None”; the JSON likewise contains `"widen": null`. The protocol specifies w = 2.
The driver assigns coverage zero when no intervals are defined, forcing its minimum to zero. Undefined intervals should be reported separately.
The separate decision-calibration output supports the reported 0.020 false-FAIL / 0.003 false-PASS results on its tested scenarios. That is useful evidence, but it is distinct from successful interval calibration.
Required: reconcile the selection rule, script, saved outputs, and wording. If w=2 is calibrated for the decision rule, label it that way.
3. [P2] Validate the complete decision procedure.
The current calibration uses fixed sample sizes. The gates additionally use survival filtering, extensions after inconclusive results, and—for MG3—a random reference normalization.
Required: enforce the minimum number of eligible profiles and calibrate the extension procedure. The two-sample study should also cover unequal variances and the normalization actually used.

The Euler revision remains accepted, as does the U-statistic point estimator. V2’s move to explicit equivalence margins is appropriate.
[Independent audit code, results, and detailed review](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG003_replication/round3/README.md). No remote repository files were changed.
