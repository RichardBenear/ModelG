# MG3b Part 3 — calibration results (for review before any lattice run)

Claude, 2026-10-29. Registered design: `GATES.md` (MG3b) and `spec/MG3b_followup_proposal.md` (Part 3, frozen at
`202c00f`). Code: `src/calib_mg3b.py`. Tests: `tests/test_calib_mg3b.py` (10 pass). Synthetic only; no lattice simulation.
Platform: cloud container (Linux, Python 3.13.17, numpy 2.5.3, scipy 1.18.1), 4 workers.

## Conclusions first
| Item | Decision (`decisions.json`) | Basis |
|---|---|---|
| Per-cell rate CI | **(a) Garwood** | The only eligible candidate. Minimum coverage over the 51 constant-hazard configurations is 0.950. The profile likelihood (0.938) and the bootstrap (0.922) fail |
| Part 1 rate-ratio CI | **(a) conditional binomial, Clopper–Pearson** | The only eligible candidate. Mid-p and Wald fail the 90% tail requirement and the equivalence-error requirement (P(EQUIVALENT) at the margin up to 0.052–0.056, against ≤ 0.05) |
| Part 1 runs per arm | **64** | P(EQUIVALENT \| r = 1) = 0.85 under startup losses (0.35, 0.05), above the 0.80 rule. It is 0.75 if τ/2 has the same startup loss as τ (0.35, 0.35) |
| Activated-slope CI | **No candidate eligible.** Fallback: **(b) profile likelihood**, reported as "approximate (observed minimum coverage 0.937 over qualifying calibration cells; not a guaranteed floor)". No claim about the slope may be stated at 95% | See "Slope" below |
| Stationary-check split | **(b) the cell's exposure median** | Split (a) at t = 200 is not eligible: its size is 0.0143 > 0.012 in one configuration. Split (b) has a maximum size of 0.009. Its power is low (see below) |
| Bootstrap replicates | R = 4000 (no fallback needed) | The benchmark projected 0.60 core-hours, against the 8-hour ceiling |

Cost: 0.62 core-hours for the main run, plus 3 benchmarks of about 0.03 h each and about 0.03 h of diagnostics. The total
is about 0.74 core-hours, well under the 8 core-hour Part 3 ceiling (ledger: `results/MG3b/budget_ledger.csv`).
The benchmark ran three times; the second and third runs followed a fix to the GLM fit, below.

## Slope: why no method is eligible
- Coverage averaged over the 57 (configuration × branch) cells that have ≥ 1000 replicates is good: 0.952 for both Wald
  and profile.
- Within particular production branches, the tails are lopsided. Wald fails in 12 cells and profile in 11 (out of 57).
  The worst cell for each is `cap_ps0.6`, 64 runs, r = 1, branch τ/2-after-INCONCLUSIVE: coverage 0.936 and 0.937. The
  usual failure is a one-sided miss above 0.035, up to 0.047, with the true slope above the interval (Ŝ too large). It
  is concentrated in branches entered by selection (r = 0.71 → τ; r = 1 → τ/2 after INCONCLUSIVE; Part 1b) and in the
  capped configurations. This is the selection effect the review asked to check.
- Monte Carlo noise contributes. One failing cell (S = 0.073, 96 runs, r = 0.71, branch τ) was re-run after the main run
  with three new seeds at R = 8000. Wald gave coverage 0.951–0.956 with an upper miss of 0.034–0.035; the registered run
  had 0.941 and 0.044. An independent implementation (`independent/indep_branch_check.py`, which uses scipy's binomtest
  and a BFGS Poisson fit) gave 0.955 with an upper miss of 0.032 (n = 2973 in the branch). So that cell's true upper miss
  is about 0.034, right at the 0.035 limit. With 57 cells, several near the limit will fail by chance.
- **The registered rule decides the outcome, and these re-checks do not change it.** I did not re-run the calibration to
  try for a pass. The registered consequence is that the slope is reported with the profile-likelihood interval,
  labelled approximate, with calibrated minimum coverage 0.937. Any change, such as a larger R, different tolerances or a
  different interval, would be a PI decision recorded as an amendment.
- Bootstrap (c): minimum coverage 0.927. It fails in 39 cells, mostly by under-covering the upper tail.

## Part 1 operating characteristics (Clopper–Pearson, 64 runs per arm)
| r = k(τ/2)/k(τ) | 0.5 | 0.71 | 0.8 | 1 | 1.25 | 1.41 | 2 |
|---|---|---|---|---|---|---|---|
| P(EQUIVALENT), p_s (0.35, 0.05) | 0.035 | 0.42 | 0.64 | 0.85 | 0.64 | 0.44 | 0.038 |
| P(DIFFERENT) | 0.83 | 0.29 | 0.15 | 0.041 | 0.16 | 0.34 | 0.90 |
| P(EQUIVALENT), p_s (0.35, 0.35) | 0.042 | 0.37 | 0.54 | 0.75 | 0.59 | 0.39 | 0.039 |

The 95% coverage is 0.955–0.964 and the 90% coverage is 0.910–0.926 in all 28 configurations. The full table, including
96 runs per arm, is in `part3_results.json` under `part1`.

## Stationary check: the design can rarely diagnose time dependence
- Power of the chosen median split at α = 0.01: mean 0.085 against a ×2 step at t = 200 (max 0.54), 0.058 against ×0.5
  (max 0.30), 0.14 and 0.23 against Weibull shapes 0.7 and 1.5, and 0.06 against an early excess.
- The fixed split has more power against steps (mean 0.20, max 0.71), but it is slightly anticonservative in one
  high-rate, high-loss configuration.
- Consequence: in Part 2, the "informative" label (simulated power ≥ 0.5 against both ×2 and ×0.5) will usually not be
  met, and the check will mostly be reported as uninformative. That is a finding about this design, stated in advance.
- Under the non-constant hazards, k̂ differs from the exposure-weighted mean hazard by −0.1% to +9.7%.

## Other results
- **Deviance lack-of-fit (df 2):** size at nominal 0.05 is 0.037–0.065 across the activated configurations. Its power
  against Ω = 130 at a third of the activated rate is 0.78, with Ŝ shifting from 0.073 to 0.103.
- **Capped configurations work as intended.** Ω = 130 capped below 20 events in 91% (prefactor ÷ 3), 86% (p_s 0.6) and
  100% (both) of replicates. The mean slope estimate stays close to its generating value of 0.073: mean Ŝ is 0.074 in each case. No bias estimate with Monte Carlo uncertainty is claimed. [Wording corrected per `review/MG007_replication/README.md` §3.]
- **Branch frequencies** (S = 0.073, 64 runs, r = 1, tree with both arms at p_s 0.35): τ 75%, τ/2 after INCONCLUSIVE 22%,
  τ/2 after Part 1b 1.6%, stop 2.0%.

## Deviations and corrections (stated plainly)
1. **Implementation choices were not appended to `GATES.md` before the run.** The registration says details the code
   needs should be appended before the run they affect. I made them in the code before running but documented them only
   afterwards. They are listed in `GATES.md` (MG3b, "Part 3 implementation choices", dated 2026-10-29) for review.
   None was chosen after seeing results, but that can't be shown from the record, so the reviewer should judge them on
   their merits.
2. **GLM fit fix before the main run.** The first benchmark showed undamped-Newton divergence in a few replicates
   (16 in 20,000 at S = 0.10). I added step halving and a step cap before the main run, and re-ran the benchmark. The fit
   now converges in every replicate checked (score below 2e-12). This was found from warnings, not from coverage results.
3. **Branch-threshold correction after the main run.** For the bootstrap, the code scaled the 1000-replicate threshold by
   its live-replicate count, giving about 995 instead of 1000 at R = 4000. Corrected; `python src/calib_mg3b.py recheck`
   recomputes eligibility from the saved results. The selection is unchanged: none eligible, (b) chosen
   (`eligibility_recheck.json`). `decisions.json` is left exactly as the run wrote it.

## Addendum 1 (2026-10-29): post-selection checks requested in `review/MG007_chatgpt.md`
This addendum was registered in `GATES.md` before it was run. It covers Part 1 and Part 1b only, with Clopper–Pearson
and 64 runs per arm (τ/4 included): 70 configurations at R = 20 000, costing 0.005 core-hours. Results are in
`addendum1/addendum1_results.json`.

**(A) The reused Ω = 120 rate: its Garwood interval is labelled approximate in every branch.**

| Production branch | Qualifying configurations | Observed minimum coverage | Largest one-sided miss | Failing configurations |
|---|---|---|---|---|
| τ after EQUIVALENT | 50 | 0.943 | 0.057 (lower), 0.040 (upper) | 19 |
| τ/2 after INCONCLUSIVE | 70 | 0.875 | 0.090 (upper) | 30 |
| τ/2 after Part 1b EQUIVALENT | 31 | 0.841 | 0.159 (lower) | 15 |

- Even in the τ branch, where coverage is high (0.987 at r = 1), the tails are lopsided whenever selection depends on
  the arm's own count. For example, at r = 1.41 the τ arm is chosen when its count happens to be high.
- These results reproduce ChatGPT's figures. At r = 1 with 35% startup loss in every arm, I get 0.987, 0.908 and 0.862,
  against ChatGPT's 98.6%, 90.8% and 84.8%. The last branch has only 369 replicates there, so it does not qualify.
- Consequence (registered rule): the Ω = 120 rate interval is reported as "approximate (observed minimum
  branch-conditional coverage X)", with the X of whichever branch production lands in. No 95% claim is made for it. Its
  data enter the slope fit unchanged; the slope's own calibration already covers selection.

**(B) Part 1b equivalence error: within 5% for the whole study, but not for triggered comparisons.**
- At the boundaries r2 ∈ {0.5, 2}, the whole-study probability of triggering Part 1b and then wrongly declaring
  EQUIVALENT is at most **0.040** (r = r2 = 2). That is ≤ 0.05, so the registered rule adds no label.
- The **conditional** error among triggered comparisons is much higher, up to **0.161** (among configurations with ≥ 100
  triggered replicates). In ChatGPT's illustrative case (r = 0.71, r2 = 0.5), I get 9.6% conditional and 2.4%
  whole-study, against ChatGPT's 9.06% and 2.285%.
- So a Part 1b EQUIVALENT will be reported as "whole-study error ≤ 0.040 at the margin; conditional on Part 1b being
  triggered, up to 0.161". It will not carry the standalone 5% interpretation.

**Scope restriction (per `review/MG007_replication/README.md` §3):** in the slope's decision-tree calibration, all
time-step arms share the Ω = 120 cell's startup probability. Different startup losses between arms, such as
(0.35, 0.05), were calibrated for Part 1 and in addendum 1, but not for the slope. The slope calibration therefore
covers equal startup loss across arms only.

**Implementation correction:** the τ/4 arm now always has 64 runs. The main run's 96-run tree configurations used 96.
Those are hypothetical (the selected design is 64), and the main-run results are left as they are.

## Files
- `bench.json`: benchmark (the last of the three runs).
- `part3_results.json`: every configuration.
- `decisions.json`: the decisions the driver will read.
- `eligibility_recheck.json`: slope eligibility recomputed with the corrected branch threshold.
- `independent/indep_branch_check.py`: independent re-check of one failing cell.

Reproduce: `MODELG_WORKERS=4 python src/calib_mg3b.py bench`, then `run`, then `recheck`. Seeds come from
SeedSequence(20261029); results do not depend on the worker count.
