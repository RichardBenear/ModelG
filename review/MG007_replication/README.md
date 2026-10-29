# MG3b calibration review — supplementary selection checks

ChatGPT, 2026-10-29. Reviewed `results/MG3b/calibration/README.md`, `src/calib_mg3b.py`, the frozen proposal, tests, `decisions.json` and `eligibility_recheck.json`. This is not a rerun of the entire calibration, and no lattice trajectories were simulated.

## Assessment

The reported selection of Garwood for standalone rates, Clopper-Pearson for the first timestep comparison, 64 runs per arm, the exposure-median stationary split, and approximate profile likelihood for the slope follows the supplied result summaries and the registered fallback. The report correctly preserves failure of all slope candidates rather than treating diagnostic reruns as a new eligibility attempt. Low stationarity power is explicitly acknowledged. These are useful calibration results.

Two post-selection quantities still need explicit checks or qualifications: the reused Omega=120 per-cell rate, and the Part 1b comparison that reuses a selected tau/2 sample. They should not inherit the standalone calibration without checking the conditioning involved. This is a limitation of the scope of the completed calibration, not a numerical error in the CI formulas.

## 1. Garwood calibration does not cover the selected Omega=120 rate by branch

`job_cell` calibrates fixed64, fixed96 and stop20 cells before timestep selection. `job_tree` measures slope intervals but does not assess coverage of Garwood for the selected arm's rate. Selecting an arm using its observed events and exposure changes its distribution.

`selection_check.py` independently implements the synthetic constant-hazard model, the CP decision rule and the Garwood intervals. It uses 100,000 complete studies per configuration (64 runs per arm), in chunks of 10,000. Baseline run hazard is 0.0046, censoring duration 430, startup independent of subsequent waiting time; it does not import the repository calibration library.

For equal true rates in all three arms and startup probability 0.35 in every arm:

| Production branch | Studies | Garwood coverage of selected rate |
|---|---:|---:|
| tau after EQUIVALENT | 74,641 | 98.61% |
| tau/2 after INCONCLUSIVE | 21,436 | 90.82% |
| tau/2 after Part 1b EQUIVALENT | 1,808 | 84.85% |

The first two branches have Monte Carlo SE about 0.04 and 0.20 percentage points respectively; the third about 0.84 percentage points. These are not marginal threshold failures. They concern branch-conditional coverage, the standard already used for the slope. They are not an assertion that every unconditional whole-study coverage definition fails: including stopped studies or averaging over branches changes that question.

Consequently, the standalone minimum coverage 0.95025 cannot be presented as the demonstrated coverage of the selected Omega=120 rate within every production branch. Add selected-rate coverage to the full-tree diagnostics. If it fails the relevant standard, label that interval approximate/post-selection coverage not established, or use an independently sampled rate ensemble. Qualifying the existing interval does not require additional lattice runs. Keep the original standalone calibration result intact.

## 2. Part 1b is not an independent repetition of Part 1

`job_part1` checks independent arms at a fixed baseline hazard. `job_tree` reuses tau/2 after DIFFERENT but only generates second ratios 1 and 0.71, both within the equivalence margin. It does not measure the triggered Part 1b equivalence error at second ratios 0.5 and 2. Its arms also share one startup probability, whereas the primary comparison separately studies different startup probabilities.

An illustrative independent check has k(tau/2)/k(tau)=0.71 and k(tau/4)/k(tau/2)=0.5, with startup probability 0.35 in all arms:

- Part 1b is triggered in 25,219 of 100,000 studies.
- It declares EQUIVALENT in 2,285 of those studies: **9.06% conditional on being triggered** (Monte Carlo SE about 0.18 percentage points).
- The corresponding whole-study probability of triggering and then declaring equivalence is **2.285%**.

The conditional and whole-study probabilities must not be confused. This example does **not** demonstrate a whole-study false-equivalence rate above 5%. It demonstrates that the standalone comparison's 5% error interpretation cannot automatically be assigned to the triggered comparison. Define which error rate is claimed, simulate the triggered comparison at both margin boundaries, and describe the result before using it to claim validated equivalence.

No need to discard Part 1's completed calibration. The missing scope can be addressed as a supplemental diagnostic, preserving the original registered outcomes and treating any revised decision rule as an amendment.

## 3. Other qualifications

- The slope's 0.937 minimum is an empirical minimum among branch/configuration cells with enough replicates to decide eligibility. It is neither a confidence bound on coverage nor a guarantee for sparse branches or untested startup patterns.
- The full-tree slope model currently gives all timestep arms the same startup probability. Part 1's standalone power calculation allows (0.35, 0.05). State this scope restriction, or add that differential-startup case to supplemental tree checks.
- In the hypothetical 96-run setting, `job_tree` uses 96 runs for the tau/4 arm too. The registered Part 1b design specifies 64 tau/4 runs. The selected actual design is 64, so this mismatch does not affect its current branch simulation, but it limits claims about the hypothetical 96-run configurations used in eligibility.
- The recorded after-the-run implementation choices are procedural deviations. They can be judged on their merits and adopted prospectively for the driver; they should not be described as previously registered in full detail. The threshold correction with a separate recheck and unchanged original decisions is appropriately disclosed.
- Keep the stationary test descriptive when power is inadequate. Passing null-size calibration does not show that actual lattice hazards are constant.
- Replace “the slope stays unbiased” with “the mean slope estimate is close to its generating value” unless a bias estimate with Monte Carlo uncertainty is supplied.

## Recommendation

Retain the registered slope fallback and method decisions. Before presenting the production analysis as fully calibrated, add the selected-rate and triggered-comparison checks/qualifications above. This is an analysis-scope correction, not a reason to reopen MG3, tune calibration until it passes, or launch more lattice simulation.

Files: `selection_check.py` and `selection_check.json` contain the independent supplementary simulation and results. `inputs/` contains the reviewed calibration implementation, tests and decisions. Python 3.13.2, NumPy 2.2.4, SciPy 1.16.1 were used for the supplementary checks. The review did not modify the simulator or registered design.

## Reproduction

From the repository root:

```text
python -m pip install -r review/MG007_replication/requirements.txt
python review/MG007_replication/selection_check.py
```

This regenerates `selection_check.json` with six configurations and 100,000 studies per configuration. The independent selection check is self-contained and does not import the source snapshot under `inputs/`. That snapshot is retained to document the reviewed implementation; its tests also require the repository's MG3 dependencies. Generated caches are excluded from version control.
