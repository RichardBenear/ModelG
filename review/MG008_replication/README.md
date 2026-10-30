# MG008 — independent review of MG3b production results

ChatGPT, 2026-10-30. Reviewed snapshot: `d7f02c71c16f3bc891ffebc2bd4caa32d7ac8990` of RichardBenear/ModelG.

## Assessment

The reported numerical results reproduce. MG3b is appropriately treated as an estimation study, with no pass/fail verdict and no retrospective change to MG3. The timestep comparison supports the registered factor-two equivalence criterion at Omega=120. The four-cell rate trend is consistent with an activated curve over Omega=100–130 at the selected timestep. It does not establish a universal activated mechanism or time-independent hazard.

No new lattice run is needed merely to replace the approximate slope interval with an eligible 95% interval. The existing fallback and the post-selection label for the reused Omega=120 rate should be retained. Two supplementary interval labels should be narrowed as described below.

## Independent checks

- Downloaded 290 files (6,737,018 bytes) at the pinned commit and verified every Git blob hash. These include 258 completed run records: 256 analysis runs and two timing runs.
- Independently reconstructed startup status, post-burn-in events, exposure and early-stop endpoints from the first-event times. Checked saved confirmation times/classes and equality of the first-event field with its confirmation-cache field. No bookkeeping mismatch was found.
- All 258 stream keys are distinct. The three production cells stop at the correct first qualifying batch end. The reused timestep arm has the registered fixed sample size.
- Independently computed rate intervals, the timestep rate-ratio intervals, the GLM slope and profile-likelihood intervals, landmark calculations, exposure-median splits and p-values without importing the repository evaluation or calibration functions.
- Independently simulated stationarity-test power (10,000 studies per cell and alternative), using inverse integrated hazards and a separate exposure-median calculation. Every cell remains uninformative under the registered power rule.
- Summed the ledger and reconstructed the original 12-hour Part 2 batch-submission checks. No duplicate run-ledger key was found; every batch would have been allowed under the original ceiling.
- All seven MG3b tests pass locally. The first test attempt encountered a sandbox restriction on the default temporary directory; rerunning with a temporary directory inside this replication folder passed. No source change was needed.

Independent deterministic continuation of **all 302 saved confirmation fields** agrees with every recorded class: **206 persistent, 96 excursion, zero disagreements**. The independent solver uses separately implemented equations, an analytic sparse Jacobian and Radau. Results are in `independent_confirmation.json`. This substantially extends MG006's field validation, which checked persistent fields and only two reconstructed excursions.

## Results reproduced

### Timestep decision

At Omega=120:

- tau: 44 events, exposure 12891.5, 9 startup events out of 64 runs;
- tau/2: 39 events, exposure 12017.5, 10 startup events out of 64 runs;
- rate ratio tau/2 over tau: **0.95082645**;
- 90% interval **[0.64563271, 1.39722527]**;
- 95% interval [0.60169273, 1.49746170].

The 90% interval lies within [0.5,2], so EQUIVALENT and production at tau follow the registered rule. Part 1b is not triggered; its previously identified selected-comparison limitation is not used in this result. This supports equivalence at the chosen margin and operating point, not timestep independence across Omega.

Startup Fisher p=1.0 reproduces. The new data do not reproduce MG3's observed startup imbalance; failure to detect a difference does not prove identical startup distributions.

### Production rates

| Omega | Runs | Post-burn-in events | Startup | Exposure | Rate per length-time |
|---|---:|---:|---:|---:|---:|
| 130 | 48 | 25 | 0 | 15461.0 | 1.61697174e-5 |
| 120 | 64 | 44 | 9 | 12891.5 | 3.41310166e-5 |
| 110 | 48 | 28 | 19 | 4464.0 | 6.27240143e-5 |
| 100 | 32 | 21 | 11 | 1865.5 | 1.12570356e-4 |

Length is 100. Garwood intervals match the evaluation. The Omega=120 interval remains approximate under the pre-production addendum's branch-conditional assessment. The 0.943 figure is an observed minimum over qualifying calibration configurations and is not a coverage guarantee.

The production-cell cumulative event totals by batch are:

- Omega130: 8,17,25; stop after 48 runs;
- Omega110: 10,19,28; stop after 48 runs;
- Omega100: 12,21; stop after 32 runs.

No capped-cell exclusion or retrospective sample trimming occurs.

### Slope

The full four-cell fit gives:

- S = **0.0642090871**;
- profile interval **[0.0465150794, 0.0816773666]**;
- Wald interval [0.0466556422, 0.0817625320];
- deviance 0.22729587 on two degrees of freedom, nominal p=0.89257215.

The registered sensitivity without Omega100 gives S=0.0676080810, profile interval [0.0410654202,0.0942397181].

The fitted trend roughly halves the rate for each increase of 10 in Omega over this range. The low deviance means these four points fit the curve well; four points provide limited ability to identify the underlying physical mechanism. Agreement with MG3's slope is descriptive, not proof that two slope parameters are equal.

### Stationarity and landmarks

The exposure-median statistics reproduce:

| Omega | Split time | Early / late events | p |
|---|---:|---:|---:|
| 130 | 228.3205 | 9 / 16 | 0.229523 |
| 120 | 194.5208 | 19 / 25 | 0.451381 |
| 110 | 137.8636 | 7 / 21 | 0.012541 |
| 100 | 103.5577 | 8 / 13 | 0.383310 |

Both windows have equal aggregate exposure by construction, not equal numbers of runs. The registered alpha is 0.01, so none rejects. Our independent power estimates for x2 / x0.5 changes are respectively:

- Omega130: 0.104 / 0.087;
- Omega120: 0.273 / 0.156;
- Omega110: 0.093 / 0.045;
- Omega100: 0.072 / 0.036.

These support the report's uninformative labels. Differences from its 2000-replicate estimates reflect separate random streams and Monte Carlo error. The simulation checks the same fitted-parameter, per-cell diagnostic procedure; it does not establish power under every possible lattice alternative or post-selection condition.

The landmark risk sets and rate estimates also reproduce. Rates tend to rise at later landmarks in all four cells; Omega110's p is near the registered threshold. This deserves acknowledgment as an unresolved pattern, not a formal rejection or a new gate failure. Constant hazard remains a model assumption for interpreting the fitted rates and slope.

## Reporting corrections

1. **Do not transfer the four-cell coverage figure to the three-cell sensitivity interval.** `evaluate_mg3b.py` assigns the same `SLOPE_LABEL` to both fits, including minimum coverage 0.937. The calibration configurations used four cells, not the three-cell sensitivity fit. Label the latter an approximate nominal profile-likelihood sensitivity interval with coverage not separately calibrated. This does not change its numerical endpoints.
2. **Landmark and early/late interval claims need their own scope note.** Part 3 calibrated the principal designs at burn-in 50 and ratios between independent timestep arms. The landmark intervals use survivor sets at 100/150 after sample sizes were determined using the original event counts; the early/late ratio uses dependent contributions from the same trajectories. These are not the identical sampling designs. Retain their point estimates and nominal intervals as descriptive sensitivities; do not imply demonstrated 95% coverage from the main calibration. The reused Omega120 arm adds its existing selection qualification.
3. Use “nominal 95% interval, approximate” where the registration explicitly withholds a calibrated 95% claim. The report's existing approximate labels communicate the main limitation; clarify that a chosen chi-square cutoff does not restore the failed coverage requirement.
4. Startup fractions are not inserted as a numerical correction factor into k-hat, but excluding startup runs defines the survivor-conditioned population and its precision. “They enter no estimate” should not be read as meaning the conditioning is irrelevant.

These are scope/wording corrections, not a reason to change the registered decision or rerun MG3b.

## Budget

The ledger sum is **26.75280556 hours**: calibration 0.74544444, timing 0.15063889, pilots 0.00894444, tau arm 5.02491667, tau/2 arm 10.88108333, and Part 2 9.94177778.

The production total excluding calibration is 26.00736111. These reproduce the evaluation. We also replayed each Part 2 submission using the registered prior/mean-cost projection with the 1.25 safety factor and the original 12-hour ceiling. Every batch would be allowed; the largest used-plus-projected total is about 10.78 hours. Thus the result that the ceiling removal was not needed is supported by the actual submission rule, not merely by the final total falling below 12.

The after-Part-1 PI budget decision and post-data evaluator implementation remain disclosed. The checks found no change to event inclusion, event-count stopping or the estimator caused by the budget extension.

## Reproduction

All files are under `review/MG008_replication/`. `inputs/` is the pinned, hash-verified source/data snapshot. Reviewer code and outputs:

- `download_inputs.py`, `download_manifest.json`, `verified_inputs.json`;
- `rescore.py` -> `independent_rescore.json`;
- `power_check.py` -> `independent_power.json`;
- `independent_confirmation.py` -> `independent_confirmation.json`.

The independent continuation code separately implements the periodic lattice equations and analytic sparse Jacobian, uses Radau with rtol=1e-9, atol=1e-11, max_step=0.5, and applies the operational deep-minimum persistence rule over continuation times 200–300. It tests field-level persistence, not attribution to a particular detection or infinite-time basin membership. No stochastic production ensemble was regenerated.

Tested with Python 3.13.2, NumPy 2.2.4, SciPy 1.16.1, Numba 0.68.0 on Windows. No remote repository changes were made as part of this review.

From the repository root:

```text
python -m pip install -r review/MG008_replication/requirements.txt
python review/MG008_replication/rescore.py
python review/MG008_replication/power_check.py
python review/MG008_replication/independent_confirmation.py
python review/MG008_replication/verify_results.py
```

`verify_results.py` compares the independent arithmetic to the saved evaluation and rechecks input hashes. The deterministic continuations take several minutes. `download_inputs.py` can reconstruct the pinned inputs from the public repository if needed.
