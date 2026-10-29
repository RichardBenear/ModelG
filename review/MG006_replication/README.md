# MG006 independent review — MG3 results

Reviewer: ChatGPT. Date: 2026-10-29. Repository snapshot: `ed025ee4d3ccc6ac5b1023db043caf5299b79e76`.
The pre-run implementation specification was also checked at commit `5f8a7b55be15dde2375eeefe67090b3dc1ff1210`.
No repository source, criterion, or remote file was changed.

## Recommendation

Preserve MG3's recorded **FAIL on C2(a)**. The numerical criterion results reproduce. Spontaneous production of deterministically persistent particle-depth structures is supported by independent continuation, despite the gate failure. Proceed with a separately registered follow-up study; do not lower the C2(a) threshold to turn this result into a pass.

The report needs narrower interpretation of timestep convergence, stationarity, confidence intervals, and detection-level attribution. The compute overrun must remain explicit: cell-start checking does not enforce a hard cap. These qualifications do not erase the C2 failure or the observed persistent structures.

## Scope and verification

- Downloaded 224 files (2,986,202 bytes) at the pinned snapshot and verified each against its Git blob hash. `download_manifest.json` and `verified_inputs.json` record provenance.
- Independently re-scored all 168 committed completed run records, reconstructing startup/event status and exposure from first-event times. No mismatch in event, exposure, early-stop endpoint, or confirmation-time bookkeeping.
- Computed rates, registered intervals, C3 regression, C2 fractions and wavelengths without importing the repository's MG3 statistical or evaluation functions.
- Continued **all 36 Stage-B first-event fields** with a separately implemented periodic lattice ODE, sparse analytic Jacobian and Radau integrator (rtol 1e-9, atol 1e-11, maximum step 0.5). All 36 independently satisfy the registered persistent-structure rule at continuation times 200–300. This includes startup events.
- Replayed B_100 replicate 13 using the repository simulator and original streams. All three birth records reproduce, and the first-event field is array-identical. Independent continuation of the two reconstructed excursion fields relaxes to homogeneous (final X perturbations about 5e-26); the third field retains a persistent structure. The replay is not an independent stochastic simulator.
- Independently implemented the Euler mean-map N2 control: exact homogeneous start stays exact; the prescribed small random perturbation decays to 1.068e-11.
- All eight tracker tests and all fourteen MG3 tests pass, including the resume test. The pre-run text says fifteen MG3 tests, but this snapshot contains fourteen named test functions.
- No complete stochastic ensemble was regenerated. Excursion confirmation fields are not stored in the compact records: only the two replayed excursions received independent continuation. The remaining excursion labels were re-counted, not independently re-solved.

## 1. Re-scoring and exposure

For first confirmed detection time t, startup t < 50 contributes no exposure; 50 <= t <= 480 contributes t - 50; no qualifying event contributes 430. Stopping after the confirmation cannot remove exposure preceding the first event. Every completed run reaches the time needed for classification.

| Cell | Runs | Post-burn-in events | Startup | Exposure (run-time) | Rate per length-time |
|---|---:|---:|---:|---:|---:|
| B_150 | 32 | 5 | 0 | 12333.5 | 4.0539993e-6 |
| B_120 | 16 | 10 | 5 | 2158.5 | 4.6328469e-5 |
| B_100 | 16 | 11 | 5 | 694.5 | 1.5838733e-4 |
| V5_120 | 16 | 14 | 0 | 3423.5 | 4.0893822e-5 |
| C4_b29.16_120 | 16 | 6 | 10 | 394.5 | 1.5209125e-4 |
| C4_b29.3_120 | 32 | 3 | 29 | 11.0 | 2.7272727e-3 |
| N1_b28.5_100 | 16 | 0 | 0 | 6880 | 0 |

Length is 100. The published integer exposure column rounds half units; use the unrounded exposures to reproduce rates.

- **C1 PASS:** Stage A+B totals 5, 11, 15 at Omega 150, 120, 100.
- **C2 FAIL:** confirmation fractions 5/7, 15/21, 16/26. Same-Omega reference eligibility is 1/8, 0/8, 0/8, so the pre-stated replacement applies. Nucleated eligibility is 1, 3, 0. Wavelengths are 2.9413513 and 2.9599435 versus deterministic 3.0295767 (-2.9121%, -2.2984%); no eligible profile exists at Omega 100.
- **C3 PASS under the pre-run normal approximation:** slope -0.07257481, interval [-0.09357429, -0.05157532], weighted R-squared 0.99263285. The residual-based t interval is [-0.15201810, +0.00686849].
- **V5 PASS by its registered rule:** rate ratio 0.88269315, interval [0.36465263, 2.22117351].
- **C4 PASS by its registered calculation:** monotone point rates; high-b / b=29 ratio 58.86818, interval [10.41064, 228.62805].
- **N1:** zero first confirmed events and zero persistent-track outcome records; average detections per snapshot is nonzero, so this is not zero detector activity. **N2 PASS** independently.

## 2. Early stopping and the C2(a) estimand

The early stop is valid for a first-passage analysis: everything after the first qualifying event is irrelevant to its waiting time. It also leaves the specified profile window through detection + 50 available, except for the ordinary T=500 boundary.

For C2(a), the pre-run implementation explicitly chooses outcomes through and including the first persistent field. That is an intelligible **pre-first-event confirmation fraction**, including startup detections. It is not the fraction among all detections in a fixed 500-unit trajectory. The change in population comes primarily from ceasing confirmation at the first event, even though the simulation continues for another 50 units. Post-event detection quality could differ in either direction. Therefore:

- The reported FAIL is correct under the implementation specified before the run.
- It is not a full-horizon detector-quality estimate. Retain the explicit qualification wherever the 80% criterion is described.
- Neither stopping nor restricting the denominator should be retrospectively changed to rescue this gate. A new detector-quality study should choose its observation window and sampling unit separately from first-passage efficiency.

## 3. Confirmation versus transient structures

The continuation test is the accepted operational field-level rule. Independent checks support its implementation: 36/36 saved persistent fields agree, and the two replayed excursions genuinely decay. No confirmation-solver defect was found in those checks.

However, a field containing a persistent structure does not necessarily attribute that structure to the candidate detection. B_100 replicates 13 and 14 have `near=False` for their persistent outcome: no persistent endpoint falls within the diagnostic attribution radius. This is allowed by the field-level protocol and is **not grounds for retrospective rejection**. It does limit wording about every identified track crossing a separatrix. Endpoint displacement is not itself a proof of unrelated origin, either.

Prefer: “Some fields sampled after persistent detections relax to homogeneous under noise-free continuation; others retain particle-depth structures.” A detection that fails this finite-horizon test is an operational excursion, not a mathematical proof that its entire previous noisy trajectory never crossed a basin boundary.

Across these three cells, 18/54 sampled outcomes are excursions (33.3%); they are not mostly excursions. Any inherited wording equating C2(a) failure with “mostly excursions” should be corrected: below 80% confirmation permits a minority of excursions.

The reported profile-eligible fractions 1/5, 3/15, 0/16 are joint isolation, tracking and coverage qualifications, **not pure survival probabilities**. Loss, extra detections, identity changes and insufficient block coverage can each cause ineligibility. Label them profile eligibility rather than biological/physical survival without further analysis.

## 4. C3 inference and stationary hazard

The normal interval was present before the run. The t interval including zero does not automatically overturn that pre-stated decision: it estimates uncertainty differently. But the normal procedure does not have literally known variances. `1/n` is a plug-in large-count approximation for the variance of log rate; one cell has only five events, and cell extension/qualification depends on event counts. Three cells provide just one residual degree of freedom and little ability to distinguish an exponential law from another smooth decreasing curve.

Report “PASS under the specified normal approximation; consistent with activated scaling over Omega 100–150.” Do not describe this as establishing a large-deviation mechanism or a universal activated law.

At Omega 150, four early events and one late event give the reported p=0.068306. The early rate is about 6.71 times the late rate, but uncertainty is large. Not rejecting at alpha=0.01 with five events is weak evidence for constant hazard.

More decisively, B_100 and C4_b29.3_120 have **zero exposure in the late window**. Their p=1 is a degenerate comparison, not evidence supporting stationarity. The tau/2 cell has p=0.05990 (4 early, 10 late), also not a convincing validation. Retain the registered eligibility outcome but explicitly mark uninformative checks.

### Statistical “exactness” qualification

The registered Garwood and conditional-binomial numerical formulas reproduce. Their claimed exact Poisson/binomial sampling justification does not automatically transfer to these stopped first-passage data: exposures end at events and are random, and some cell sizes depend on event counts. Under a constant exponential hazard the likelihood is proportional to hazard^n exp(-hazard * exposure), so the point MLE is appropriate; identical likelihood form does not make the event count conditionally Poisson given its realized exposure.

A simple counterexample is one subject followed until its first event without censoring: its count is always one, its exposure is exponential, and conditioning on exposure cannot yield a Poisson count. Our illustrative calculation gives Garwood coverage about 97.12% in this example. That is conservative; **it does not demonstrate undercoverage in MG3**. It does show that the Poisson derivation cannot simply be asserted for this sampling design. Similarly, the conditional-binomial rate-ratio and window tests need survival-design justification or calibration to warrant an exactness claim here.

For this result preserve the registered calculations and qualify their model-based uncertainty. For a follow-up, use the individual censored first-event likelihood and calibrate the complete procedure, including any event-dependent extension. No new A1/block-independence assumption is needed to write that likelihood; time-homogeneous hazard is the substantive assumption already at issue.

## 5. C4 and severe left truncation

The three b=29.3 survivors have post-50 waits **5, 4 and 2**, totaling 11. Under a shared memoryless post-burn-in hazard, excluding startup runs and using this risk set is a valid conditional likelihood analysis. It is not a bookkeeping error, and 29 startup events need not be added to the registered post-burn-in rate.

Nevertheless only three survivors inform that rate; none provides late-window exposure. A startup exclusion rule does not establish a quasi-stationary survivor ensemble or absence of run heterogeneity. The numerical factor 58.9 is therefore an imprecise, model-dependent contrast, not a well-measured physical factor. Keep the registered PASS but report the risk-set size and waits beside it. More starts at this b, a separate shorter-window startup analysis, or a b with less severe truncation would be more informative for follow-up.

## 6. Stage-A design and selection

For an exponential first-event time with run hazard lambda from time zero,

    P(50 <= T <= 150) = exp(-50 lambda) - exp(-150 lambda).

It is maximized at lambda=log(3)/100 and equals 2/(3 sqrt(3))=0.38490018. For independent runs the maximum 4/4 probability is **0.02194787 per cell**. Thus the stop is poorly powered, not literally unreachable. Even twelve independent cells all at the optimal probability would have about a 23.4% chance of at least one 4/4 result. The 0.385 bound assumes constant hazard from the initial state; a model assumed only after burn-in need not obey it.

Selection follows the pre-run rules: only Omega 100 is partial (2/4); add screened neighbor 120, then 150; exclude 80 and insert no midpoint. The 120 startup and t=330 events do not create a [50,150] partial cell. Rates use independent Stage-B streams, not the Stage-A runs used for selection.

The crowding mean 2.257576 at Omega 80 can be checked against the excluded-cell JSON, but its four checkpoint trajectories are not committed. Thus I cannot reconstruct that screen independently from raw detections. This is a reproducibility limitation, not evidence that the reported exclusion is wrong.

## 7. V5 and startup asymmetry

The Fisher calculation for startup 5/16 versus 0/16 reproduces: two-sided p=0.04338154. Different streams and small samples can produce a chance difference; finite-tau changes to rare-event paths or startup relaxation are also plausible. The current data do not identify a cause. This exploratory p-value does not retrospectively alter the post-burn-in V5 criterion.

However “insensitive to halving tau” is stronger than this experiment supports. V5's interval admits substantial rate changes (approximately 0.36–2.22 times), and the startup distributions disagree enough to warrant follow-up. Prefer “V5 passed its registered post-burn-in comparison; startup sensitivity remains unresolved.” The exact-Euler deterministic limit does not guarantee finite-Omega rare-event convergence.

## 8. Budget and follow-up

The pre-run specification raises the cap to 18 hours but does not explicitly replace the hard-budget rule with permission to overrun an entire cell. The implementation checks the cap only at cell start; a cell can overrun without bound relative to the cap, and the evaluator detects only a `[BUDGET]` log marker. Consequently the absence of that marker is not proof that the cap was obeyed.

The cell-start convention was in the driver before the data, so this was not an outcome-dependent, post-hoc decision to continue. But the report acknowledges 19.03 hours: reconcile the hard-cap wording with that cell-start convention, or record the overrun as an exception. It does not turn the already observed C2 failure into a pass or remove the prior results. Under a literal hard-cap reading, the budget outcome could instead be INCONCLUSIVE with the observed C2 failure retained. Do not silently equate these conventions. My recommendation is to preserve FAIL with the overrun explicitly documented, rather than silently re-score the gate on an ambiguous resource rule.

Committed completed runs, references and pilots sum to **18.95763 hours**. The small difference from 19.03 is consistent with additional work in excluded/partial runs, but those wall-time records are not committed, so I cannot independently reconcile the full total. Publish a compact budget ledger including excluded work.

Recommended separately registered follow-up:

1. Keep MG3's failure record. Distinguish existence of persistent structures, detector confirmation fraction, persistence after nucleation, and rate scaling as separate scientific questions.
2. Use a fixed grid around 120–160 with additional interior values (for example 130 and 140). Omega 100 is useful as a diagnostic boundary but currently gives no isolated profile-eligible nucleated runs. Replace the poorly powered 4/4 trigger; with current information a new adaptive scan may be unnecessary.
3. Predefine sample size or censoring from a precision/budget calculation, preferably aiming for substantially more than five events per fitted cell. Use censored first-event likelihood fits and calibration covering any remaining adaptive rules. Treat `1/n` as estimated information rather than known variance.
4. Investigate burn-in and timestep together: report complete startup survival curves; predefine sensitivity landmarks such as 50/100/150 using the same first-event histories; compare tau and tau/2 with adequate independent samples, adding tau/4 if disagreement persists. A common initial-state ensemble can separate initial-condition variation from subsequent dynamics, but matching initial states alone does not couple tau-leap paths.
5. Record informative risk-set sizes and exposures for hazard checks. No late exposure means the late-hazard comparison is unavailable. Choose observation windows that can actually diagnose time dependence.
6. Study C2(a) separately with a clearly defined pre-first-event or fixed-window population. Save every confirmation field, including excursions and unresolved cases, plus compact tracker summaries and screen/budget records. Keep field-level confirmation distinct from attribution to a particular detection.
7. A revised confirmation threshold is a new scientific design decision, not a re-score. Decide whether the scientific target is detector precision or particle existence before choosing it; avoid picking a threshold solely because it would admit these fractions.

## Artifacts and reproduction

- `rescore.py` -> `independent_rescore.json` (raw-record arithmetic and all registered numerical statistics).
- `independent_confirmation.py` -> `independent_confirmation.json` (36 independent deterministic continuations).
- `replay.py` -> `replay.json`, `replay_birth_*.npz` (one stochastic replay and three independent classifications).
- `additional_checks.py` -> `additional_checks.json` (own Euler control, Stage-A bound, attribution and statistical diagnostics).
- `download_inputs.py`, manifests and `inputs/` pin the source material.

All reviewer code and outputs remain under `review/MG006_replication/`.

### Running the replication

Tested with Python 3.13.2 on Windows, NumPy 2.2.4, SciPy 1.16.1 and Numba 0.68.0 (`requirements.txt`). From the repository root:

```text
python -m pip install -r review/MG006_replication/requirements.txt
python review/MG006_replication/rescore.py
python review/MG006_replication/independent_confirmation.py
python review/MG006_replication/replay.py
python review/MG006_replication/additional_checks.py
```

The committed `inputs/` contains the hash-verified source/data snapshot; `download_inputs.py` can reconstruct it from the pinned public commit if needed. Generated Python and Numba caches are excluded. Floating-point roundoff and stochastic replay identity may vary across software versions or platforms.
