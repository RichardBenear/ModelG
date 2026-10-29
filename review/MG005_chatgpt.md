# MG005 independent review: MG2 result

Disposition: retain the registered **FAIL**. Independent re-scoring preserves C1/C5/C6/N1/N2 PASS,
C2 INCONCLUSIVE, C3/C4 FAIL. Recommend proceeding to MG3 by explicit PI decision after the measurement
bugs below are corrected and regression-tested; do not reselect b to turn this MG2 result into a pass.

## Scope and provenance

Downloaded 25 saved cell files, the deterministic reference NPZ, and reviewed source/spec/result
snapshots (40 files, 5,413,202 bytes). Each download was verified against its Git blob SHA;
`download_manifest.json` and `verified_inputs.json` record exact inputs. Repository main was read only.
Pickles were loaded with a whitelist permitting only NumPy array/dtype/scalar reconstruction.

The scoring formulas were implemented independently in `rescore.py`, rather than running Claude's
evaluator. This re-scores saved profiles, variances and event records; it cannot regenerate discarded
raw snapshots. Independently implemented Radau continuation in `independent_confirmation.py` checked
all 60 stored fields in the required strong-seed cells. One reported C4 replication was replayed with
the registered seed and timestep, then continued by the independent solver.

## 1. Independent scoring

| Criterion/control | Independent decision |
|---|---|
| C1, b=29 and b=29.16, Omega=1e8 | PASS |
| C2, both b | INCONCLUSIVE |
| C3 | FAIL |
| C4 | FAIL |
| C5 | PASS |
| C6 | PASS |
| N1 | PASS |
| N2 | PASS |

C1 upper bounds in percent (G/X/Y): b=29: .0700213/.0119820/.0122413;
b=29.16: .0654645/.0115064/.0117188. C6: .0686098/.00928088/.0102075.
N1 has no saved detections and no births; the tracking bug below therefore does not change its result.
N2 fails against the wrong reference for every field, so the discriminating-power control passes.

C5 strong-seed counts are 19/20,20/20,20/20; weak-seed count 0/20. Independent Radau continuations
agree with all 60 strong-seed saved field classifications (59 persistent, one excursion). Weak-seed
and N1 snapshots were not retained, so absence of detections is verified from saved summaries only.
These are the registered empirical-frequency decisions, not 95% guarantees that the population
probabilities exceed/fall below their thresholds.

### Correction to the reported C4 fraction

At b=29.16, Omega=300 there are **12 replication events in 11 of 24 attempted runs**, not replication
in 12 of 24 runs. At Omega=100 all 24 runs have at least one violating event (121 total). No recorded
violations occur at Omega>=1000. The correction changes the descriptive event fraction, not C4 FAIL.

## 2. C2

All usable ladder cells have bias lower bound zero. No cells are selected by the registered rule,
so alpha is undefined and INCONCLUSIVE is exactly correct. Insufficient power to resolve the bias
is supported; neither alpha=1 nor an alternative exponent is established. A zero truncated point
estimate is not proof of zero bias. C1/C6 establish high-Omega profile agreement, not the unmeasured
asymptotic exponent.

## 3. C3 interpretation

For b=29, Omega=100, first 12 attempted runs:
Omega*Var(X)=131.922129, run-level SE=1.625544, versus LNA target112.69 (+17.07%).
The C3 scoring failure is reproduced. Those 12 runs contain 16 births and 23 replications;
the cited 24 births/31 replications refer to all24 runs, not the C3 measurement sample.

The records do not establish that contamination dominates the discrepancy. Finite-Omega nonlinear
corrections to the leading-order homogeneous LNA are also expected; Omega*Var=C0+O(1/Omega)
in a regular weak-noise expansion. The shrinking discrepancy is consistent with that possibility.
Events, imperfect detection, selection by a changing far-field mask, finite-time effects and numerical
error can coexist. Saved farvar aggregates cannot disentangle them. A descriptive control would compare
homogeneous/event-free segments, mask-radius sensitivity, and finite-Omega covariance or structure
factors, without changing C3's registered scoring. Do not attribute causation from event counts alone.

## 4. Genuine structure production versus detector artifacts

Replaying b=29.16, Omega=300, replicate0 at the exact registered tau reproduces the replication event
at t=23. Independent Radau continuation from the t=43 field yields two persistent structures:
minimum phiX=-9.843503 and max|dphi/dt|=7.38253e-5 at continuation t=300, agreeing with the recorded
diagnostic. Recorded event time agrees exactly; this audit does not claim bit-identical complete
trajectories (the exact-array block-profile comparison did not pass across environments).

Together with the other two documented continuations, this supports genuine multiparticle production
in the numerical Model G system for the checked cases. It is not evidence about real-world particles,
does not validate all12 labels (or the crowded Omega100 labels), and does not uniquely distinguish a
split from nearby noise-driven nucleation. The amplitude and two-particle persistence are more
persuasive than merely reproducing a detector event list.

## 5. Design versus mean-field convergence

Finite-Omega transitions are compatible with the previously verified deterministic mean-map limit
and C1/C6. Thus C4 is not a counterexample to mean-field convergence. The deterministic parameter
distance b_r-b_mid is not a stochastic escape barrier, however: b is fixed during a run; noise moves
the state through configuration space. A narrow stochastic persistence region is a plausible
interpretation, not a mapped result. Choosing the deterministic window midpoint did not guarantee
noise-robust persistence, and a different b may still fail at Omega100. Do not assert that b_mid alone
is the established cause or that a b re-selection will rescue C4.

## 6. Concrete implementation defects to fix before MG3

### First homogeneous birth is omitted (P1)

`measure.track(...,x0=None)` assigns id0 to its first newly detected particle. Event generation then
treats every id0 as an original seeded particle and omits its birth. A regression input with one
detection at x=0 every half-time-unit from .5 to30 gives **no birth**, although persistence exceeds20.
If MG3 uses this event list, its first-event rate and startup-event exclusion will be wrong.

Only designate an original particle when x0 was supplied. Test homogeneous one/two-particle births,
seeded starts, burn-in boundaries and censoring. This is an implementation correction to the intended
definition, not a new scientific criterion. C5 bypasses this event path, and N1 has no detections,
so their recorded results do not change.

### Pooled audit cannot flag these data sizes (P1)

`evaluate_mg2.py` standardizes every12-run cell by its own sample mean/SD before pooling. For13
such cells (156 total), this bounds the largest possible subsequent leave-one-out statistic by
3.429853. The applied critical value is4.648985. Thus the pooled test cannot flag ANY dataset of this
shape, even one with an arbitrarily extreme within-cell outlier. A synthetic regression confirms it.
Its 'ok' result is consequently not additional evidence for A1.

Fix and validate the audit rather than treating this as accepted residual tail risk. An illustrative
within-cell external-studentization audit with Bonferroni over all1404 run/direction/field tests flags
none of the saved data (max statistic4.21751), but this is descriptive and not a retroactive registered
replacement. Record the original audit defect and any corrected audit transparently. Per-cell C1/C6
diagnostics were independently rechecked and passed; the numerical FAIL remains unchanged.

### MG3 stopping-rule cleanup

The gate text still says 'down to Omega=1' despite the amended grid ending at50 and the crowding stop.
Specify the no-event/no-bracket/too-few-unscreened-cells outcomes on the actual allowed grid before
execution. Retain the registered crowding exclusion, hazard checks, deterministic confirmation, and
rare-event timestep validation V5.

## Recommendation

Choose option(a): preserve MG2 FAIL, record its successful high-Omega agreement and low-Omega
limitations separately, and authorize MG3 after the technical corrections/tests above. An exploratory
b/Omega persistence map could be useful later. It must not replace this registered failure with a
post-hoc success. No reaction-rule change is indicated by these results.

## Reproduction files and limitations

- `rescore.py` -> `independent_rescore.json`
- `independent_confirmation.py` -> `independent_c5_confirmation.json`
- `check_pipeline.py` -> `pipeline_checks.json`
- `replay_replication.py 0` -> `replication_replay.json`, confirmation field NPZ

An initial replay trial used1/1500 rather than the registered floating-point literal
6.666666666666667e-4. Its random trajectory diverged; it is preserved separately as
`replay_rounded_tau_trial.json` and associated NPZ files and is not used for the registered-stream
comparison. The corrected replay used the exact registered value. This distinguishes stochastic
path sensitivity from changes in ensemble dynamics.
