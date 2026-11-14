# MG176 acquisition and independent noise-scale pilot

The new method passed 95/96 intact-core paths and 166/168 negative-control paths in this small conditional pilot. It gives the tracker a fixed period to acquire a core, then measures persistence after acquisition, while estimating noise strength from a separate background sample. The registered primary screen FAILED. This tests a more practical measurement procedure, but it still depends on a valid background reference, known noise structure and a supplied initial center. It does not establish that a Model G particle persisted, qualify the procedure for actual trajectories, or change the earlier failed results.

## Method and frozen scope

Source, criteria, seeds and input bindings were published at `4f9694f` before new synthetic draws. Seed 17600; one worker and numerical-library thread; 1200-second cap and 250 MB local reserve. There were 264 primary paths across two noise regimes, two temporal correlations and noise amplitudes 0.8, 1 and 1.25. Every path had 24 observations and its own 16 independent background-reference fields. Five separate sensitivity conditions added 10 paths. No Model G trajectories were generated.

The reference estimates one scalar covariance multiplier using all nonzero Fourier modes, with a 95% chi-square interval and eight frequency-band adequacy checks. The reference is independent of the tracking fields and the estimate stays fixed during tracking. Thus a damaged or extra core in a tracking field cannot directly inflate its own noise estimate. This protection depends on access to a valid external reference; the method does not infer a signal-free region from an unknown scene.

The tracker retains MG174 stationary/moving modes and structural alternatives. It plugs the estimated scale into likelihoods and uses a finite-reference F threshold for the innovation guard. That threshold has its stated predictive interpretation only at fixed true geometry under the Gaussian model; selecting geometry/history and using plug-in class weights do not produce calibrated posterior probabilities. The 0.95 structural cutoff and geometry spread limits were unchanged.

Acquisition requires two consecutive single labels ending by observation 5 (zero based). Positive paths count acquisition failure as failure, then require at least 90% single labels and 90% geometrically accurate single labels after acquisition. Absence must never produce a two-single acquisition anywhere in 24 observations. Other negative paths must acquire by 5, then avoid a single label at the damage deadline 15. Damage starts at 10. These endpoints differ from MG174; results are not a direct superiority comparison. Two paths per case must both pass the conditional screen.

## Primary results

| Noise regime | Temporal correlation | Noise amplitude | Intact paths successful / 8 | Negative paths successful / 14 | Reference adequate / 22 |
|---|---:|---:|---:|---:|---:|
| 0.1 | 0 | 0.8 | 8 | 14 | 22 |
| 0.1 | 0 | 1 | 8 | 14 | 22 |
| 0.1 | 0 | 1.25 | 8 | 14 | 22 |
| 0.1 | 0.95 | 0.8 | 8 | 14 | 22 |
| 0.1 | 0.95 | 1 | 8 | 14 | 22 |
| 0.1 | 0.95 | 1.25 | 8 | 14 | 22 |
| 1 | 0 | 0.8 | 8 | 14 | 22 |
| 1 | 0 | 1 | 8 | 14 | 22 |
| 1 | 0 | 1.25 | 8 | 13 | 22 |
| 1 | 0.95 | 0.8 | 8 | 14 | 22 |
| 1 | 0.95 | 1 | 8 | 13 | 21 |
| 1 | 0.95 | 1.25 | 7 | 14 | 22 |

Intact-path acquisition endpoints: {1: 75, 2: 17, 4: 2, None: 1, 3: 1} (None means acquisition failure). All 264 primary paths remain in the denominator. Reference adequacy passed 263/264. Estimated variance divided by injected variance ranged from 0.970412 to 1.03358; the stated reference intervals contained injected variance in 252/264 paths. This observed coverage is descriptive, not a population qualification.

## Separate assumption checks

| Sensitivity, two paths each | Successful paths | Reference adequate | Acquired by deadline |
|---|---:|---:|---:|
| weak_static | 1 | 2 | 2 |
| weak_moving | 2 | 2 | 2 |
| wrong_rho | 0 | 2 | 0 |
| shape_X | 2 | 2 | 2 |
| contaminated_reference | 2 | 2 | 2 |

Weak cores have amplitude 0.8 while the assumed amplitude remains 1. Wrong-rho fields have actual correlation 0.8 with assumed 0.95. Shape-X multiplies X-noise amplitude by 1.5 in reference and tracking fields. Contaminated-reference fields add a stationary core to each reference field while the tracking field contains its ordinary single core. These conditions are separate from the primary denominator and test different assumptions. A successful positive tracking endpoint does not mean contamination was detected; inspect the reference-adequacy column. A noise-unqualified output withholds a measurement and does not imply a particle disappeared.

Failures retained without threshold changes:
- regime 1, rho 0.95, noise 1.25, static, replicate 1: acquisition None, reference adequate True, emitted labels {'uncertain': 10, 'single': 14}.
- regime 1, rho 0.95, noise 1, split, replicate 0: acquisition None, reference adequate False, emitted labels {'noise_unqualified': 24}.
- regime 1, rho 0.95, noise 1, weak_static, replicate 0: acquisition 2, reference adequate True, emitted labels {'uncertain': 4, 'single': 20}.
- regime 1, rho 0.95, noise 1, wrong_rho, replicate 0: acquisition None, reference adequate True, emitted labels {'single': 1, 'uncertain': 23}.
- regime 1, rho 0.95, noise 1, wrong_rho, replicate 1: acquisition None, reference adequate True, emitted labels {'uncertain': 24}.
- regime 1, rho 0, noise 1.25, split, replicate 0: acquisition 2, reference adequate True, emitted labels {'uncertain': 10, 'single': 14}.

## Verification and storage

Separate recomputation passed 6576 label checks, 274 acquisition/path outcomes, 4384 real-space reference quadratics, 2192 band checks, and 548 selected-history quadratics. Maximum selected-history quadratic difference was 1.84433e-06. The driver checked 822 explicit full-prior endpoints. The checker shares templates/precision and the same author; it is not external replication.

Evaluation took 271.389 seconds; checking took 8.09103 seconds. Full fields and reference samples remain in `local_raw/MG176_acquisition` (84.5464 MB), with a public hash manifest. Only compact compressed decisions, source and summaries are published.

## What the failures mean

The three primary failures have different meanings. One high-noise stationary path missed the fixed acquisition deadline. One nominal split-control reference failed the band gate, withholding all measurements and therefore failing the acquisition requirement; this was a rejection of a valid generated reference, not a false particle-survival label. The other split path, at late independent noise and amplitude 1.25, incorrectly emitted single at observation 15: single weight 0.966513 versus multiple 0.0334867. Its estimated variance was 1.54919, with interval [1.51576, 1.58374] containing the injected 1.5625, so inaccurate scalar noise estimation is not necessary to explain that false-single decision.

The altered X-noise structure and contaminated references each passed the reference gate in both examples. Their successful tracking endpoints are therefore not adequacy successes. The current gate checks frequency-band energy summed over channels; it does not directly test species covariance or temporal correlation. The wrong-correlation tracking cases both failed despite acceptable reference gates, since those independent reference snapshots contain no information about tracking-time correlation. These results block general noise-adequacy qualification.

## Interpretation and next decision

The method estimates scalar noise strength, not an arbitrary covariance spectrum or temporal model. Its band test does not prove reference purity. It still assumes the template family and ordinary amplitude, known temporal correlation, supplied initial position, spatially homogeneous Gaussian noise and compressed geometry histories. It does not integrate scale uncertainty through all structural likelihoods. The conditional F safeguard alone cannot cure those assumptions.

Prioritize a saved-data audit of the false-single split path: distinguish insufficient secondary-core information from structural-prior and compressed-history effects at the correctly estimated noise scale. In parallel, specify reference checks that include species covariance and temporal behavior, and challenge them with contaminated samples. Then design a defensible background-reference selection test using held-out saved homogeneous Model G data. Do not expand production or relax thresholds before these gaps are addressed. Keep acquisition failure and false acquisition visible. Any revised estimator or scoring rule needs separate prospective validation. No automatic application to MG163 or new kinetic production follows from this pilot.

Administrative deviation: the GATES summary was added with the results rather than before execution. Complete criteria, source, seed, scope and budget were frozen and published prospectively in the specification and plan; no criterion changed after outcomes.

MG174 and MG163 remain unchanged. G is unresolved. Real-persistence readiness is false, and PI closure is pending.
