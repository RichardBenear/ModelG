# MG033 — ChatGPT review of MG9d results

Reviewed revision: `d3bcfd35ad5b85a88ac59a6cdfd1e4cd78f06d98`. Production code: `9407f59`. Registered design: MG9d Draft 3, frozen at `69ad83e`.

**Recommendation: accept MESH EFFECT SUPPORTED on D02, with C4 DECREASING and the finer-mesh comparison unavailable. No result-blocking error was found.** The raw data support the registered classification. They do not establish eventual pinning, a continuum limit, or stationarity of the original MG9c patch.

## Independent replication

The replication script imports no ModelG production code and performs no trajectory integrations. It independently reconstructs the starting fields, Fourier alignment, envelope crossings, amplitude-weighted centres, minima, reaction/diffusion residuals, and the rule branches reached by these data.

- All 77 snapshot hashes match the manifest and per-run records. Independently reconstructed starts agree within `8.9e-15`.
- All seven records contain the expected 201 diagnostic samples and 2001 tracker samples, ten completed integration segments, and one successful solver attempt per segment. No numerical retry or interrupted attempt is recorded.
- At the 77 saved field times, independently recovered alignment agrees within `9.46e-10`; centre and width measurements within `7.2e-15`; final residuals within `7.2e-14`. Fields are stored every 1000 time units: the intervening diagnostic and tracker values were checked as records, not independently reconstructed from unavailable fields.
- Recomputed all six validity conditions from the full records and checked the field-based quantities at saved times. All D02 arms are VALID. R4 fails only item 5.
- Independently reproduced all endpoint speed comparisons and component labels, including the single reflection sign adjustment and the ordered overall reading.
- The ledger contains seven finished attempts totaling `0.0842152947` recorded core-hours (about 303 seconds), comfortably within the registered budget.

| Arm | Independent full-window speed | Validity |
|---|---:|---|
| R0 baseline | −3.449735e-7 | VALID |
| R1 tight solver | −3.449732e-7 | VALID |
| R3 reflected, sign adjusted | −3.449735e-7 | VALID |
| R2a, shift +0.05 | −1.410306e-6 | VALID |
| R2b, shift +0.10 | +6.042454e-7 | VALID |
| R2c, shift +0.15 | +1.177605e-6 | VALID |
| R4, h = 0.1 | approximately −1.4e-11 | INVALID: shape change |

The tiny R4 speed depends slightly on alignment optimization details; it is far below the `1.5e-8` resolution threshold either way. Its displayed digits should not be interpreted as a resolved velocity.

## Registered conclusions

**C1 NUMERICS-AGREE:** baseline and tight-solver slopes differ by at most `3.28e-13`. **C3 SYMMETRIC:** reflection reverses the measured motion as expected. Neither check shows evidence that ordinary solver tolerance error causes the observed baseline drift.

**C2 ALIGNMENT-DEPENDENT:** every shifted arm exceeds the registered DIFFER threshold. A quarter-cell shift gives about 4.09 times the baseline leftward speed; half- and three-quarter-cell shifts reverse direction. This is not merely a boundary-midpoint artifact: weighted centres and matched interior minima corroborate the alignment displacement, with small registered shape changes and nearly constant widths.

**C4 DECREASING:** the baseline late/early speed ratio is `0.84546174`, a decrease in magnitude of about 15.45%; the tighter solver gives the same label. This describes the baseline, not every shifted arm. Mesh dependence and relaxation can coexist.

**C5 R4 INVALID / LOST-IDENTITY:** use the registered composite label, but explain that R4 actually fails the shape criterion, not the cell-count/width identity criterion. Its `S_ref(T) = 0.00568151` exceeds `0.001`; all 18 cells persist and the full-run width change is `0.12306`, below the identity limit. Thus R4 cannot provide registered evidence that refinement suppresses the same patch's drift.

The overall rule is therefore **MESH EFFECT SUPPORTED**, regardless of the unavailable finer-mesh comparison.

## Interpretation and limits

The strongest conclusion is that the D02 patch's translation depends on its position relative to the lattice. This makes a lattice-induced relaxation or pinning interpretation considerably more plausible. A specific pinning landscape and a preferred final position have not been measured: the phrases “slides toward a preferred position” and “toward the next node” in GATES should remain explicitly hypotheses. The nonlinearly deformed shifted states and the finite observation horizon also leave room for lattice-dependent transient dynamics rather than eventual arrest.

R4's final residual independently reproduces at `2.1465e-9`, and its centre is effectively unchanged over the scored window. This is interesting descriptive evidence for a nearly stationary, adjusted finer-grid state, but not a valid matched-profile refinement result. Do not move the validity reference time to 6000 after seeing these data or rescore R4.

An additional descriptive observation: R0 now ends at residual `7.19518e-6`, below MG9c's old `1e-5` threshold, while its motion remains resolved. This is consistent with slowing, but MG9d starts from a filtered state and is a separately registered experiment. It neither retroactively changes MG9c's classification nor establishes exact stationarity.

## Implementation follow-up and minor documentation correction

The MG032 fixes are present: empty-output failure handling, saved incomplete records, interrupted-ledger blocking, in-run deadline checks, and separate runtime transfer verification. All 22 MG9d tests pass locally. The registered proposal differs from its frozen version only in the registration-status header; the reviewed production driver/evaluator have no post-run source changes. The repository records a PI launch go after the fixes.

One wording correction remains in `review/MG032_response.md`: an RHS deadline check does **not** guarantee overshoot of at most one RHS evaluation. Linear algebra, Jacobian evaluation, diagnostics, or file writing can occur between checks. Describe it as cooperative deadline enforcement, or add an external hard timeout if a strict wall-time guarantee is needed in future. This does not affect these completed runs, which consumed only about five minutes and never approached the cutoff. No rerun is warranted.

## PI recommendation

Sign off MG9d with the registered outcome and its R4 limitation. The study accomplished its bounded purpose: it identifies strong lattice-position dependence and separates that evidence from solver sensitivity and a mere centre-definition artifact. Keep MG9c closed as previously scored. Proceed to designing MG9b on explicitly lattice-conditional terms, with its own deterministic event-confirmation checks; further investigation of asymptotic pinning can be optional, separately registered work rather than a prerequisite.

Replication: `review/MG033_replication/audit.py`, `audit.json`, and `README.md`. No production trajectories were rerun and no repository scientific source or result files were edited.
