# MG10-C1 — Depleted homogeneous stability and resolved-minimum measurement

MG044, 2026-11-06. **REGISTERED prospectively (MG044). Only declared qualification/timing data preceded registration.**
PI instruction: “go ahead with the next study”. Scope is the follow-up recommended
in MG043. MG10's original UNRESOLVED result remains frozen. This is a new diagnostic
study, not a rescoring or a broader reservoir-capacity claim. No external review
has been requested for this follow-up; the lead performs the required internal
challenge and separate recomputation. PI closure remains a separate decision.

## Questions

1. Does a symmetry-preserving reduced ODE agree with uniform Euler evolution?
2. Does the depleted homogeneous trajectory amplify specified spatial disturbances,
   with the linear prediction supported by nonlinear runs before saturation?
3. Can a prespecified numerical contrast requirement reject false minima on flat
   depleted backgrounds while retaining known resolved motifs?

All physical parameters are MG10 P17, FINITE C=1, N=500, h=.2, T=300. No noise,
capacity extension, new reaction or energy/entanglement interpretation is introduced.

## Reduced reference and Fourier response

Reduced state is G,X,Y,uA,uB with the MG10 equations and no Laplacian; initially
(140/9,140/9,29/(140/9),0,0). Two independently expressed RHS implementations use
DOP853 and Radau, rtol=1e-11, atol=(1e-12,1e-12,1e-12,1e-14,1e-14), max_step=.1.
Save integer times 0..300 and dense solution coefficients for reproducibility.
Primary reduced-reference maximum normalized disagreement must be <=1e-8 using
scales (G0,X0,Y0,1,1). Exact A availability error must be <=1e-10.

Fourier modes m=1..250 have discrete Laplacian eigenvalue
-lambda_m=-4*sin(pi*m/500)^2/.2^2. uA perturbation is initially zero and decoupled.
For physical perturbation coordinates (dG,dX,dY,duB), the variational matrix is:

    [-1-lambda, .1, 0, 0]
    [1, -1-29*rB+2*X*Y-lambda, X^2, 29*X]
    [0, 29*rB-2*X*Y, -X^2-12*lambda, -29*X]
    [0, 29*rB/B0, 0, -29*X/B0]

Use normalized coordinates (dG/G0,dX/X0,dY/Y0,duB) and initial vector (0,1,0,0).
Integrate all 250 directions with BDF (primary) and Radau (reference), rtol=1e-9,
atol=1e-12,max_step=.2. The primary uses the DOP853 base; the reference uses the
independent Radau base. A block-sparse analytic Jacobian is checked against finite
differences of the original full-lattice drift on synthetic states before use.
This tests this X-directed disturbance family, not the full optimal singular gain.

Freeze the selection rule now: choose the smallest m attaining the largest
integer-checkpoint max |dX/X0| in the primary linear calculation. This planned
selection is based on the registered linear scan, not nonlinear outcomes. Report
all modes, not just the selected one. If missing/nonfinite, stop; never substitute
a hand-picked mode. Both references must agree in peak log gain to absolute .02
for modes whose primary peak gain is >=1; otherwise label numerical uncertainty.

## Fixed nonlinear matrix and measurements

Six N500 trajectories: homogeneous zero perturbation and two X-only cosine
perturbations, epsilon in {1e-6,1e-8}, each with E2 and O+.
Initial X_i=X0*(1+epsilon*cos(2*pi*m*i/500)), other fields unchanged; phase zero,
using the selected m, no RNG. Zero perturbation is structurally uniform initially.
E2 uses the unchanged MG10 synchronous Euler dt=1/6000. O+ uses the existing
independent MG10 ODE RHS/Jacobian with tighter BDF settings: rtol=1e-10,
atol=1e-12 active,1e-13 depleted fractions,1e-9 cumulative fluxes,max_step=.25.
All 301 states and cumulative ledgers are saved; accounting and A checks remain
the adopted MG10 tests. No output clipping, mode replacement or extra run.

Uniform E2 must agree with the primary reduced ODE within 1e-4 normalized active
RMS over all integer times. Its per-species spatial range must be exactly zero.
Zero-perturbation O+ departure is descriptive: reproducibility of roundoff onset
is not required, and departure is not evidence of spontaneous physical nucleation.

For seeded trajectories compute the complex cosine Fourier coefficient of X/X0,
subtracting its spatial mean first. Factor is 2/N, or 1/N at Nyquist. The prediction
is epsilon*v_X(t). The linear comparison window consists only of integer times
with 10<=|v_X(t)| and epsilon*max(|v_G|,|v_X|,|v_Y|)<=1e-3. Require >=3 such
checkpoints for each epsilon. At every eligible point, relative complex prediction
error must be <=.05 for BOTH E2 and O+. No shift or phase fit. Report later
saturation descriptively; no whole-trajectory equality claim after instability.

## Measurement qualification (a diagnostic, not a particle definition)

Keep the old detector unchanged. The new resolved-minimum detector uses its
reference threshold and merging distance, but requires periodic peak prominence
of smoothed -phi_X >= eta, where eta=1e-8*max(1,X0,|A_s|), fixed from the committed
reference. Center before FFT smoothing, then restore the mean. Find prominences
on three periodic copies and keep central-copy minima, then merge within lambda/2
as in MG2. This is a numerical contrast rule. It does not prove particle identity
and is not applied to rescore MG10's frozen recovery outcomes.

Qualification, before production: zero detections on every flat field at X in
{0,.1,1,5,10,X0,30}, including constant roundoff-sized ripples <=1e-12; retain one
committed reference soliton and one/two imposed Gaussian troughs on background 5,
depth1,width.8, centers 0 and25; translations/seam invariance; sinusoidal ripple
amplitude eta/100 rejected and eta*100 retained on background5 at mode16. Record
all fixtures and both old/new detector counts. Missing tests or failed required
conditions stop production. Known motifs are imposed validation, not discoveries.

## Outcomes and resources

Report separate clauses rather than invent a recovered MG10 verdict:
- HOMOGENEOUS REFERENCE RESOLVED if reduced agreement, uniform Euler checks and
  mandatory ledgers pass; otherwise UNRESOLVED or QUALIFICATION FAILED as applicable.
- CONTROLLED AMPLIFICATION SUPPORTED if peak selected gain >=100, linear numerical
  checks and both nonlinear epsilon-window checks pass. If complete and peak <100,
  NO >=100x AMPLIFICATION IN THE TESTED FAMILY/HORIZON. Otherwise UNRESOLVED.
- RESOLVED-MINIMUM QUALIFICATION PASS only if every prescribed fixture passes.

Local <=8 workers, one numerical-library thread each. The lead allocates a new,
conservative ceiling under this execution instruction: qualification/timing <=360
worker-seconds; production <=1800 worker-seconds (.5 core-hour), per-job600 seconds,
poll <=.5 seconds, owned jobs only. This is not an invented carryover approval of
the earlier four-hour cap. No automatic extension. Publish a measured conservative
projection and require it <=75% of the production ceiling before registration.
Development timing consists of six T5 cases covering both reduced solvers, both
250-mode linear solvers, and both seeded N500 solvers (development mode7,epsilon1e-5).
These do not enter scientific outcomes. Forecast each class from twice its measured
compute-only time times60, plus max(5s,twice maximum observed startup/serialization/
poll overhead); include all remaining prerequisite/dependent jobs in reserve checks.
Primary/reference reduced and linear jobs plus the six nonlinear jobs are the
fixed production denominator (10 jobs). Run prerequisites first; an unavailable
selected mode leaves all dependent nonlinear jobs INCOMPLETE, never replaced.

Reuse the already tested MG10 owned-process controller where possible; qualification
must cover its adapter and source/spec hashes. Retain all raw fields, reduced dense
coefficients, linear responses, fixture inputs, hashes and complete worker ledger
in persistent results/MG10C1/, with no blob >=100MiB, and commit/push the evidence.
Separate recomputation must read raw files without importing the evaluator or new
model implementation; shared original detector/reference/ledger helpers are declared.
No production begins until this specification, qualified implementation, measured
projection and registration are committed. No new physical particle-rate claim,
fundamental tick law, quantum model or finite stochastic sampler is in scope.

Implementation note before timing: retain the actual adaptive DOP853 degree-7 and Radau degree-3 dense polynomials in power-coefficient arrays with their interval knots. A coarse Hermite attempt failed the unchanged 1e-10 midpoint-accuracy check during qualification and was replaced; its failed test record is retained and charged. Test the serialized representation against the solver at interval midpoints. No threshold or scientific criterion was relaxed.
