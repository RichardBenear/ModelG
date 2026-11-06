# MG10 first run plan — numerical definition and budget
## Draft 1, 2026-11-06 — MG040

**NOT REGISTERED. No MG10 trajectory has run.**
The PI approved the model direction and extended the capacity grid to include 100,000 and 1,000,000 in the project chat. This document finalizes a candidate deterministic plan; external-review disposition and registration are recorded separately before production.

## 1. Scope and units
Only MG10-A, deterministic MG10-B qualification, and MG10-C baseline recovery belong to this first batch.
Use P17, b=29, h=0.2, L=100, N=500, T_ref=T=300.
Initial active fields: (H) homogeneous steady state; (S) the committed h0.2_b29.0 reference soliton plus S0.
A/B availability starts at one; Z=W=0. FINITE has no replenishment.
Modes: CHEMOSTAT and FINITE at C_A=C_B in {1,10,100,1000,10000,100000,1000000}.
There are 16 mode/start pairs. Capacity changes stored values, not the lattice size or a loop over individual tokens.

## 2. Fixed primary measurement
Save all seven species (or equivalent scaled stocks plus cumulative transfers) at physical t=0,1,...,300.
Evaluate every integer checkpoint in [50,300], endpoints included. Use exact requested physical times, not rounded step indices.
For soliton starts, define s_s=RMS(phi_reference_s) for each of G/X/Y on all 500 nodes. The actual s values and reference SHA-256 are frozen by design/static_design.json.
For homogeneous starts, s_s=S0_s.
At common time t:
  e_s(t)=RMS(c_FINITE_s(t)-c_CHEMOSTAT_s(t))/s_s.
The primary error is max over all three species and all checkpoints 50..300. There is no shift, Nyquist filter or background subtraction in this primary measure.
Report absolute RMS and maximum absolute errors, spatial-mean offsets, and aligned/background-corrected differences descriptively.
Since all frozen s_s are positive, no denominator is estimated or substituted after outcomes.

CHEMOSTAT soliton must remain single at each checkpoint under the original fixed MG2 smoothing/detection threshold. FINITE must meet the same requirement to count as recovery. The homogeneous control and recovering finite arm must have zero particle detections at every scoring checkpoint.
A numerical failure, missing checkpoint or unavailable control is never labeled a physical particle loss.

## 3. Numerical methods and qualification
Solve each mode/start pair four ways:
- E0: synchronous expected-count Euler map, dt=1/1500 (exact physical units);
- E1: dt=1/3000;
- E2: dt=1/6000;
- O: an independently coded ODE RHS/reference, BDF with sparse Jacobian, rtol=1e-9, max_step=0.5, absolute tolerances 1e-11 for G/X/Y, 1e-12 for A/B availability ratios, and 1e-8 for Z/W concentrations and cumulative concentration-unit fluxes. Freeze the exact state ordering and Jacobian construction in the pre-run implementation note.
The Euler maps use section 4 of MG10_proposal.md, all from one start-of-step snapshot. The reference RHS must not call the production Euler increment.
dt here differs slightly from the historical rounded 6.667e-4; first prove the IDEAL Euler path agrees with the baseline mean map at these exact step values. This is disclosed and does not rescore MG8.
Use availability variables or an equivalent scaled representation to avoid subtracting small fluxes from enormous bath stocks during integration. Preserve the actual initial stock normalization in outputs.

Mandatory qualification before batch dispatch:
1. Exact integer stoichiometry and transfer-ledger unit cases.
2. Pointwise active CHEMOSTAT drift equals the baseline lattice ODE for homogeneous and nonuniform synthetic states.
3. The unchanged IDEAL stochastic regression uses only development seeds; no stochastic FINITE equivalence claim is made.
4. Deterministic A availability agrees with exp(-t/(C*T_ref)); scaled absolute error <=1e-8 at all saved checkpoints.
5. Independently reconstructed per-species changes from cumulative reaction and diffusion fluxes agree to <=1e-8*(1+absolute integrated throughput for that species), checked per site and globally. Do not normalize away errors by a huge initial bath stock.
6. Correct accounting of CHEMOSTAT external inputs and FINITE zero external input.
7. Numerical invalidity and missing data, mirror/check failures where applicable, checkpoint endpoints, interruption and budget states tested with synthetic records.

Small-system SSA and a finite-reactant tau-leap remain mandatory for any later stochastic FINITE study; they are not prerequisites for claiming a deterministic ODE comparison. No such sampler is yet qualified.

## 4. Decision rules
For each capacity/start, compute the primary error separately for E0/E1/E2/O.
Numerical agreement uses E2 versus O: normalize their per-species RMS differences by the same frozen s_s, over checkpoints 50..300.
- QUALIFICATION FAILED: mandatory pre-run check fails. Stop.
- NUMERICALLY UNRESOLVED: any required E2/O path fails, is incomplete/nonfinite, disagrees in detector count, or normalized E2/O error exceeds 0.001 (one tenth of the recovery margin).
- RECOVERS BASELINE for that capacity/start: numerically resolved, required detector conditions hold, and both O and E2 primary errors are <=0.01.
- BASELINE NOT RECOVERED at that capacity/start: numerically resolved and both primary errors are >0.01, or both agree on failure of the required FINITE detector condition while CHEMOSTAT remains valid.
- UNRESOLVED AT MARGIN: O/E2 lie on different sides of 0.01, or another outcome ambiguity remains.

The batch reports RECOVERS BASELINE ON THE TESTED RANGE only if the two largest prespecified capacities (100000 and 1000000) recover for BOTH starts and all mandatory controls pass.
If all required runs are numerically resolved but that clause fails, report BASELINE NOT RECOVERED ON THE TESTED RANGE. This does not deny a larger-capacity limit.
Otherwise report UNRESOLVED, with reasons. No data-selected capacity extension.
E0/E1/E2 ladders and error ratios are descriptive convergence evidence, not a claim of uniform convergence for arbitrary histories.
Missing administrative runs and budget stops retain their denominators and carry BUDGET-SHORT/INCOMPLETE.

The tight 1% soliton normalization intentionally sees G background offsets. The static calculation predicts C=10000 may miss it; this is why the PI approved the two larger capacities. It is not a measured result.

## 5. Dispatch, caps and compute location
Location: PI's local PC. At most eight worker processes; OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1, OMP_NUM_THREADS=1 and NUMEXPR_NUM_THREADS=1.
No affinity changes or global system settings are required. Respect a lower worker count if memory or foreground responsiveness requires it.
Qualification and two short timing cases use development keys and do not enter recovery outcomes.
- Qualification/timing cap: 0.25 core-hour.
- Fixed production matrix: 16 mode/start pairs x 4 solvers = 64 trajectories. No replacement trajectories for physical outcomes.
- Production ceiling: 4 core-hours; per-job timeout 600 seconds. No automatic extension.
- Dispatch only when the projection for the remaining fixed matrix fits the remaining ceiling with 25% reserve; otherwise stop and report.
- Projection sums costs by start/mode/method class. Use measured throughput and the three Euler step counts; do not multiply cost by reservoir capacity.
- Check the ledger before each batch. Charge all elapsed worker time, including interrupted/failed attempts. A complete record cannot be skipped without its cost being reconciled.
- Poll cumulative charged worker-seconds, including in-flight elapsed time, at least once per second. At the ceiling or a job timeout, stop dispatch and terminate only the affected owned jobs (all owned jobs for the global ceiling), retaining available checkpoints. Record termination and the bounded polling overshoot; do not silently let an entire batch finish beyond the ceiling. A timeout is not particle loss.
- Publish the measured projection before production. If it cannot fit this ceiling, return the concrete choices to the PI rather than silently changing the matrix or thresholds.

Planning estimate only: a few minutes per trajectory is a provisional allowance, not a measured benchmark. At eight workers, four core-hours corresponds nominally to about 30 minutes of fully occupied compute, plus serial/setup work; it is a ceiling, not a wall-time promise.
The exact solver and serialization can dominate; measure before making a firm forecast. Local compute does not eliminate assistant usage.

## 6. Retention and reproducibility
Keep full unit-time deterministic decision fields for every trajectory, source/parameter hashes, checkpoints, solver settings, flux checks, labels and ledger.
Conservative storage bound for 64 x 301 x 7 x 500 float64 field entries: 539392000 bytes (~514 MiB), before compression, flux records and metadata.
Local durable store: a named results directory under the PI's ModelG working area, verified before production; do not place the only copy in temporary cloud storage.
Commit a manifest and sufficient compressed decision arrays/recomputation inputs; verify and record the location of anything retained outside Git. Do not commit a single >100 MiB Git blob.
Two analysis paths must recompute errors, availability, detector outcomes and labels from raw checkpoints, with shared detector dependencies declared.

## 7. Review and authorization
Model direction, the extended grid and local eight-worker preference are PI-approved.
The PI determines external review per run/batch. There is no external-review decision recorded for this first scientific batch yet.
This remains a candidate plan until implementation qualification, measured projection, retention verification and prospective registration are recorded. Do not invent completed qualification or claim an approved measured runtime.
No MG10 simulation was run while preparing this plan.
