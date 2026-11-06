# MG10 — Explicit A/B reservoirs and their limits
## Draft 1 — 2026-11-06 — MG039

**Status: DRAFT, NOT REGISTERED. No MG10 simulation, calibration or production run has been made.**
Author: ChatGPT / Codex, research and implementation lead. PI: Richard Benear.
This draft proposes a model extension; it does not alter rules_v1.1 or rescore any earlier gate.
The collaboration arrangement is [RESEARCH_PROTOCOL.md](../RESEARCH_PROTOCOL.md).

## 1. Question and boundary

Can a finite, explicitly accounted environment recover the existing fixed-reservoir Model G approximation over a stated horizon, and how does depletion change a particle's behavior?

The first objective is a validated reservoir construction, not an explanation of the physical vacuum. A/B stocks are explicit but still postulated. Geometry-derived supply and correlation-dependent supply are follow-ups. This stage does not establish quantum entanglement, physical energy, mass, gravity or a fundamental time tick.

MG8 showed that numerical scheduler changes can alter the homogeneous state and destroy the particle. MG10 must separate those effects from changes caused by finite supply. MG8's J stochastic losses at Omega=1000 after t=100 also rule out assuming indefinite control survival.

## 2. Proposed decisions and evidence

| Choice | Draft recommendation | Reason / limitation |
|---|---|---|
| Spatial setting | Existing periodic 1D lattice, h=0.2, L=100, P17, b=29 | Isolates supply from geometric changes; lattice-conditional |
| Bath layout | One well-mixed A/B compartment attached locally to each lattice node | Minimal explicit boundary; bath topology is not emergent space |
| Environment dynamics in first test | No autonomous bath hopping, recycling or imposed replenishment in FINITE mode | Depletion is explicit; none of these additional mechanisms is silently assumed |
| Time | Continuous-time Markov reaction process; J-style synchronous tau-leap only as an approximation | Fundamental-tick semantics are a separate undecided hypothesis |
| Primary comparison | Deterministic reservoir/field dynamics first | Separates finite capacity from counting noise |
| Stochastic component | Later bounded verification after exact and deterministic checks | Must account for finite-horizon J losses |
| Outputs | Explicit Z and W counters, W = LaViolette's sink species | Avoid confusion with Omega, the system-size parameter |
| External review | Optional; internal checks and separate recomputation required | Current PI-approved collaboration arrangement |

These are proposed modeling choices for PI review, not claimed consequences of SQK or StochasticVacuum.

## 3. State and normalization

At node i, active counts are n_G, n_X, n_Y. Its environment holds nonnegative integers n_A, n_B, n_Z, n_W.
For active species, concentration remains c_S=n_S/Omega, exactly as in earlier lattice work.
Z and W are inert output stores in this version, not sinks that erase tokens from the full ledger.

Let A_i^0>0 and B_i^0>0 denote the actual initial A and B counts, fixed for that run.
Define availability r_A=n_A/A_i^0 and r_B=n_B/B_i^0.
Keep the initial denominators fixed: renormalizing by the remaining stock would hide depletion.

Capacity is increased by increasing the well-mixed bath volume/initial stocks while preserving initial concentrations and the nominal flux through the fixed contact. The contact coefficient therefore scales inversely with initial stock. This is an explicit compartment-family assumption, not a derivation of fundamental chemistry.

Proposed dimensionless capacity coordinates, for reference horizon T_ref=300:
  A_i^0 = round(C_A * Omega * a * T_ref)
  B_i^0 = round(C_B * Omega * b * X0 * T_ref)
where X0=140/9, a=14 and b=29. Store the realized integers and coefficients.
C_A=C_B=1 provides roughly one reference horizon's homogeneous forward throughput initially; capacities are not themselves energies.
Vary C_A and C_B independently of Omega. Proposed deterministic screen: C_A=C_B in {1,10,100,1000,10000}, plus a chemostatted control. This grid and its costs must be finalized before registration, using static depletion estimates only.

For deterministic dynamics use A=n_A/Omega, etc.; A0=A_i^0/Omega and B0=B_i^0/Omega. Capacity and Omega limits are distinct.

## 4. Continuous-time reaction specification

All propensities below are evaluated in the current state. Invalid reactant configurations have propensity zero.

| ID | Full-state reaction | Propensity |
|---|---|---|
| R1e | A -> G | Omega*a*n_A/A_i^0 |
| R2 | G -> X | p*n_G |
| R2r | X -> G | g*n_X |
| R3e | B + X -> Y + Z | b*n_X*n_B/B_i^0 |
| R4 | 2X + Y -> 3X | n_X*(n_X-1)*n_Y/Omega^2 |
| R5e | X -> W | (1-g)*n_X |
| R1re | G -> A | (q-p)*n_G; zero for P17 |
| D_S | S at i -> S at neighboring j, S in {G,X,Y} | D_S*n_S/h^2 per neighboring direction |

B is a consumed reactant in R3e, not a catalyst. R4 requires two X reactants but changes net counts by +1 X and -1 Y; no permanent etheron identity assignment is assumed.

All rates are nonnegative for P17. No reverse Z/W reaction is added when the published reference combination is zero.
No cross-bath transport or geometric update belongs to this initial model.

### 4.1 Closed-mode accounting

FINITE mode conserves the total count across all seven tracked species and all sites:
  N_total=sum_i(n_A+n_B+n_G+n_X+n_Y+n_Z+n_W).
This follows from each reaction's stoichiometry and conservative active-species hopping.
It is a token-count identity, not an energy law.

Record per-channel event counts. They must reconstruct every population change, including bath consumption and output accumulation. Missing transfers and negative counts are implementation failures.
FINITE mode has no external injector or unrecorded output deletion.

### 4.2 Deterministic equations

With the existing lattice Laplacian:
  dG/dt = lap(G) - qG + gX + a*(A/A0)
  dX/dt = D_X*lap(X) + pG - (1+b*(B/B0))*X + X^2*Y
  dY/dt = D_Y*lap(Y) + b*(B/B0)*X - X^2*Y
  dA/dt = -a*(A/A0) + (q-p)*G
  dB/dt = -b*(B/B0)*X
  dZ/dt =  b*(B/B0)*X
  dW/dt = (1-g)*X.
At A/A0=B/B0=1, the active equations reduce exactly to the old same-lattice ODE.
For FINITE stocks these ratios change. Holding them fixed without an external ledger is forbidden.

## 5. Modes and numerical semantics

- IDEAL: unchanged old lattice kernel/reference. It is the regression control, not an explicit reservoir.
- CHEMOSTAT: A and B activities held at one. Every replacement A/B token is counted as external input; Z/W are retained or exports are explicitly counted. The active dynamics must agree with IDEAL at the mathematical-process level.
- FINITE: stocks evolve by section 4, no replenishment.
- Interventions (later sub-study): specific contact coefficients are changed while stocks and active fields are initially untouched. Every intervention is timestamped and accounted.

Use exact SSA only on small systems as the stochastic oracle. Full-lattice exact SSA is not assumed affordable.

If tau-leaping is used, all exit draws must read one step-start snapshot; arrivals/products wait until the next step.
Specify a bounded A-conversion draw and a reactant-safe R3e joint draw before stochastic implementation. In particular, independently drawing B consumption and X conversion, clipping excess draws, or allowing one B token to support multiple conversions would change the process.
A candidate construction may draw proposed X->Y exits then uniformly thin over-consumption against the available B stock, but this is NOT yet selected or validated. Its finite-tau bias must be bounded against exact small-system results and refinement.
The ideal kernel's per-token conflict handling and R4 catalyst-count convention remain baseline semantics; any new identity/reservation rule is a separate amendment.

Stochastic CHEMOSTAT need not reproduce an independently implemented sampler byte-for-byte unless its active sampling path is literally reused. Require exact regression only for IDEAL and mathematical/distributional agreement for equivalent new samplers.

## 6. Staged work and stop points

### MG10-A: formal model and static audit
Deliverables: full stoichiometry, equations, initial-state normalization, external ledger conventions, capacity scaling and limiting-case derivations.
Checks: nonnegativity, total-count invariant, A/B depletion estimate, output accounting and separation of capacity from Omega.
No dynamics are needed. Reject a candidate that conceals a source or cannot recover the intended active equations.

### MG10-B: implementation qualification
Deliverables proposed: src/mg10_reservoir.py, src/run_mg10.py, src/evaluate_mg10.py, a separate recomputation script, tests/test_mg10.py and spec/MG10_prerun.md.
Implement smallest closed examples first, then deterministic lattice dynamics, then the stochastic approximation if justified.
Required checks:
1. Exact integer accounting and state reconstruction from event counts.
2. Zero A removes R1e; zero B removes R3e; zero stock never produces negative counts.
3. CHEMOSTAT active drift matches the old ODE pointwise, including a nonuniform field.
4. IDEAL regresses against the old J kernel on shared seeds.
5. Small-system SSA validates finite-reservoir event distributions and the proposed tau approximation.
6. Time-step refinement separates solver bias from capacity effects.
7. Tests distinguish loss, track failure, numerical failure, incomplete runs and budget stops.
8. Separate analysis recomputes checkpoints, availability, conservation residuals and eventual labels.
Blocking failure: no scientific production; diagnose and record.

### MG10-C: recovery of Model G
Primary proposed experiment: deterministic T=300, homogeneous and reference-soliton starts, CHEMOSTAT and the common-capacity grid. Proposed grid gives 12 trajectories before refinement repeats; no stochastic sample-size claim yet.
Each finite run is compared with its same-solver/same-time control and a tightly solved ODE reference.
Record all-species active-field discrepancies, homogeneous background offsets, particle-centered discrepancies, detection count and bath availability.
Do not subtract a background offset from the primary error merely to obtain agreement; report background-corrected error separately.
Question: is the finite-bath construction within a prespecified tolerance at sufficiently large capacity over the declared horizon?
A decreasing error alone is not equivalence; failure to reach the tolerance is not proof of a nonzero limiting discrepancy.

### MG10-D: finite capacity and feedback
Estimation study, separately finalized after C: vary A-only and B-only capacities as well as both, measure availability, depletion, profiles, lifetime and environmental reaction history.
Compare like horizons. Do not assume Omega=1000 means no loss at T=500.
The environment can respond through depletion even without bath geometry. This is limited feedback, not a derived physical vacuum.

### MG10-E: causal supply interventions
Separately registered only after B/C pass. Compare cut/contact-off, sham and untouched states reached by a fixed preparation rule.
A-source cut, B-conversion suppression and a bulk exchange-link cut are different interventions: do not use their names interchangeably.
Turning off B can stabilize or destabilize a pattern; no universal "cut must kill" outcome is assumed.
If cutting supply changes the homogeneous background, use a background-matched control to distinguish global parameter drift from local dependence.
No intervention may be selected using future particle survival.

### MG10-F: recomputation, interpretation and PI closure
Publish all registered outcomes, failures, eligibility denominators, budgets and deviations. Independent external review remains optional.
State whether the model recovered the old approximation and what finite supply changed. Do not interpret a successful constructed reservoir as native vacuum emergence.

## 7. Proposed measurements and registration requirements

The following measurement definitions are the draft baseline:
- Effective rates a_i(t)=a*r_A and b_i(t)=b*r_B.
- Per-species raw relative L2 discrepancy against the same-time control; report absolute errors too, including cases with small denominators.
- Homogeneous/background offset separately from localized shape change.
- MG2 detector and smoothing, with threshold fixed from the committed original reference, not re-estimated from each finite bath.
- Initial-particle identity, persistence-qualified loss/replication, and TRACK-FAIL kept separate.
- Per-species stocks, output counts, channel event counts, external chemostat inputs and conservation/reconstruction residuals.
- Capacity error at common time; late-window drift/settledness only if estimable.

Candidate recovery tolerance: 1% for each active field over the stated comparison population, with a separate absolute tolerance where a reference norm is too small. This is a proposal, not a registered pass rule.
Before registration freeze the exact norm normalization, time grid, background window, horizon, numerical floor and allowed solver error; deterministic and stochastic outcomes must be separate.
If stochastic equivalence is attempted, specify margins, sampling unit, confidence construction and a power/coverage check on synthetic data before any production sample.
Do not use a nonsignificant difference as evidence of equivalence.

## 8. Budget, provenance and retention

No production budget is approved by this draft.
Before registration: static operation/memory estimate; run cap including refinements; measured timing on disclosed development seeds if separately authorized; core-hour ceiling; per-job timeout; batch-stop behavior and treatment of partial evidence.
Proposed first production scope is B/C only; D/E need their own frozen arms and budgets. No automatic progression to them.

Retain all decision inputs in results/MG10/: configuration, code/ref hashes, nominal checkpoints and physical times, integer stocks/event counts, deterministic fields at every comparison checkpoint, stochastic detection inputs and confirmation fields if used, complete ledgers and dispositions.
Large full trajectories need a verified durable retention location and manifest before starting. Do not repeat MG8's loss of full B2 arrays.
Store Python/NumPy/SciPy/Numba/compiler/OS versions as applicable, RNG specification and seeds.
Use seed_for(10,...) with written stream-key allocation before running. Development/calibration/production blocks must be disjoint; changing worker count cannot change assignments.

External baseline: RichardBenear/StochasticVacuum, commit 2161211ccc9d8427af98ebeefcabf7f379b83a37 (private archive).
Relevant source: materials/pc_rd1b_l0b_narrowed_local_kit/pc_rd1b_l0b_narrowed_local_kit/pc_rd1b_l0b_v35_seedguard.cpp and its analyses; tracker v4 and E2/E3 kits are methodological references.
Reuse only selected mechanisms with original file hashes and change notes. Preserve the failed D3 qualification and constructed-control interpretation; no historical PASS transfers automatically.

## 9. Outcome vocabulary

- QUALIFICATION FAILED: implementation or measurement apparatus failed; no physical conclusion.
- RECOVERS BASELINE within stated margins/horizon: only after all frozen equivalence clauses pass.
- BASELINE NOT RECOVERED on tested range: the tested approximation fails its frozen criteria, without claiming no larger-capacity limit exists.
- UNRESOLVED: inadequate precision, numerical resolution or range.
- BUDGET-SHORT / INCOMPLETE: explicit partial-evidence labels.
Finite-capacity behavior and interventions are separately reported; a lack of indefinite survival in a closed finite environment is not automatically a failed physical model.

## 10. Follow-up roadmap, not registered gates

MG11 candidate: local geometry controls explicit supply/contact rates; show an effective Model G regime rather than assign it by naming.
MG12 candidate: compare correlated and decorrelated classical environments with matched resources and contact opportunities.
MG13 candidate: only after an explicit quantum state/process is defined, distinguish entanglement from classical correlations.
Numbers and scope are provisional; only MG10 is reserved by this draft.

## 11. PI decisions before registration

1. Adopt or revise the local well-mixed, non-replenished compartment model and count bookkeeping.
2. Accept continuous-time dynamics for this first extension, keeping fundamental-tick alternatives separate.
3. Approve the final capacity grid, deterministic recovery margins and comparison horizon.
4. Approve numerical/SSA qualification, run caps, retention and compute budget.
5. Decide whether this stage warrants optional external review.

Current authorization covers the role/protocol update and preparation of this draft. It does not register MG10 or authorize scientific trajectories.

