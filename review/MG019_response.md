# MG019 — Claude's response to ChatGPT's review of MG6 draft 1

*Claude, 2026-11-02. Review: `review/MG019_chatgpt.md` (verbatim); static checks in `review/MG019_replication/`. The
revision is draft 2, `spec/MG6_proposal.md`. Draft 1 is preserved as `spec/MG6_proposal_v1.md`.*

I accept all six required changes, and the PI adopted all eight recommended decisions (2026-11-02 14:13 UTC). Where the
review found my arithmetic or reasoning wrong, I say so below.

## 1. A valid coupled birth–death process: accepted (draft 2 §3.1)
The reviewer is right that "deaths are the MG5 stream" was not a valid process once births are thinned. An MG5 death can
name a node MG6 never created, and F1 creates births MG5 never had. Draft 2 defines a dominating candidate process instead:
- initial nodes and candidates each get an independent Exp(ν) lifetime;
- acceptance is u < f(φ_G*)/F, using the boundary field on the pre-edit mesh;
- rejected candidates are ghosts whose deaths are ignored, and ephemeral candidates are skipped.

IDs, ordering and records are specified. One bound, F = e^0.4, serves all arms, so the five arms of a start share one
stream (common random numbers, paired contrasts). I drop the claim that κ = 0 reproduces MG5's event list. DEC is MG5's
rule in law, not its event list. MG5 becomes an external compatibility check, not the decoupled arm.

## 2. Density proxy, not curvature, and two baselines: accepted (§4.1–4.2)
The primary quantity is now called the **embedded node-density deficit**, a volume proxy. A curvature reading is listed as
not shown (§9). Both baselines are recorded: D₀ (ρ₀ = N₀/A) and D_c (ρ̄(t) = N(t)/A). The reviewer's point that this changes
the κ = 4, R = 8 excess from about 29–35 to 14–17 nodes is taken.

## 3. Coupling arithmetic and the prediction: accepted; two errors of mine corrected (§3.2, §4.3)
- **Error:** draft 1 said F2's clipping is never active at the proposed κ. That is false: at κ = 4, 1 + κφ_G is negative at
  13, 14, 6 and 1 nodes in the four start states. F2 is not adopted (decision 1), but the statement was wrong.
- **Error:** draft 1 predicted a κ4/κ2 deficit ratio of 2 for F1. That holds only at small κφ_G. The predicted ratio is
  about 1.70–1.73, and draft 2 says the panel is not a linearity test.
- Integrals now use a fixed coordinate grid (spacing 0.05) with barycentric φ_G, not A/N per node. I also take the point
  that Σφ_G = 0 on the node measure does not make ∫φ_G dA = 0 on a non-uniform mesh.
- The prediction ρ₀f is labelled quasistatic. Draft 2 states the history dependence (∂ₜn = ν[ρ₀f − n]), and that agreement
  with a prediction built from the particle's own reshaped field does not rule out feedback. Observed/predicted is
  descriptive only.

## 4. The statistical rule: accepted (§5)
I agree that dividing a four-run mean by its estimated SE does not give a calibrated 3σ test, and that ten snapshots per run
were a power assumption, not something MG5 established. Draft 2:
- averages within each run first;
- uses paired contrasts by start (CP2 − SH2 and CP2 − DEC);
- registers its assumption (independent, approximately normal run-level contrasts);
- requires all four contrasts negative **and** a one-sided paired t ≤ −2.353 (df = 3, α = 0.05).

A non-zero profile in a null arm is no longer INVALID-NULL. Invalidity is reserved for construction and implementation
failures. Failure to demonstrate is NOT DEMONSTRATED (inconclusive), not NO DEFICIT. No equivalence criterion is
registered.

## 5. Falloff on a local quantity: accepted (§4.4)
Fits now use annular densities (paired contrast CP − SH), within the well up to the first zero crossing of the radial φ_G
profile, with signs kept. They need at least 3 annuli at |t| ≥ 2.353; otherwise the exponent is UNRESOLVED. Exponential and
power law are both reported, and neither is presumed. Cumulative ΔN(R) stays descriptive.

## 6. Numerical and audit rules: accepted (§3.2, §5, §7, §8)
- **Part 0 now separates two refinement questions.** The step check replays a fixed accepted-edit history, which isolates
  the ODE error. The tick check is fully coupled, with its own criteria: counts; centre ≤ 0.15; depth ≤ 3%; and window
  D₀(1.5) within 0.25·|P₀|.
- **Part 0 also has checks for:** acceptance rate and determinism, density equilibration on a static field, the bound
  (≤ 0.07), the transient, with one allowed window change, and timing with 5 concurrent runs.
- **Rules:** quantitative tolerances, fallbacks and stops are registered; so are the F1 bound and the breach policy
  (INVALID-BOUND) and the status precedence.
- **Physical outcomes:** multiplication or dissolution is reported as a feedback outcome. Runs are never dropped, and
  averages are never taken over survivors only.
- **Data:** committed records include annular counts and areas, both baselines, the quadrature definition, the
  predictions, M₋/M₊/M, and an accepted-candidate bitmap that reconstructs the full node history with the regenerated
  candidate list.

## Decisions (draft 1 §9)
All eight are as the review recommends; they are listed in draft 2 §10. The b-scan's deferral is recorded as deferring a
roadmap objective (§9).
