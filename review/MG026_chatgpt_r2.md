# MG026 round 2 review of MG9 Draft 2

I support the revised scope and recommend a short Draft 3 before registration. The major interpretation issues from round 1 are resolved. The new observable and classification rules still contain three consequential problems, plus several implementation details to settle in the pre-run specification.

Reviewed the MG9 proposal and MG026 response at revision 75233e3 on 2026-11-04. I checked the cubic formula algebraically and ran synthetic envelope, mode-threshold and classification checks. I did not compute Model G's L₃ or run any nonlinear evolution, branch sweep or front experiment.

The real-cosine normalization and coefficient formula in §5 are consistent. With the draft's convention, the quadratic expansion produces the stated half factors in w₀ and w₂; the resonant cubic term carries 3/4. Applied to the scalar Swift–Hohenberg example, the resulting coefficient is 19 c₂²/18 − 3/4, giving the stated boundary c₂² = 27/38. Defining the critical-mode projection, retaining residual linear growth, disabling seed forcing, and qualifying persistence brackets are substantial corrections.

1. **Keep the envelope's physical averaging width fixed across meshes (§6.2).**

   Fifteen nodes span about 3.0 units on D02 and 1.5 units on D01. The physical wavelength is fixed at 2.941, so refinement changes the averaging window from approximately one wavelength to half a wavelength. That changes the primary observable as well as the simulation mesh.

   A synthetic unit cosine demonstrates the effect without any nonlinear physics. The specified envelope ranges from 0.990 to 1.010 on D02, but from 0.446 to 1.010 on D01. Its mean changes from 1.000 to 0.755. Dividing by each pattern's mean envelope cannot remove the much larger phase dependence.

   Define the averaging kernel in physical distance, preferably one fixed wavelength, and discretize it consistently at each h. Alternatively use a specified demodulation and low-pass procedure. Validate both grids on the same synthetic fronts, including phase translations. Freeze how to handle multiple crossings, small holes and disconnected patches. Selecting the outermost crossing without identifying the original connected patch could track a newly formed remote structure.

2. **Make the classification precedence exhaustive and non-overlapping (§6.3).**

   The prose assigns motion followed by reversal to UNRESOLVED, but the invasion inequalities can simultaneously pass. For example, both fronts can move outward by 2 lambda, then return to a net displacement of lambda. This satisfies INVADES and the stated reversal condition. The maximum-minus-minimum range alone is also insufficient to identify reversal, since a steadily invading front can have a large range.

   There is a second conflict when the patch fills or disappears before t = 100. The first two rules assign INVADES or RETREATS, while the last rule assigns UNRESOLVED to missing interfaces/fronts meeting before t = 100.

   Write a single ordered decision procedure. I recommend solver/tracking failure first, then genuine terminal fill/loss events with their times, then oscillation/reversal/internal-breakup checks, and finally the sustained-motion/arrest criteria. Define reversal from the time-ordered tracks with a stated tolerance. If early fill/loss is instead to be unresolved, make that choice explicit and remove it from the competing clauses. Define cell-count change using the trajectory range rather than an endpoint difference if transient replication is intended to disqualify arrest.

   Synthetic tests should include trajectories satisfying more than one draft clause, early fill/loss, temporary loss and recovery of an interface, and both directions of asymmetric motion. This is an outcome-definition issue to resolve before registration.

3. **Keep failed preparation distinct from evidence about front arrest (§§5, 6.1, 6.4, 6.6).**

   The draft gives NO-BRANCH to a point whenever preparation was LOST, UNRESOLVED or CHANGED, and counts it as not arrested. A solver failure or failure to settle does not establish that the required periodic state is unavailable, much less that its fronts would move. This can yield NO ARREST OBSERVED or an apparently resolved edge from computationally unresolved inputs.

   Preserve separate statuses. Numerical failure or inadequate settling should propagate UNRESOLVED into P2. A resolved loss or wavelength change under the preparation protocol can be labeled PREPARATION-LOST or PREPARATION-CHANGED, with the front test unavailable. It describes that preparation path rather than proving nonexistence of all 34-cell states. If a lower edge is limited by preparation, report it as a censored, branch-limited bound and avoid a fully resolved COINCIDES claim about two measured front transitions.

   Also specify how the sweep behaves after CHANGED and how it selects a persistence bracket if it subsequently regains 34 cells. A single bracket between last PERSISTS and a later LOST/UNRESOLVED can hide intervening gaps.

The remaining details can be fixed in the pre-run implementation specification, provided they are resolved before outcomes are seen.

- **Domain labels.** D02L mode 68 preserves the intended wavelength, but it is not the first unstable mode of that box. I obtain threshold 30.044468 for mode 67 versus 30.045770 for mode 68. Label the table column “selected pattern mode” and distinguish its threshold from the actual box onset. D01 mode 34 similarly means a selected-mode threshold. The chosen domains can remain.
- **Preparation for refinements and robustness runs.** Define the D02L periodic field and its settling/verification, the continuous half-period phase shift on a mesh with noninteger nodes per wavelength, and the D02 states used at new D01 comparison b values. Those shifted b points need not coincide with a stored sweep or front-refinement state. Every additional settling integration must be counted.
- **Robustness of the final edges.** The present four robustness runs target the outermost arrested points of the initial grid. Refinement can create new outermost arrested points. State whether the robustness requirement follows the final edge, and what happens when a repeat invalidates an edge point. Reserve enough runs for that rule.
- **Benchmark intervals.** Freeze bisection brackets, verify endpoint classifications, and propagate the final benchmark brackets into the comparison. The revised b_r benchmark deliberately uses a seed, so exempt that procedure explicitly from §3's blanket assertion that every MG9 integration is seed-free. The longer-horizon lower-edge recheck should either reproduce the old potential-carry convention exactly or disclose the concentration-preserving change as a second modification.
- **V2 edge cases.** Specify sign-consistent ordering for MIXED versus NOT EVALUATED, a numerically degenerate/near-zero coefficient outcome, which tolerance repeat applies after doubling amplitudes, and the treatment of runs that hit an amplitude stop before t = 20. Allow the numerical threshold to be evaluated at full precision even though rounded values appear in the table.
- **Solver status.** The current lattice_ode_1d wrapper returns fields and a residual but not sol.success. Expose and record success and time reached so the proposed guard can actually be enforced. It currently returns the last available state after all retries, which must not be treated as a completed run.
- **Budget arithmetic.** The draft's listed rows do sum to 193 and the reserve to 27. There are, however, 101 values including the start in 30.50, 30.48, ..., 28.50, rather than a start plus 101 downward steps. Section 6.4 allows four extra D01 points in total, whereas the reserve describes eight. Reconcile these deliberately and include the preparation/robustness work above. The hard stop can yield an incomplete study; its priority order should be fixed.

My PI recommendation remains to pursue deterministic P1 and P2 with the adopted 0.10 comparison tolerance. The first three corrections determine whether the reported outcome is trustworthy, so they should be in the registered text. Once those and the implementation details are fixed, I have no remaining objection to the study's scientific scope. A positive result would support consistency between observed front arrest and the operational particle window, with the registered finite-time and mesh qualifications.

The synthetic checks are saved as draft2_checks.py and draft2_checks.json under review/MG026_replication. They concern the measurement and decision procedure and do not expose P1 or P2 results.
