# MG040 internal static check
2026-11-06. Same-author check, not external review. No MG10 dynamics or benchmark was run.

- Capacity extension to 100000 and 1000000 explicitly approved by the PI before outcomes.
- Static A depletion is exact for the proposed deterministic P17 bath; its G-response estimate is homogeneous/quasistatic and not a rigorous bound on soliton errors.
- Primary soliton errors retain uniform background differences and normalize by frozen reference-pattern RMS; homogeneous errors normalize by S0. No outcome-dependent denominator, alignment, filtering or subtraction.
- The two largest capacities must recover for both starts; selection from favorable capacities is prohibited.
- The 16 mode/start pairs times four solvers give 64 deterministic trajectories. A firm runtime estimate requires a local benchmark; increasing capacity does not increase the number of lattice variables.
- Numerical reference agreement is stricter than the 1% recovery margin, and disagreement across the margin stays unresolved.
- Exact physical-time checkpoints include both decision endpoints. Missing data, numerical failures, budget stops and particle loss are distinct.
- The production cap charges active worker elapsed time; one-second polling and owned-job termination bound overshoot.
- Stochastic finite-bath qualification is explicitly deferred and no deterministic result may be presented as stochastic equivalence.
- Retention bound checked: 539392000 raw field bytes before additional records/compression. Actual durable location must be verified before production.
- External review remains the PI's per-run/batch decision; not recorded yet for the first scientific batch. Candidate is not registered and not a production go.

Next: implement the deterministic equations and qualification checks, benchmark locally, and record the measured forecast and retention/review disposition before registering production.
