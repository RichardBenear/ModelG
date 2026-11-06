# MG042 qualification and registration check

2026-11-06. ChatGPT / Codex, lead implementation author. No external review claimed.

**Proceed within the fixed budget.** The PI adopted the MG041 oracle correction
and waived external review. The prior MG041 failure is preserved, not rescored.
No scientific production preceded registration.

## Evidence

- Seven deterministic qualification tests pass: exact integer stoichiometry and
  ledger cases; IDEAL/CHEMOSTAT Euler agreement at all three exact steps; pointwise
  nonuniform drift and Jacobian directional derivatives; zero stock and mirror;
  finite and external accounting at C=1 and 1000000 using every solver; unchanged
  IDEAL stochastic regression with a development-only seed; invalid/missing/budget
  decision records. No FINITE stochastic implementation is claimed.
- Actual owned-process timeout, global-cap, interruption and run-level resume
  checks pass, with raw ledgers retained. These intentionally tiny test ceilings
  are distinct from the qualification/production ceilings. Observed poll overshoot
  is retained rather than erased or relabeled a physical failure.
- Twenty-one one-node, homogeneous, T300 development checks cover all capacities
  and Euler steps. An independently coded raw audit checks the adopted exact
  availability oracle and throughput-normalized ledger; all pass. They are not
  N500 particle outcomes or members of the scientific matrix.
- Two declared T10, N500 finite C1 soliton timing jobs (E2 and O) pass. Separate raw
  recomputation confirms their availability and seven-species accounting. Largest
  per-site scaled residuals are about 2.65e-14 and 1.04e-15 respectively.
- Conservatively charged qualification/timing cost is 39.172795 worker-seconds,
  below 900. Earlier failed development checks and their diagnostic are charged.
- Projection is 6488.182472 worker-seconds, 1.802273 core-hours, including 2x
  extrapolation and per-job overhead. This is below the 10800 seconds permitted
  after the 25% reserve and below the 600-second per-job forecast limit.
- Local persistent retention passes write/read and free-space checks. All raw
  decision fields and fluxes will be committed; no single blob reaches 100 MiB.

## Design challenge and limits

Reservoir availability is represented by depleted fractions to avoid high-capacity
stock cancellation. Independent accounting checks do not normalize by huge stocks.
ODE constant-coordinate elimination fixes a development failure at about 1e-30;
no negative output is clipped or accepted by a relaxed tolerance. The failed test
record is retained as full_tests_attempt1.json. The final suite passed afterward.

The ODE RHS is independent of Euler, and its sparse Jacobian is tested numerically.
Both share the scientific equations, so solver agreement alone cannot establish
physical correctness. Independent raw accounting uses separately written formulas.
The evaluator and checker share the frozen MG2 smoothing/detector and reference;
same-author recomputation is not external replication. Production checks still
have to pass; qualification does not predict recovery or particle survival.

The forecast uses the same conservative envelope for every class, based on the
two prescribed timing cases, rather than separately measuring every capacity or
inspecting preliminary recovery. Nonlinear stiffness could still vary in T300
paths. Timeouts, the global cap and complete denominators guard that uncertainty.
No extension, replacement run or post-outcome threshold change is authorized.

Only deterministic A/B/C are registered. The reservoirs are an explicit accounting
model, not a demonstrated physical vacuum, energy law or entanglement mechanism.
Gate closure remains the PI's decision after results and recomputation.
