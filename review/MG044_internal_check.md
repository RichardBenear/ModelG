# MG044 — MG10-C1 prospective internal check

2026-11-06. ChatGPT / Codex, lead author. No external review claimed.

PI instruction authorizes the next study recommended in MG043. The lead fixes a
small new ceiling (.5 production core-hour, .1 qualification/timing core-hour),
not an automatic extension of MG10's old budget. No external review was requested
with the execution instruction. The original result stays UNRESOLVED.

Six qualification tests pass: known detector fixtures; independent reduced drift
and dense serialization; full Euler/BDF adapter and ledgers; modal Jacobian against
full-lattice finite differences; decision boundaries; actual owned-process cap,
timeout, interruption and resume. The independent checker reproduces all 23 fixture
counts using periodic contour walks rather than the production prominence routine.
It verifies all six development timing files, hashes, initial/endpoints, dense
background coefficients, exact availability and full-field ledgers.

Failed development records are retained and charged. An initial coarse background
interpolation failed the unchanged 1e-10 check; exact adaptive dense polynomials
replaced it. An initial checker incorrectly demanded bit-exact initial linear
interpolation; it now uses the prespecified 1e-12 linear absolute tolerance. Neither
changes a scientific margin or a registered outcome; registration follows these fixes.

The six T5 timings cover every solver class; projection is 462.838314 worker-seconds
(0.1285662 core-hours), below 1350 allowed after reserve. Production is ten jobs:
two reduced references, two full 250-mode linear scans, six nonlinear trajectories.
The mode-selection rule is frozen before the scan; both full perturbation amplitudes
and the linear comparison window are frozen before nonlinear outcomes. No extra
mode or amplitude is substituted on failure. Source/spec hashes and fixed job IDs
are registered before production. The owned-process controller includes all future
prerequisite/dependent forecast costs in reserve checks.

Challenges and limits: exact homogeneous evolution and evolution with finite
disturbances are different questions. A large instantaneous eigenvalue is not a
finite-time gain; the registered calculation integrates the time-dependent response
instead. Only X-directed perturbations are scanned, not optimal singular vectors.
Early agreement cannot imply equality after nonlinear saturation. A prominent minimum
is numerical contrast, not a demonstrated particle. The threshold/merging rule is
retained only for this diagnostic, with a new fixed prominence floor; MG10 is not
rescored. Full raw responses and fields support checking rather than selecting
favorable snapshots. Closure and any broader physical claim remain PI decisions.

Same-author recomputation shares the original reference/smoothing and the earlier
independently written MG10 ledger helper, but not the new model/evaluator or the
prominence implementation. This is not independent external replication.
