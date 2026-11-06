# MG041 internal qualification — STOP

2026-11-06. Author: ChatGPT / Codex. This is not external review.

The PI authorized proceeding and explicitly waived external review for this batch.
The earliest isolated exact-case qualification exposes an error in our own MG040
candidate specification. Mandatory check 4 fails for four of the 21 finite Euler
availability paths. Stop is required by MG040 section 4 and RESEARCH_PROTOCOL.md.

Two unit tests pass: saved endpoints/discrete recurrence and expected first-order
refinement. A separate program uses 60-digit Decimal exponentials and integer powers
to recompute all 21 labels from raw saved checkpoints; all agree, including the four
failures. Shared dependencies are the raw arrays and NumPy reductions; no evaluator
or qualifying implementation function is imported by the recomputation.

No full MG10 dynamics, particle detector result, benchmark matrix, production
registration or scientific production has run. The remaining implementation and
qualification checks have not been claimed complete. The result is QUALIFICATION
FAILED, not BASELINE NOT RECOVERED and not a physical failure of the reservoir model.

One local worker was used; numerical-library thread variables were set to one for
the qualification invocation. Qualification recurrence took 3.49591410002904 worker
seconds. Tests took 0.255 seconds. Conservatively charging the entire first command
wall time (4.9153401 s) plus the entire recomputation command (0.8208497 s) gives
5.7361898 worker-seconds, 0.0015933861 core-hours, below the 900-second qualification
ceiling. No timeout, budget stop or production job occurred. These tool wall times
include startup/command overhead; they are not a forecast for full trajectories.

The evidence is retained in this persistent local checkout under
results/MG10/qualification and committed to Git; no sole temporary copy is involved.
Production's separate durable-store verification remains unattempted.

Recommended correction is in spec/MG10_qualification_correction.md, explicitly
proposed rather than silently adopted. It distinguishes the Euler discrete oracle
from the continuous ODE oracle and preserves scientific recovery thresholds.
