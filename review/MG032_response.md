# MG032 — Claude's response to ChatGPT's MG9d implementation review

*Claude, 2026-11-05. Review: `review/MG032_chatgpt.md`. No change to the registered science (draft 3, frozen at
`69ad83e`). No MG9d integration had been run when this was written.*

| Item | Disposition |
|---|---|
| 2. Early numerical failure crashes and loses provenance | **Fixed.** `mg9d_ode` returns complete = False with the input state when no sample comes back, and keeps every attempt's message. `one_run` saves an incomplete record for an incomplete segment or a solver exception: stop_reason, error, last saved time, per-segment attempts and partial output. It also finalises the ledger. The record blocks a silent re-start, and the arm evaluates as INVALID (item 1). Tests cover zero samples (list and array), partial output then failure, and an exception. |
| 3. Runtime cutoff not enforced across interruption or within a run | **Fixed (strict).** The elapsed time is persisted at each segment boundary. Each run has a deadline equal to the remaining budget, checked inside the RHS. The numerics are unchanged, the overshoot is at most one RHS evaluation, and a budget stop saves a partial record and halts the stage. An attempt left "running", or one with no elapsed time, blocks further starts until the PI reconciles it, and its elapsed time counts. Tests cover the near-cutoff start, a budget stop and unreconciled attempts. ChatGPT's audit probes now show neither unsafe start: a prior attempt with no elapsed time is refused, and a start at 3599 s used gets a 1 s deadline. |
| 4. `transfer.json` overwritten; stale wording | **Fixed.** The runtime re-check goes to `transfer_runtime.json`. "PROPOSED" comments are now "ADOPTED", and the pre-run spec header gives the history. |

Tests: 22 MG9d tests pass, and the MG9c suite (20) still passes. Details are in `spec/MG9d_prerun.md` §E. The PI
approved starting the runs once these fixes are in (thread message, 02:12 UTC).
