# Model G Hypergraph Study — standing instructions for Codex

This file is read automatically at the start of every Codex session in this repository. Current state and
next steps are in `HANDOFF.md`; read that next.

## What the project is
Mapping Paul LaViolette's Model G (Subquantum Kinetics) onto discrete, Wolfram-style token/hypergraph rewriting
rules, and testing step by step whether the rules reproduce Model G's particles (Turing-type dissipative solitons). Run
as a pre-registered research programme with formal gates. See `README.md` (charter), `InitialPlan.md`, `GATES.md`,
`STAGE_LOG.md`, and `spec/GLOSSARY_AND_ASSUMPTIONS.md` for every symbol and assumption (Ω, b, h, τ, λ*, etc.).

## Active roles (effective 2026-11-06, MG039)

- **PI: Richard Benear.** Scientific direction, modeling choices, registration, budgets, amendments and closure.
- **ChatGPT / Codex: research and implementation lead.** Planning, code, validation, authorized runs, analysis, recomputation and documentation. Not restricted to read-only review.
- **Claude / other reviewer: optional external review when requested.** No standing requirement to execute or approve each stage.

Read [RESEARCH_PROTOCOL.md](RESEARCH_PROTOCOL.md), the authoritative prospective role and review policy. Historical role descriptions in old specifications and reports do not override it.

## Rules of the protocol
1. Gate criteria, predictions and null controls are written into `GATES.md` and committed **before** the run that tests
   them. Anything a gate leaves open is written down as a pre-run implementation specification, also before the run.
2. A registered result is never re-scored after the data are seen. Deviations, design flaws and budget overruns are
   recorded explicitly, not smoothed over. A new idea becomes a new, separately registered study.
3. Amendments are versioned (`spec/phase2_amendment_N.md` → `spec/phase2_protocol_vX.Y.Z.md`). Earlier versions and
   records are preserved and marked superseded, never deleted.
4. Follow `RESEARCH_PROTOCOL.md`: internal checks and a separate recomputation path are required; external review is optional unless the PI explicitly requires it. Unresolved blocking issues stop production. Preserve historical reviews verbatim; same-author recomputation is not independent external replication.
5. Phase 2 protocol **v1.3.4** (`spec/phase2_protocol_v1.3.4.md`) governs its historical scope. Later gates have their own frozen specifications. MG10 first deterministic batch was registered as MG042 and completed as MG043: UNRESOLVED, PI closure pending. See results/MG10/MG10_results.md; further production needs its own registered scope. The registered MG10-C1 follow-up completed as MG045 and was PI CLOSED as MG046 with reference/amplification/measurement clauses supported; see results/MG10C1/MG10C1_results.md. Original MG10 remains UNRESOLVED.

## Working conventions
- Commit and push regularly with descriptive messages prefixed by the stage (e.g. `MG007: ...`). **Pull before you
  push.** The PI also commits from VS Code.
- Every number in a record should be reproducible from committed files. Commit raw run records (they are small), budget
  ledgers, and any checkpoint a decision depended on. `results/*/partial/` is gitignored scratch.
- Before stating a result, recompute it from raw files through a separate analysis path (not only through the evaluator), stating shared dependencies.
- Write for the PI in plain language: the conclusion first, then the evidence. Say plainly when something failed or
  when an earlier claim of yours was wrong.

## Compute
- Drivers are resumable at run level and use a process pool. Set the worker count with the environment variable
  `MODELG_WORKERS` (default 2). On the PI's 10-core machine use at most 8 workers, leaving capacity for background activity (MG039b). Set numerical-library threads to one per worker.
  - PowerShell: `$env:MODELG_WORKERS=8; python src/run_mg3.py`
  - cmd: `set MODELG_WORKERS=8 && python src/run_mg3.py`
- Cost model (measured on the 2-core cloud container, Linux): about 1.07 core-s per simulated time unit per run at
  τ = 6.667e-4, about 2.4 at τ/2, plus about 7 s per deterministic confirmation. Measure it once on the local machine
  before budgeting.
- Budgets are summed per-run wall time, which depends on the machine. For new studies register **run caps** and report
  core-hours alongside.
- Reproducibility: random streams come from `runner.seed_for` (numpy SeedSequence, entropy 20261027) and do not depend
  on the machine or the worker count. A Windows replay of an MG3 run reproduced its first-event field array-identically
  (`review/MG006_replication/`). Small floating-point differences across platforms remain possible; record the platform.
- Shared cache `results/MG003/ref_profiles.npz` (deterministic reference solitons). Compute a reference for a new b in
  the **main process** before starting the pool, so workers don't race to write it.
- Numba compiles on first call in each worker (a few seconds).

## Code map (src/)
- `lattice.py`: the tau-leap kernel `tau_leap_euler_1d`, lattice ODE, LNA, mean map. `rules_v1.1` in `spec/`.
- `measure.py`: smoothing, detection, tracking (`track`; A4 fix: only an `x0` track is the original particle), MG2 analysis.
- `confirm.py`: field-level deterministic confirmation (protocol §5.4).
- `runner.py`: `reference(b, h)`, `seed_for`, `one_run`.
- `stats_mg2.py`: calibrated bounds, Welch bound, A1 diagnostics, `audit_within_cell` (A4).
- `mg3.py`: segmented simulation with online detection and confirmation; first-passage statistics; Stage-B selection.
- Drivers and evaluators: `run_mg2.py`/`evaluate_mg2.py`, `run_mg3.py`/`evaluate_mg3.py`, `validations.py`.
- Tests: `python tests/test_track.py` (8) and `python tests/test_mg3.py` (14; one test takes about a minute).
