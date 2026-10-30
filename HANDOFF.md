# Handoff — state of the project and next steps (2026-10-29)

Written by the Claude session that ran MG1–MG3 in the cloud, for the Claude Code session that continues on the PI's own
machine (10 cores). Everything referred to here is committed. Read `CLAUDE.md` first for the rules and conventions.

## Where things stand
| Gate | Question | Result |
|---|---|---|
| MG1 | Continuum PDE reproduces the published stationary soliton | **PASS by PI exception** (the registered FAIL on 1D stationarity is preserved) |
| MG2 | Token rules converge to the lattice mean field as Ω → ∞ | **FAIL as registered** (C3 at Ω = 100; C4 noise-driven replication at b_mid). The core convergence tests C1 and C6 PASS; C5 PASS |
| MG3 | Particles nucleate spontaneously from intrinsic noise | **FAIL as registered** (C2(a): 0.71, 0.71, 0.62 confirmed vs ≥ 0.80). C1, C3, C4, V5 PASS. Budget overrun 19.03/18 h recorded as an exception |
| MG3b | Follow-up: time step, rates with more events, calibrated uncertainty | Estimation study, no verdict. Part 1 EQUIVALENT (τ/2 vs τ ratio 0.95); S = 0.064 per unit Ω, approximate 95% 0.047–0.082 (MG3b record) |
| MG4–MG7 | Irregular graphs, evolving hypergraph, gravity-like deficit and mass, multiway | PENDING (not designed in detail yet) |

Key physics so far, on the h = 0.2 lattice at b = 29:
- The deterministic limit is reproduced to better than 0.1% (MG2 C1, C6).
- The single-particle existence window is narrow (b_s ≈ 28.93, b_r ≈ 29.40), and noise at Ω ≤ 300 drives loss near the
  lower edge and replication near the upper edge (MG2).
- Spontaneous nucleation exists. Its rate rises steeply as Ω falls: about 4 × 10⁻⁶, 5 × 10⁻⁵ and 1.6 × 10⁻⁴ per unit
  length and time at Ω = 150, 120, 100. This is consistent with activated scaling, with S ≈ 0.073 per unit Ω (MG3).
- About a third of persistent noisy detections are transient (they relax without noise).

Full records: `GATES.md` (each gate's "result" section), `STAGE_LOG.md`, `results/MG2/MG2_evaluation.md`,
`results/MG3/MG3_evaluation.md`. Reviews: `review/` (most recent: `MG006_replication/README.md`, ChatGPT's full MG3
review, and `MG006_response.md`).

## Immediate next steps (updated 2026-10-30, MG007d)
0. **MG3b: done and signed off 2026-10-30** (`results/MG3b/MG3b_evaluation.md`; review `review/MG008_*`). The PI
   decides what comes next (MG4 design, or deferred items D1–D3).
   **MG4 design: draft 4 written 2026-10-30 (MG009b, after two ChatGPT reviews; noise-free only, PI decision)**, `spec/MG4_proposal.md`; reviewer is comfortable with registration; awaiting PI approval.

### Earlier steps (2026-10-29, MG007)
1. **Local environment: done.** Python 3.13.2, numpy 2.2.4, scipy 1.16.1, numba 0.68.0; tests pass locally (8 and 14).
   Measured cost at Ω = 150, 9 workers: 0.89 core-s per time unit at τ, 2.43 at τ/2.
2. **MG3b follow-up: REGISTERED** in `GATES.md` (MG3b), design frozen as `spec/MG3b_followup_proposal.md` at commit
   `202c00f`. Reviews and responses: `review/MG006b_*`. Estimation study, no pass/fail; MG3's FAIL is unchanged.
3. **Part 3 calibration: done 2026-10-29 (MG007a)**; see `results/MG3b/calibration/README.md`. Awaiting ChatGPT review.
   Original plan: (synthetic, no lattice). Write `src/calib_mg3b.py` exactly to the proposal's Part 3:
   the frozen configurations, candidates, eligibility and selection rules, the full decision tree, and the runtime
   benchmark first (200 replicates per candidate and configuration, then R = 4000 or the registered R = 1000 fallback).
   Output to `results/MG3b/calibration/`, including `decisions.json`. Send for ChatGPT review before going further.
4. **Driver: written 2026-10-29 (MG007c)**, `src/run_mg3b.py`, with tests passing. Run it on the PI's machine: `$env:MODELG_WORKERS=9; python src/run_mg3b.py`. It is resumable, and exits with code 2 on a budget stop, 4 if a clip pilot fails, 5 if `decisions.json` is missing. Original plan: the MG3b driver, reusing `src/mg3.py`, with the ledger, the per-batch budget check, full confirmation-field
   storage and the `decisions.json` gate. Tests first. Then the τ/4 timing runs, Part 1, Part 1b if triggered, and Part 2.

Workflow note: ChatGPT reads the public repo's `main`. Claude works on a branch and opens a PR, and the PI merges it, so
review material must be merged before it is sent for review.

## Also pending (not started; PI decisions)
- **Exploratory E1** (registered in protocol v1.3.3): structure factor S(k) of φ_X at Ω ≤ 100, compared with the LNA
  spectrum. A candidate stochastic-Turing-pattern interpretation of the crowded regime.
- **Deferred by the PI:** D1 (noisy existence-window edge, δ ∝ Ω^(−2/3) test), D2 (h = 0.1 repeats; b relative to each
  lattice's window), D3 (MG2 cells at Ω = 10, 30). See `STAGE_LOG.md`, "Deferred work". The local machine makes these
  affordable, but the PI decides when.
- Phases 3–5 (MG4–MG6). The definition of "mass" is a PI decision at Phase 5.

## Things the previous session learned the hard way
- Don't match process names with `pkill -f`/`pgrep -f` patterns that also match your own shell.
- A `str.replace` that silently finds nothing can drop text. Assert that the anchor exists before editing records.
- Stopping early after a first event is valid for first-passage rates, but it changes what a "fraction of detections"
  means. Define the population before running.
- Record every PI decision (budgets, deferrals, approvals) in `GATES.md`/`STAGE_LOG.md` the moment it is made, before
  the data it affects.
