# Handoff — state of the project and next steps (2026-10-29)

Written by the Claude session that ran MG1–MG3 in the cloud, for the Claude Code session that continues on the PI's own
machine (10 cores). Everything referred to here is committed. Read `CLAUDE.md` first for the rules and conventions.

## Where things stand
| Gate | Question | Result |
|---|---|---|
| MG1 | Continuum PDE reproduces the published stationary soliton | **PASS by PI exception** (the registered FAIL on 1D stationarity is preserved) |
| MG2 | Token rules converge to the lattice mean field as Ω → ∞ | **FAIL as registered** (C3 at Ω = 100; C4 noise-driven replication at b_mid). The core convergence tests C1 and C6 PASS; C5 PASS |
| MG3 | Particles nucleate spontaneously from intrinsic noise | **FAIL as registered** (C2(a): 0.71, 0.71, 0.62 confirmed vs ≥ 0.80). C1, C3, C4, V5 PASS. Budget overrun 19.03/18 h recorded as an exception |
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

## Immediate next steps
1. **Check the local environment.**
   - Python 3.13 with numpy, scipy and numba. The cloud used numpy 2.5.3, scipy 1.18.1, numba 0.68.0; ChatGPT's Windows
     replay used numpy 2.2.4, scipy 1.16.1, numba 0.68.0. `requirements.txt` lists the packages.
   - Run the tests: `python tests/test_track.py`, `python tests/test_mg3.py`.
   - Measure the cost per time unit with 9 workers. A short run of one cell into a scratch output folder works:
     `MG3_OUT` can point the MG3 driver elsewhere.
2. **Revise the follow-up proposal** `spec/MG3b_followup_proposal.md` (currently a DRAFT written before the full review)
   to ChatGPT's §8 recommendations in `review/MG006_replication/README.md`:
   - a fixed Ω grid around 120–160 with interior values (e.g. 130, 140); Ω = 100 only as a diagnostic boundary; no 4/4
     adaptive trigger;
   - sample sizes from a precision/budget calculation, well above 5 events per fitted cell, with run caps;
   - **censored first-event likelihood** fits; calibrate the whole procedure, including any event-dependent extension,
     by simulation; treat 1/n as estimated information;
   - burn-in and time step studied together: full startup survival curves; landmarks at 50/100/150 from the same
     histories; τ against τ/2 with adequate samples, adding τ/4 if they disagree;
   - risk-set sizes and late-window exposure reported, with windows chosen so the hazard check can actually diagnose
     time dependence;
   - C2(a)-type detector questions as a **separate** study with a defined population (pre-first-event or fixed window).
     Save **every** confirmation field (including excursions and unresolved), tracker summaries, screen checkpoints and
     a budget ledger;
   - any new confirmation threshold is a new design decision, not a re-score of MG3.

   With 9 local workers the budget is much less of a constraint; propose run caps and a core-hour estimate. Then: send
   it for ChatGPT review (the PI relays), revise, get PI approval, register it in `GATES.md` before running.
3. Write the follow-up driver, reusing `src/mg3.py` (simulation, confirmation, first-passage bookkeeping). Add
   per-run storage of all confirmation fields and a budget ledger. Tests first.

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
