# Current handoff — 2026-11-08 (MG076)

- **Roles:** Richard Benear is PI; ChatGPT / Codex leads planning, implementation, validation, authorized runs, analysis and documentation. Claude/other external review is optional. [Active protocol](RESEARCH_PROTOCOL.md).
- **Last closed gate:** MG10-E1, PI closure MG072 (2026-11-08): Richard: “I sign off”. CLOSED as COMPLETE INTERVENTION ESTIMATION;32 parents/128 branches, separate160-file audit and full published archive verified. Direction-consistent seed-by-location field interaction supported within registered pilot scope; no physical mortality or reservoir-origin claim.
- **Standing result convention:** Every new stage result report includes a one-paragraph plain-language project-level summary near the beginning, explaining its contribution and limits; use that style in results presented to the PI. See RESEARCH_PROTOCOL.md / AGENTS.md (MG068). No automatic backfill/rescore of earlier frozen reports.
- **Current result:** MG10-E1 CLOSED as COMPLETE INTERVENTION ESTIMATION (result MG071, closure MG072).32/32 parents/128/128 branches, no failure/timeout/stop;160 raw audits PASS. Primary seed-by-local-vs-remote fixed-origin X interaction coarse-0.833949 CI[-1.102872,-0.565027], fine-0.856842 CI[-1.001798,-0.711885], n8/step; registered direction-consistent pilot evidence. Local seed-excess nearly0 vs remote~.86-.87, no physical mortality claim. All intervention S arms fail endpoint reference classifier including remote/uniform, so matching loss cannot define death.19.808474core-hours under36, elapsed~153.19minutes including audit. Completed Nov7~22:37 local. See results/MG10E1/MG10E1_results.md. No new production authorized; do not restart E1. Richard explicitly approved public publication of E1 results/full raw archive on2026-11-08, resolving the earlier approval block. Full results and160 raw archives (290504596bytes) published at18adea684d2025e2f7433044ba9ec67d165cc5d1; all published blobs match audited SHA256 and local bytes. See results/MG10E1/retention_verified.json. PI scientific closure recorded2026-11-08 (MG072). E2 is separately approved/registered below.
- **Active stage:** MG10-E2 APPROVED/REGISTERED MG074. Richard approved design,40-core-hour cap and public results/full-raw publication.32 fixed-T5 parents then128 postT25 forks, local continuous/local restore/remote restore/sham, post5 switch, H/S/two steps/8blocks. Primary late20..25 restoration seed interaction; movement/shape diagnostics secondary, no physical revival/quantum memory claim.20 actual-writer held replays and complete/missing/prefix/empty raw-file/dense statistical audits PASS; total preparation allowance336.239390seconds<900. Exact matrix/source/reference/seeds frozen; publication before sampling required. Forecast29.558218core-hours/~3h42 on8. Check results/MG10E2/run_state.json and parents/branches ledgers for live status. Automatic analysis/audit after simulation, full raw upload/verification explicitly approved; PI closure separate. No automatic restart/retry.
- **E2 stopped/incomplete:** Automatic analysis/audit finished2026-11-08T15:32:56.627145+00:00 (~08:32 Mountain), no active workers.32/32 parents complete;1/128 branches complete,7 retained prefixes,120 not started. BUDGET-SHORT: FORECAST stop, not40core-hour exhaustion;5.452379core-hours spent. No complete paired blocks; primary n0/step, no recovery inference. Raw audit40files PASS within incomplete scope. Do not restart or alter frozen criteria automatically; prospective budget/controller amendment and approval needed. Existing run_state/evaluation/independent_check and all raw files retained. Prior session97858 is finished.
- **Active continuation MG076:** Richard approved144-core-hour aggregate amendment (includes original5.452379),5400/job,8workers,forecast admission stop drains active jobs. E2-C1 separately registered in results/MG10E2C1; original stopped E2 remains incomplete/immutable.32parents+1complete reused;7original-seed full replays require byte-exact saved-prefix proofs;120unstarted same seeds. Full mixed-provenance complete/missing/prefix audit PASS (77.777412seconds<300), original tamper rejected, budget sum preserved. Exact original job matrix/source/original-file hashes frozen. Forecast~11h17remaining8workers,95.763055aggregate core-hours; not guaranteed. Commit/push registration before launch. Check new run_state/branches ledger/progress; never restart original E2. Public full raw retention already approved; PI closure separate.
- **E2-C1 live execution:** Started 2026-11-08T16:07:54.373455+00:00 (Nov8~09:07 Mountain), registrationeacd1ac pushed first. DriverPID17180, persistent execution session46751. Eight active workers all advance checkpoints;7replays match original prefixes so far,1fine unstarted. No stop. Original5.452379core-hours is prior charge within144aggregate cap. Use results/MG10E2C1/run_state.json, branches/ledger.json and progress; never restart original E2 or current continuation. Forecast~11h17remaining at8, uncertain. Automatic combined analysis/separate audit, then approved public archive verification and PI closure.
- **Archive baseline:** [StochasticVacuum](https://github.com/RichardBenear/StochasticVacuum), private, commit `2161211ccc9d8427af98ebeefcabf7f379b83a37`. Reuse selected machinery with provenance and new validation.
- **Time semantics:** MG10 proposes the existing continuous-time model for its first extension. A fundamental-tick law remains a separate modeling decision.
- **Q5:** MG8's frozen record remains open as scored. A possible reframing of Q5 was discussed with the PI, but no new Q5 study or rescoring is registered here.
- **Authority:** this current handoff and RESEARCH_PROTOCOL.md supersede the operational next-step and role directions in the historical snapshot below. Frozen scientific criteria remain binding within their scopes.

---

## Historical handoff snapshot (preserved; not the current task list)

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
| MG4–MG8 | Irregular graphs, evolving hypergraph, gravity-like deficit (MG6), mass scaling (MG7), multiway (MG8) | see GATES.md; MG7 CLOSED 2026-11-04 (NOT SCORED, feedback) |

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
   **MG4 design: draft 4 written 2026-10-30 (MG009b, after two ChatGPT reviews; noise-free only, PI decision)**, `spec/MG4_proposal.md`; **registered 2026-10-30 (MG009c)**. **MG4 CLOSED 2026-10-31: FAIL.** Every Delaunay and RGG run at b = 29 multiplied (33–37 particles fill the torus); the T control stayed single. Both C3 parts fail (the C3(b) design flaw is recorded). ChatGPT accepted it (MG010), and the PI signed off. Next: PI decision on the follow-up. **MG4c CLOSED 2026-10-31: W1 + O1.** **MG5: REGISTERED 2026-10-31 (MG014e; draft 3 frozen at `e2a25e3`).** Code and pre-run specification done (MG015). **Part 0 passed all checks (MG016): proceed with dt = 0.05, Δ = 0.2; restated budget 33.3 core-h, under the 40 core-h cutoff.** First launch interrupted; amendment MG017 (controls first, 5 workers, cutoff 80 core-h). **MG5 CLOSED 2026-11-02: Gate PASS (survives at ν = 0.002, 0.02, 0.2; PINNED, PINNED, MOBILE), Part B A1**; ChatGPT MG018 accepted, PI signed off. **MG6: REGISTERED 2026-11-02 (MG019c; draft 3 frozen at `106b3e2`).** Code, tests and pre-run specification done (MG020). **MG6 Part 0 STOPPED on the tick check (2026-11-02); amendment MG021 (Δ = 0.1, part (b) governs; PI decision) in `spec/MG6_amendment_1.md`, draft 2 after ChatGPT review MG021 (narrow centre-clause waiver; guard fixes).** Part 0 completed under MG021 (checks 5–8 pass; projected 90.7 core-h); PI raised the cutoff to 110 core-h. **MG6 CLOSED 2026-11-03: gate DEMONSTRATED, feedback SURVIVES, falloff UNRESOLVED (MG022; ChatGPT MG022 accepted; PI signed off).** Next: PI decides the follow-up (deferred: mass scaling b-scan, falloff, curvature, shuffled null, ν = 0.02).

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
