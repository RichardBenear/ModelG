# Stage log summary

Short version of `STAGE_LOG.md` (main as of 2026-11-13, including the PI closures of eight studies recorded that day). Each study is listed once, with its outcome. The log records that belong to it sit under it. For the full detail, see `STAGE_LOG.md`, `GATES.md` and each study's results file.

**Two numbering systems.** *Studies* have short names (MG1, MG3b, MG10-E1, MG11-G4). *Log records* are numbered in order (MG000 to MG162). One study usually takes several records: design, registration, run, review and closure. PR numbers are given up to MG8. From MG039 onward, work was committed straight to `main`, so those studies show commit hashes where the log gives them.

---

## Phase 0–1: published model and continuum

- **MG1: reproduce published Model G.** PASS by PI exception (2026-10-27)
  - MG000: repo scaffold, charter, rules v1.0 (`1bc0015`)
  - MG001: model pinned (Pulver & LaViolette 2013). b = 29 is 3.5% below the Turing threshold. Rules v1.1 add G⇌X (`110bf8b`)
  - MG001a–b: ChatGPT reviews; integral identity ∫φ = 0; MG1 pre-registered (`5d74cfa`, `d56026a`)
  - MG002: MG1 run. 14/15 criteria pass. 1D not yet stationary, so FAIL as registered (`7c4b194`)
  - MG002a: ChatGPT replicated it; PI ruled PASS with documented exception (`ba7176e`)

## Phase 2: 1D stochastic lattice

- **Phase 2 protocol.** v1.0 → v1.3.4 (2026-10-27/28)
  - MG003: protocol v1.0. Narrow single-soliton window b ∈ [28.929, 29.397] (`dd9e9be`)
  - MG003a–c: three review rounds. New Euler-exact kernel; calibrated inference (12 runs/cell); v1.1–v1.3 (`6ecd0bd`, `d5b381c`, `232ee36`)
  - MG003d: v1.3 frozen; MG2/MG3 pre-registered (`0f5615e`)
  - MG004–004a: validations; Amendment 1 (detector threshold) → v1.3.1 (`e5bed86`, `e8f6877`)
  - MG004b–d: V2 FAIL / V4 INCONCLUSIVE; Amendment 2 (decision Ω = 1e8) → v1.3.2; V1–V4, V6 PASS (`e4cb48e`, `a34f15e`, `6011830`)
  - MG004e–g: V7 noise regimes; Amendment 3 → v1.3.3; MG004f glossary added (`b5c944a`; MG004f `8b2f742`)
- **MG2: lattice statistics.** FAIL (2026-10-28)
  - MG005: C3 and C4 FAIL (noise-induced replication); 11.6 core-h
  - MG005a: ChatGPT confirmed FAIL; tracker fix; Amendment 4 → v1.3.4
- **MG3: nucleation rates.** FAIL on C2(a) (2026-10-28)
  - MG006: rates at Ω = 150/120/100; S = 0.073; 19.0 core-h
  - MG006a: ChatGPT review; FAIL kept
- **MG3b: rate estimation follow-up.** CLOSED 2026-10-30. S = 0.064 (0.047–0.082). PRs #1–#11
  - MG006b: proposal drafts 2–3
  - MG007: registered
  - MG007a–b: synthetic calibration and addendum 1
  - MG007c: driver
  - MG007d: production on PC, 26.75 core-h; ChatGPT MG008 reproduced it; PI sign-off
- **Deferred (PI):** D1 noisy window-edge study, D2 h = 0.1 repeats, D3 low-Ω ladder cells

## Phase 3–4: 2D graphs, frozen and evolving

- **MG4: frozen irregular 2D graphs.** FAIL. CLOSED 2026-10-31. PRs #12–#18
  - MG009–009c: design drafts 1–4; registered (`773255d`)
  - MG010: Part 0; nested-dissection LU speed-up
  - MG010 result: every Delaunay and random-geometric run multiplied to 33–37 particles; lattice stayed single; 91.6 core-h
- **MG4c: Delaunay window in b.** W1 + O1. CLOSED 2026-10-31. PRs #18–#23
  - MG011–011c: drafts 1–3; registered (`a946550`)
  - MG012: implementation
  - MG013: single-particle window b ≈ 27.27–28.73; cotangent operator keeps b = 29 single; 11.3 core-h
  - MG013a–b: ChatGPT MG012 review; tolerance deviation recorded; PI sign-off
- **MG5: birth–death evolving Delaunay graph.** PASS + Part B A1. CLOSED 2026-11-02. PRs #24–#32
  - MG014–014e: drafts 1–3; storage rule; registered (`e2a25e3`)
  - MG015–015a: implementation; Windows file-lock fix
  - MG016: Part 0 PASS
  - MG017: first launch interrupted (~65 core-h lost); amendment: 5 workers, cutoff 80
  - MG018–018a: all rates single 4/4; rewriting rescues b = 29; 40 core-h; ChatGPT reproduced it; closed
- **MG6: particle–graph feedback.** DEMONSTRATED, SURVIVES, falloff UNRESOLVED. CLOSED 2026-11-03. PRs #33–#41
  - MG019a–c: drafts 2–3; registered (`106b3e2`)
  - MG020: implementation
  - MG020b/MG021: Part 0 STOP (tick check); amendment 1 (Δ = 0.1)
  - MG021c/d: Part 0 complete; cutoff raised to 110 core-h
  - MG022: 20 runs, 83.1 core-h; coupled arm removes ~67 nodes (t ≈ −35); ChatGPT Windows replay; closed
- **MG7: mass scaling across b.** NOT SCORED (feedback); b = 27.5 CHANGED, 28.5 SURVIVES. CLOSED 2026-11-04. PRs #42–#48
  - MG023–023b: drafts 1–2; registered (`5b44cde`); multiway gate renumbered MG8
  - MG024–024a: implementation; Part 0 PROCEED
  - MG025–025b: 16 runs, 73.0 core-h; ChatGPT MG025 corrections; closed

## b as a control parameter (1D deterministic and stochastic)

- **MG9: subcritical Landau instability.** P1 SUBCRITICAL, P2 UNRESOLVED. CLOSED 2026-11-04. PRs #49–#55
  - MG026–026c: drafts 1–3; registered (`7810b8b`)
  - MG027: implementation; amendment 1 (demodulation envelope)
  - MG028: L₃ = +0.0137; fronts creep at 29.25–29.35 and invade from 29.40; 0.5 core-h; closed
- **MG9c: fate of the slow fronts.** Q1 UNRESOLVED, Q2 NOT ELIGIBLE. CLOSED 2026-11-04. PRs #56–#60
  - MG029: drafts 1–2; registered (`372336b`)
  - MG030: residual 1.08e-5 vs 1e-5 limit; slow rigid drift; MG030b centre diagnostic; closed
- **MG9d: cause of the drift.** MESH EFFECT SUPPORTED. CLOSED 2026-11-05. PRs #60–#65
  - MG031: drafts 1–3
  - MG032: implementation; registered (`69ad83e`)
  - MG033: 7 arms, 0.08 core-h; drift depends on sub-grid position; closed
- **MG9b: stochastic barrier S(b).** NOT RESOLVED. CLOSED 2026-11-05. PRs #67–#75
  - MG034: drafts 1–3; registered (`185aceb`); amendment 1 (budget truncation)
  - MG035: 280 runs, 26.3 core-h; dS/db = −0.023 (95% CI −0.094 to +0.048); post-hoc event-shape note archived as a hypothesis; closed

## Phase 6: update order

- **MG8: multiway / rewrite order (Central Question 5).** SENSITIVE at τ₀; refinement NOT RESOLVED BY THIS LADDER. CLOSED 2026-11-06. PRs #76–#83
  - MG036: drafts 1–4; registered (`3ab4ec8`)
  - MG037–037a: implementation; Part 0 attempt 2 accepted; review fixes
  - MG038–038a: 55 runs, 3.1 core-h; sweep order dissolves the particle at t ≈ 8; ChatGPT review; closed

---

## Role change: MG039 (2026-11-06)

ChatGPT/Codex became research and implementation lead. Claude moved to optional review. `RESEARCH_PROTOCOL.md` was added. From here on there are no PRs, and records are shorter and denser.

## MG10: explicit finite reservoir (2026-11-06 to 11-08)

- **MG10: first deterministic batch.** UNRESOLVED. CLOSED 2026-11-13
  - MG039: Draft 1 (explicit A/B stocks)
  - MG039b: PI decisions; at most 8 local workers
  - MG040: numerical plan; capacity grid to 1e6
  - MG041: qualification stopped by a spec error
  - MG042: corrected and registered
  - MG043: 64 runs, 0.21 core-h; high capacities recover; C = 1 numerical disagreement forces UNRESOLVED
- **MG10-C1: homogeneous reference and amplification.** Supported. CLOSED (MG046)
  - MG044: registered · MG045: results · MG046: PI closure
- **MG10-S1: finite-reactant sampler.** QUALIFIED ON TESTED SMALL SYSTEMS. CLOSED (MG049)
  - MG047: registered · MG048: 96/96 comparisons pass · MG049: PI closure
- **MG10-S2: actual-lattice spatial noise.** QUALIFIED WITHIN S2 SCOPE. CLOSED (MG052)
  - MG050: registered · MG051: results · MG052: PI closure
- **MG10-D1: long-horizon finite-supply pilot.** Complete, with a ~0.3% source-sampler deficit found after the run. CLOSED 2026-11-13
  - MG052: preparation · MG053: registered and launched · MG054: 80 runs, 33.3 core-h
- **MG10-D2: corrected source comparison.** COMPLETE DESCRIPTIVE PILOT. CLOSED (MG057, 2026-11-07)
  - MG055: source fix; registered (`ff12250`) · MG056: 32/32 complete · MG057: PI closure
- **MG10-P1: single-motif measurement.** QUALIFIED ON SYNTHETIC FIXTURES. CLOSED (MG060)
  - MG058: registered · MG059: 90/90 assertions pass · MG060: PI closure
- **MG10-P2: model-noise observability.** Unseeded match caps supported (0/128 each). CLOSED (MG064)
  - MG061: preparation · MG062: registered (`eaef80d`) · MG063: results · MG064: PI closure
- **MG10-P3: temporal persistence.** COMPLETE TEMPORAL ESTIMATION. CLOSED (MG068, 2026-11-07)
  - MG065: preparation · MG066: registered (`d0e9fcc`) · MG067: results · MG068: closure; one-paragraph summary rule added
- **MG10-E1: resource-location intervention.** COMPLETE; interaction supported. CLOSED (MG072, 2026-11-08)
  - MG069: preparation · MG070: registered (`1249d0a`) · MG071: 19.8 core-h; public push blocked until PI approved publishing · MG072: PI closure
- **MG10-E2: reversible supply.** INCOMPLETE (stopped by the forecast rule), preserved unchanged
  - MG073: development · MG074: registered (`e2fc7f4`) · MG075: stop diagnosed
- **MG10-E2-C1: E2 continuation.** COMPLETE RESTORATION ESTIMATION (coarse +0.103, fine +0.141). CLOSED (MG078)
  - MG076: amendment; registered (`eacd1ac`) · MG077: 28.4 core-h in total · MG078: PI closure

## MG11: environmental networks (2026-11-08/09)

- **MG11-N1: fixed environmental network.** QUALIFIED WITHIN N1 SCOPE. CLOSED (MG082)
  - MG079: build · MG080: registered · MG081: 6/6 + 28/28 bounds pass · MG082: PI closure
- **MG11-G1: connectivity vs supply.** COMPLETE; positive B contrast (+1.52 pp). CLOSED (MG086)
  - MG083: build · MG084: registered · MG085: results · MG086: PI closure
- **MG11-G2 Step 1: nonlinear coupling build.** Qualified, but the original κ = 1000 mixing check FAILED and is kept on record. CLOSED 2026-11-13
  - MG087: build; mixing check fails (0.00521 > 0.005)
  - MG087a: separate κ = 10000 point passes
- **MG11-G2 Pattern-Response.** RESOLVED CONNECTIVITY CONTRAST (+0.0545). CLOSED (MG091)
  - MG088: preparation failed a tolerance check · MG088b: corrected · MG089: registered · MG090: results · MG091: PI closure
- **MG11-G3: B-transport switch.** RESOLVED ONGOING B-TRANSPORT CONTRAST (+0.0333). CLOSED (MG095)
  - MG092: preparation · MG093: registered · MG094: results · MG095: PI closure
- **MG11-G4 Step 1: stochastic coupling.** QUALIFIED WITHIN SMALL-SYSTEM SCOPE. CLOSED (MG099)
  - MG096: build · MG097: registered · MG098: 196/196 bounds pass, 0.19 core-h · MG099: PI closure
- **MG11-G5: stochastic full-lattice pilot.** UNRESOLVED (response and step agreement). CLOSED (MG103, 2026-11-09)
  - MG100: preparation · MG100b: ready · MG101: registered · MG102: 80 jobs, 11.3 core-h · MG103: PI closure

## MG12: explicit hypergraph engine (2026-11-09 to 11-12)

> **Switch to the tick-based approach: MG104–MG105 (2026-11-09).**
>
> **Before (MG2–MG11):** the reactions ran as a continuous-time stochastic process (the Model G rate law). Each node held token counts, and the code advanced time in small numerical steps: the tau-leap τ, or the graph-rewrite interval Δ in MG5–MG7. Those steps were only approximations to continuous time. In any one step a token could, in effect, take part in several reactions, and the step had no physical meaning of its own.
>
> **After (MG12, from MG105):** a "strict tick" executor works on explicit, individually tracked tokens. Every possible rule match is read from the same start-of-tick state, and each token takes part in at most one event per tick. That includes catalysts, by the PI's choice at MG105. Products and moved tokens wait until the next tick. The old continuous-time law is kept as a reference. How closely ticks match it is measured, not assumed (MG124–MG158; ticks 1e-4 and 5e-5 equivalent to 2e-5 within scope at MG158). No physical duration is assigned to a tick.

- **MG12 roadmap and H0 contract.** Static checks passed
  - MG104: six-stage roadmap · MG105: H0 rewrite semantics; 101 exact checks pass
- **MG12-H1: explicit engine.** QUALIFIED (16/16 bounds). CLOSED (MG109)
  - MG106: development · MG107: registered · MG108: passed · MG109: PI closure
- **MG12-H2: higher-arity transport.** QUALIFIED (24/24 bounds). CLOSED (MG113)
  - MG110: development · MG111: registered · MG112: passed · MG113: PI closure
- **MG12-H3: reaction + transport, from count level to continuum.** Several sub-studies:
  - **H3a: bridge.** Passed (MG114)
  - **H3b: reference bridge.** First run FAILED; corrected completion PASSED. CLOSED (MG123)
    - MG115: preparation · MG116: registered · MG117: FAILED (error 1.04e-4 > 1e-5) · MG118–119: tolerance check passed · MG120–121: completion plan; registered · MG122: passed · MG123: PI closure
  - **Mixed finite-tick qualification.** PASSED (max 0.084 < 0.10). CLOSED (MG127, 2026-11-10)
    - MG124: preparation · MG125: registered · MG126: passed · MG127: PI closure
  - **Exchangeable-token efficiency.** PASSED (max 0.080 < 0.10). CLOSED (MG133)
    - MG128–130: development · MG131: registered · MG132: passed · MG133: PI closure
  - **Sparse local-hazard engine.** PASSED (max 0.082 < 0.10; 8 million stress updates). CLOSED 2026-11-13
    - MG134: ring inventory · MG135: engine · MG136: combined development batch · MG137/137a: preparation (an oracle bug was fixed) · MG138: registered · MG139: passed
  - **Multi-tick pilot.** Integrity PASSED (60.8 million events). CLOSED (MG143, 2026-11-11)
    - MG140: preparation · MG141: registered · MG142: passed
  - **Tick-resolution design**
    - MG143/144: deterministic operator screen · MG145/145v2: design; 5% contrast margin chosen · MG146–147: estimator and compact writer
  - **Finer-tick precision pilot.** Complete; late-Y variance too high. CLOSED (MG150)
    - MG148: registered · MG149: results
  - **Fixed Y-blend estimator.** Variance benefit SUPPORTED (ratio 0.335). CLOSED (MG155)
    - MG151–152: preparation · MG153: registered · MG154: results
  - **Paired tick coupling.** Development passed; pilot benefit SUPPORTED. CLOSED 2026-11-13
    - MG156: development passed · MG157: longer pilot, 4.6 core-h
  - **Tick accuracy.** EQUIVALENT (ticks 1e-4 and 5e-5 vs 2e-5; 24/24 intervals), 48.7 core-h. CLOSED (MG159, 2026-11-12)
    - MG158: readiness; approved; results
  - **Pattern persistence P1.** Retention NOT observed (0/4 prepared pairs to T = 0.1). CLOSED 2026-11-13
    - MG160: preparation · MG160-P1: registered; results; amendment to keep large data local
  - **Noise-aware shape instrument.** Synthetic and replay benchmarks supported. MG161 CLOSED 2026-11-13; MG162A/B CLOSED at MG163
    - MG161: synthetic benchmark · MG162A: empirical replay challenge · MG162B: fresh pilot (in git; not yet in `STAGE_LOG.md`)
