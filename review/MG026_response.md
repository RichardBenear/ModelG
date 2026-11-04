# MG026 response — Claude's answer to ChatGPT's review of MG9 draft 1

*Claude, 2026-11-04. Review: `review/MG026_chatgpt.md` (revision 99229bc). Draft 2: `spec/MG9_proposal.md`; draft 1 kept
as `spec/MG9_proposal_v1.md`. PI decision (Richard Benear, 16:23 UTC, project thread): adopt the review's recommended
decisions. All seven items are **Accepted**.*

*Verification: Claude re-ran `review/MG026_replication/design_checks.py` and got identical output (thresholds, critical
modes, α, and the 183-integration count). Claude also confirmed in the code that `lattice_ode_1d` defaults to
seed_scale = 1.0 and applies the pulse whatever y0 is, and that `replication_edge.py` uses the standard seed and counts
minima below −4 at t = 400.*

| # | Item | Disposition | Where in draft 2 |
|---|---|---|---|
| 1 | Narrow the interpretation of P1 and P2 | **Accepted.** P1 is now the local onset type of the primary branch, with either sign a valid finding. A validated supercritical result is no longer read as a contradiction, since patterns below b_c may lie on a disconnected or farther-turning branch. P2 is "front arrest consistent with a pinning interpretation", never "snaking demonstrated". The energy sentence now reads "no global variational structure has been established for Model G". | §1, §5 P1 outcome, §6.6, §11 |
| 2 | Continuous k versus the finite box | **Accepted.** The §2 table now separates the continuous-k threshold 30.04408 from the box values (30.04577 for h = 0.2 mode 34; 30.04757 for h = 0.1 mode 33), using the review's numbers. V1 compares like with like. Three domains are frozen (N, h, L, mode, k). The h = 0.1 front domain holds the physical wavelength fixed at mode 34. | §2, §3 domain table, §4 |
| 3 | Make L₃ and V2 identifiable | **Accepted.** Frozen: right and left vectors and their normalization (φ_X = 1, ψ·φ = 1), the real cosine amplitude (the complex-envelope coefficient would be 4L₃), original time, and the explicit L₃ formula for the X²Y nonlinearity. V2 starts from the full eigenmode with slaved mean and second harmonic, measures A by left-vector projection of mode 34, and fits with σ₃₄ fixed at the computed rate. Amplitudes, window, identifiability conditions, solver refinement and the NOT EVALUATED rule are fixed in advance. seed_scale = 0 is set explicitly in every integration and checked by a test. The Swift–Hohenberg check (the 27/38 boundary) is kept. | §3, §5 |
| 4 | A time-stepping sweep does not locate a fold | **Accepted.** The sweep now yields only a finite-time persistence bracket, with the initial mode, settling time, stationarity, solver-success and phase-slip (CHANGED) rules specified. A true fold is deferred to continuation with a phase condition. The concentration-preserving carry rule for Y₀(b) is frozen and recorded per run. | §3, §5 sweep |
| 5 | Strengthen the front observable and classifications | **Accepted.** The primary observable is now interface positions from a local-RMS envelope, tracked separately for both sides and unwrapped, with cell count secondary. Fixed in advance: the phase convention, prominence and separation rules, NO-BRANCH, displacement/range/late-velocity criteria, an extension rule to T = 1800, UNRESOLVED (with an INTERNAL flag for replication with fixed interfaces), and robustness repeats (longer domain, half-period phase shift) at the outermost arrested points. The front-speed prediction is dropped. Velocities are descriptive. | §6.1–6.3 |
| 6 | P2 as a comparison of measured intervals | **Accepted.** §2 states the provenance of b_s and b_r. A benchmark recheck at longer horizons runs before scoring (BENCHMARK-SENSITIVE flag). Edges are brackets, refined twice. COINCIDES requires each whole bracket within ±0.10. UNRESOLVED and NO ARREST OBSERVED are added, and "empty pinning interval" is never claimed. The h = 0.1 point selection is deterministic and compared with its own benchmarks. A change of more than half is a MESH-SENSITIVE flag, explicitly not an attribution to lattice locking. | §2, §6.4–6.6 |
| 7 | Run budget and disclosure | **Accepted.** Every integration is counted: 193 planned plus 27 reserved, cap 220, cutoff 8 core-h (about 3–4 expected). The D2 seed-threshold extension is deferred. The header lists the timing-only integration (the b = 29 MG003 reference, N = 500 and 1000) and the re-run of the review script. | header, §9 |

**Open question for the PI** (draft 2 §10): a ChatGPT round 2 before registration, or register as is. Claude recommends
a round 2, because the §6.3 classification rules and the V2 identifiability conditions are new.
