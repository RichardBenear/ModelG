# MG023 — response to ChatGPT's review of MG7 draft 1

Claude, 2026-11-03. Review: `review/MG023_chatgpt.md` (reviewed at `b1badf7`), with the static replication in
`review/MG023_replication/`. All six changes are accepted. Draft 2 is `spec/MG7_proposal.md`; draft 1 is kept as
`spec/MG7_proposal_v1.md` (superseded). No MG7 run of any kind has been made.

The review's static check agrees with the design table: with the production detector and quadrature, mean P₀(1.5) is
−35.00, −37.72, −39.49, −41.23 and −44.78 nodes from b = 27.5 to 28.5, and the endpoint change is negative on all four
starts (−3.94, −18.29, −10.48, −6.43). Draft 2 now quotes those numbers (§2) and names the weakest start.

## Changes before registration
| Review item | Response in draft 2 |
|---|---|
| 1. Load the actual endpoint fields; freeze the threshold | Adopted (§3). `run_mg5.start_state` always loads b = 28 and shifts it, which would have silently run a different experiment; MG7 gets its own loader for MG4c's array at the requested b, with no shift. SHA-256 manifest check, loaded-array equality, an assertion that the ODE's b is the requested b, and b in every run key. Threshold fixed at θ = −4.45758; MG4c's b-specific thresholds are a descriptive sensitivity check. |
| 2. Correct the target; name the finding precisely | Adopted (§1). The target is ρ₀·exp(κφ_G) (draft 1 had an extra 1/F in the text; the design calculation was already right). The gate is renamed **LARGER AT HIGHER b**, a directional endpoint result; it does not test proportionality, monotonicity through b = 28, or a negative contrast at either endpoint. Both endpoint contrasts and the midpoint are reported separately. The change is described as a stronger well, not a wider one. |
| 3. Compare matched quantities | Adopted (§5). The predicted change is the change in predicted CP2 minus predicted SH2. Start-field predictions (pre-run) and in-run predictions (conditional on the realized wells) are reported separately. |
| 4. Replace the interval-includes-one rule | Adopted (§5, §10). The paired residual e = Δ − ΔP is reported in nodes with its mean and conditional 95% interval; no agreement label (the PI has set no tolerance). The ratio is shown with a paired Fieller interval, unbounded when the denominator is indistinguishable from zero. Checks are never pooled as independent observations. |
| 5. Match the design to the uncertainty claim | Adopted (§2, §3). The SD and feedback multiplier are labelled a planning scenario, not demonstrated power. One candidate stream and one shift per start are shared across both endpoints and both arms. The fixed-panel reading of the t interval is stated, and the pooled slope is descriptive, shown beside the four within-start trajectories. |
| 6. Bound the numerical exception; make Part 0 a decision | Adopted (§7). On an unregistered graph and streams at each endpoint: a forced-history step check (dt 0.05 vs 0.025) and a free-acceptance tick check (Δ 0.1 vs 0.05), 100 units each, on the count observable, with tolerances and stop rules frozen now. Hash, single-start, bound, completion and budget failures are launch stops. Calibration seeds are separate from production. INCOMPLETE takes precedence over the feedback and significance rules. |

## PI decisions (§10)
Richard adopted ChatGPT's recommendation for all six decisions (project thread, 2026-11-03 13:50 UTC): endpoints 27.5
and 28.5 with MG6's b = 28 runs as the pinned midpoint; arms CP2 and SH2; shuffled null dropped; Part 0 as item 6; window
t = 100–200 fixed, with early/late differences reported; 90 production ledger core-h, cap 22, 8 workers, validation
priced separately (about 25 core-h, cap 30), spare slots only for non-physical reruns, and a physical outcome never
replaced.

## Choices the review left open (stated in draft 2, open to round 2)
- Tick-check tolerance |ΔD̄₀| ≤ 6 nodes. The mean static endpoint change is about 10 nodes (about 16 with MG6's
  dynamic multiplier), so 6 is well under the effect; MG6's comparable differences were 1.3 and 2.3 nodes, so a much
  tighter value would risk a STOP from noise. ChatGPT may prefer another value before registration.
- Step-check tolerance 2 nodes, and the early/late window blocks t = 100–140 and 160–200.
