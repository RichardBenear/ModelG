# MG12-H3 initial bridge and pattern-reference readiness — MG114/MG115

H3 has begun by checking that the previously qualified reaction and transport layers still work together, rather than assuming that their separate successes automatically establish the combined model. That bridge passed across all reachable states of four tiny systems, and the known ring equations were recovered with the expected finite-count correction. We also built an incidence-derived population-limit solver and verified its short trajectories against the earlier closed pattern calculations. The next deterministic pattern-reference run is ready. This establishes implementation readiness, not recovery of a particle in the strict finite-tick engine, a new geometry effect or a quantum substrate. H2 is signed off and closed; no full H3 late-window trajectory has yet been run.

**Status: mixed-law development PASS; full deterministic pattern-reference experiment READY FOR CONCRETE PI APPROVAL.** H3 is not complete.

## Mixed reaction-and-transport bridge

The frozen H2 identity-rich engine was used unchanged, now with all six P17 reaction channels enabled alongside transport cells. A reaction can occur only when its reactants are present; a globally enabled channel need not be accessible in every fixture. Finite resources and inactive retained Z/W remain unchanged. Cases were A+X on three sites, B+X on three sites, X+X+Y on three sites with overlapping cells, and A+B+X+Y on two sites. Active and resource transport share the same source snapshot and token reservations as reactions, including catalyst participation.

| Fixture | Reachable count states |
|---|---:|
| A+X, three sites | 72 |
| B+X, three sites | 36 |
| X+X+Y with overlapping cells | 300 |
| A+B+X+Y, two sites | 488 |

All **896 states** matched the independently written count-generator reaction/transport rates within the frozen 1e-12 relative/absolute bounds. Reference probability checks passed. The largest generator row-sum residual was about 2.45e-15.

Exact and strict maxima 0.005/0.001 used 128 independent paths per fixture/method: **1,536 development paths**. All four exact-executor endpoint CDF smoke comparisons passed the registered bound **0.161591351**. Strict-tick statistics are descriptive; no mixed-law statistical accuracy or convergence claim is inferred from that small sample. No setting or replication count was chosen from outcomes.

All paths passed literal chemical/cell-role, source-snapshot, catalyst-reservation and token-total checks. The saved-data reread reconstructed **919 accepted events** and all observation frames, and separately recomputed endpoint CDF distances. This shares the role-audit routine with sampling and the same author/model constants; the count generator is separate from primary incidence matching/execution. It is not independent external replication.

Development took **13.505212 seconds**, below its 300-second single-worker cap. Source, fixture and seed plans were published at ce25977 before sampling/reference integration. No failure, retry, fixture change or accuracy-threshold amendment occurred.

## Ring kinetics: what the count/continuum bridge established

An eight-site ring probe with positive G/X/Y/A/B concentrations was evaluated at token scales 1, 10, 100 and 1,000, without expanding individual token identities. The count-generator conditional drift differs from the existing continuum RHS only by its known finite-count feedback correction: subtract x*y/Omega_count from X drift and add it to Y drift.

| Token scale | Largest difference from continuum | Residual after the specified correction |
|---|---:|---:|
| 1 | 2.000000 | 7.11e-15 |
| 10 | 0.200000 | 9.94e-15 |
| 100 | 0.020000 | 1.03e-14 |
| 1,000 | 0.002000 | 9.55e-15 |

Global seven-species drift conserved the closed total to floating-point precision. This is conditional generator algebra at specified states, not convergence of full stochastic trajectories or of a mean-field expectation closure. It also does not establish that finite-tick catalyst reservations preserve the continuum dynamics at every population. Those are separate limits/tests.

## Deterministic population-limit solver preparation

The new solver builds conservative transport operators directly from incidence cells and independently implements the population-limit P17 reaction increments. It does not call the old reaction/solver RHS. It evolves seven concentration fields plus A/B external-import accounts. Closed cases have zero imports; maintained controls are explicitly open and record the replenishment required to hold A/B fixed. Z/W remain retained outputs, with no invented regeneration law.

All 24 candidate jobs were run only to **T=1**, before any new late-window trajectories. They comprise eight finite prepared contexts at two numerical settings, and maintained/isolated H/S controls at two settings. The 500-site ring, donor states, original references and old fine archives are fixed. Base BDF is rtol1e-9/atol1e-11/maxstep0.1; fine is 1e-11/1e-13/0.05. Source/reference hashes were published at 056cd61 before these integrations.

Preparation results:

- Largest scaled field difference from the closed fine reference: **7.37234316261e-07**, below **1e-5**.
- Largest base/fine scaled field difference: **7.38810960854e-07**, below **1e-5**.
- Largest relative population-minus-import balance residual: **1.70828175704e-15**, below **1e-10**.
- Largest initial RHS difference from the separately evaluated old equations: **4.42810232926e-11**, within the frozen 1e-9 absolute/1e-12 relative check.
- Deliberately corrupt balance and negative-field records were rejected; incomplete/missing outputs did not produce a complete/qualified result.
- Total preparation time: **5.523560 seconds**, below its separately frozen 300-second cap. No numerical failure, retry or criterion change. A syntax typo was corrected during uncommitted preflight, before registration/integration; no scientific outcome was affected.

The complete 24-job analysis path also passed using separately labeled **held archived reference fields**, with no new ODE integration or random sampling. It exercised full late-window projection, independent dense projection, paired-setting and continuity reporting. The held interaction was 0.033337305924; new and archived inputs here are intentionally the same fields, so this validates the reporting pipeline rather than a new physical outcome. This added 4.491327 seconds, for 10.014887 aggregate preparation seconds, below the 300-second preparation allowance. Held artifacts are labeled in population/preparation/held_pipeline_results.json.

## Concrete proposed reference run

| Item | Fixed scope |
|---|---|
| Finite contexts | Two preparation origins x H/S x LOCAL/CROSS B destination |
| Controls | Maintained and isolated H/S; imports metered for maintained cases |
| Numerical settings | Base and fine, giving **24 deterministic jobs** |
| Horizon | Finite branches post15, absolute10–25; controls T0–25 |
| Output | Quarter-unit seven-species fields and A/B import accounts |
| Reference criterion | Every new frame scaled field error <=1e-5 vs its closed fine archive |
| Refinement criterion | Every base/fine paired field difference <=1e-5 |
| Pattern-response criterion | New fine late X interaction agrees with recomputed closed G3 interaction within1e-5 projection units |
| Account criteria | Finite fields, min>=-1e-9; relative global population/import residual<=1e-10 |
| Resource location | One local worker, one numerical-library thread |
| Measured forecast | **293.65 seconds (~4.89 minutes)** including conservative compute/IO scaling |
| Requested hard cap | **3,600 elapsed seconds (one single-worker hour)**; no automatic retry |
| Retention | Full raw/reference provenance, metadata, prefixes, analysis and report in public repository |

The forecast separately scales numerical work and repeated-prefix IO, including the growing checkpoint volume; it remains an extrapolation. Existing prepared parents are not reformed or re-rounded. Exact archived fine reference hashes and all source/configuration choices are frozen. There is no RNG and no new stochastic seed selection in this deterministic batch. No external review is requested by the lead; Richard retains that decision.

The execution and evaluation code is built: src/run_mg12h3_population.py registers only with approval and passing unchanged preparation, refuses an existing attempt, saves prefixes during integration, checks each case and compares paired settings/independent projections. Concrete production registration has **not** been created. The start of H3 authorized development; this packet makes the new mathematical-limit pattern experiment, criteria and budget reviewable before launch under RESEARCH_PROTOCOL.md workflow3.

## What success would and would not establish

A pass would recover the earlier ring pattern trajectories and G3 geometry interaction in the incidence-derived deterministic limit. It would establish that the translation retains known mathematical-limit behavior on those prescribed substrates. Sampled reference continuity would be reported descriptively, including the maintained seeded control; larger projection after loss of matching would not be called a longer-lived particle.

It would not establish a particle in the strict finite-tick hypergraph process, stochastic robustness of G3/G5, a pattern on a new higher-arity/3D substrate, emergent geometry or physical vacuum supply. Formal mixed-law accuracy, efficient exchangeable-token sampling with catalyst reservations, population/tick limits and a sufficiently precise stochastic pattern study remain further H3 work. G5 remains closed and unresolved.

Read spec/MG12H3_bridge_development.md, spec/MG12H3_pattern_bridge.md, bridge/results.json, saved_evidence_check.json and population/preparation/results.json for reproducible inputs. Full preparation/raw-reference retention verification is recorded in retention_verified.json. No full reference run or new scientific closure is implied by this readiness report.
