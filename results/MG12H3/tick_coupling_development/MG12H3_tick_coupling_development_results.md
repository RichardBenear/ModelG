# MG156 — paired tick-resolution development results

The new method passed its development checks and gave a promising first indication that we can compare tick resolutions with much less random noise. It runs two copies of the model using shared events wherever their rates agree, while retaining each copy's own rates, token reservations and tick boundaries. Exact small-system checks and comparison with the existing engine found no discrepancy. On a very short 500-site test, paired differences were much less variable than differences between independently simulated runs. This supports developing the method further; four pairs and a horizon 100 times shorter than the intended study do not establish the eventual variance benefit, tick equivalence, particle persistence, or a physical interpretation of the reservoirs.

## Outcome and scope

**Development checks passed.** Both bounded phases completed; no failures or retries. Source/specification/seed plan were published before draws at `4e81d77`. Frozen marginal engines were not edited. PI closure remains pending.

This is classical statistical coupling between two simulations. Each model remains a separate strict-tick system: arrivals and products cannot participate again until its next tick; catalysts are reserved for the rest of its tick. Sharing a mark does not transfer tokens between runs or add an interaction to the model.

## Why the individual laws are preserved

For a given channel with current rates a and b, the common-event rate is min(a,b), the left-only rate is a-min(a,b), and the right-only rate is b-min(a,b). Each side therefore retains its own accepted-event intensity. Equivalently, the joint generator acting on an observable of either side alone equals that side's original generator. Other-side events leave that observable unchanged. Only the side whose tick has ended resets its availability; the calendar uses exact rational versions of the decimal tick values. Exponential memorylessness permits redrawing the joint clock when either boundary interrupts it. These arguments support the algorithm; the following checks test this implementation.

## Phase A: correctness controls

- 896 saved-state checks verified the shared/residual rate identities.
- Six scripted cases covered the R4 catalyst reservation and reset, a one-sided reaction, exact-boundary exclusion, empty populations, deadline and event-cap guards.
- 128 identical-tick paths matched the frozen standalone engine's final counts, event-channel sequence and final random-generator state. Event times agreed within 1e-12 model-time units. Both paired copies matched each other.
- 8,192 pairs across four mixed-reaction/transport fixtures were compared against independently constructed, identity-level, exact one-tick probability distributions. The partner initially had one X token moved to another site, exercising unequal rates. All eight marginal null-consistency screens passed. Observed total-variation distances ranged from 0.00915143 to 0.0324014; their predeclared bounds ranged from 0.0687933 to 0.115599. These conservative screens are **not** a qualification that TV is below 0.10.
- 4,096 pure-transport pairs used nonnested ticks 0.005 and 0.002, observed at 0.01 and 0.02. All four marginal endpoint checks agreed with the analytic law within the predeclared bound 0.0331220. The largest absolute discrepancy was 0.0159928.

Phase A took 16.3826 seconds, below its 60-second ceiling.

## Phase B: exploratory 500-site test

The canonical Omega=10 state was simulated to model time 0.0001. Each comparison used four coupled pairs and four independent reference pairs with separate seeds. The reference tick was 0.00002. Only direct projected G/X/Y contrasts were measured; the compensator estimator and fixed Y blend are not implemented in this prototype.

| Compared tick | G variance ratio | X variance ratio | Y variance ratio |
|---|---:|---:|---:|
| 0.0001 | 4.95263e-8 | 0.00558865 | 0.00981796 |
| 0.00005 | 2.45759e-11 | 0.0000633074 | 0.0000317990 |

The ratio is sample variance of the coupled tick difference divided by sample variance of the independent tick difference. Small ratios are encouraging, but four observations do not support a reliable gain estimate or a formal benefit verdict. In particular, the exceptionally small G ratios should not be extrapolated to the full time horizon. Coupling reduces sampling noise in a difference; it does not remove a genuine difference caused by tick size.

| Compared tick | Mean seconds per coupled pair | Mean seconds per independent pair | Cost ratio |
|---|---:|---:|---:|
| 0.0001 | 3.06636 | 2.19126 | 1.39936 |
| 5e-05 | 3.11934 | 2.19401 | 1.42175 |

Phase B took 43.1083 seconds including retention and audits, below its separate 60-second ceiling. Total bounded execution was 59.4909 seconds (0.0165252 core-hours), using one local worker and one thread per numerical library. Pairing costs more than two independent runs here, so a future study must evaluate variance reduction **per unit of compute** as well as raw variance reduction.

## Separate recomputation and retention

`src/check_mg12h3_tick_coupling.py` recomputed exact reference probabilities using a dense matrix exponential of the saved identity generator, checked the transport law using a two-state transition-matrix power, independently aggregated channel changes to reconstruct 16,400 chain endpoints, and replayed ten retained event/boundary sentinels. It recomputed actual contrasts and variance ratios from saved counts. All checks passed. This is a separate analysis path by the same author, sharing frozen rule and projection definitions; it is not external replication.

The archive contains every sampled endpoint and seed, exact augmented states/generator/probabilities, branch histograms, selected complete joint event/boundary records, source hashes, environment, attempts and budget ledger. Full event traces are retained for preselected sentinels rather than every replicate, consistent with the approved compact storage choice. No scientific retries or post-result threshold changes occurred. A pre-run local documentation encoding error was corrected before freezing and involved no stochastic test execution.

## Recommended next step

Prepare a prospective qualification packet covering mixed reactions with unequal tick calendars, longer numerical histories, and fresh actual-ring pairs at the intended 0.01 horizon. Evaluate the fixed direct contrasts first, before adding optional compensator/blend calculations to the joint engine. Register the variance/cost decision, denominators, confidence method, storage and measured budget in advance. The present microtest is insufficient to size a final 5% contrast-equivalence study from its four-pair variance estimates alone. No further stochastic batch has been launched.

Raw decisions: `stage_A_results.json`, `stage_B_results.json`, `separate_recomputation.json`, `budget_ledger.json`, and `plan.json` in this directory. Specification: `spec/MG12H3_tick_coupling_development.md`.
