# MG11-G2 Step1 — coupling qualification results (MG087/MG087a)

Step1 built and checked the connection between the environmental network and the full nonlinear Model G continuum reactions on fixed small systems. Explicit ledgers reconstruct resource consumption, transport, outputs and maintained-reservoir imports; the established equations are recovered when environmental transport is off or supplies are maintained. Increasing finite capacity approached the maintained control, and sufficiently rapid mixing approached a shared finite bath. The original mixing range fell slightly short of its fixed target; a separately recorded stronger-mixing test met that same target. This prepares the next study of pattern response, while leaving particle persistence, stochastic accuracy, emergent vacuum structure and quantum connections untested.

**Build complete; coupling qualified on the tested small deterministic systems, with the original mixing-range limitation retained.** This is Step1 implementation qualification only. No pattern-response or stochastic production study has run. No original registered result is rescored, and no tolerances were relaxed. PI acceptance of this combined record remains separate.

## What was built

- Reused frozen N1 continuum reaction/transport law; historical sources unchanged.
- Integrated all seven populations, seven reaction fluxes, six directed active-hop fluxes, both resources' directed environmental-hop fluxes, and explicit maintained-reservoir imports.
- FINITE mode conserves total tracked concentration over equal-volume sites. CHEMOSTAT holds uniform reference stocks with an explicit import ledger; imports do not become hidden supply.
- Added the appropriate shared FINITE-bath nonlinear reference. Fast mixing redistributes resources and does not replenish them.
- Added bounded deterministic qualification, frozen configurations/source hashes, failure retention and a separate equation/incidence/Radau recomputation path. No stochastic kernel, random sampling or fundamental-tick rule was added.

## Fixed tests and outcomes

| Test | Outcome |
|---|---|
| Local derivative reconstruction, closed/import balances | Pass |
| Zero environmental transport versus established nonlinear drift | Pass |
| Maintained supplies versus established active equations | Pass |
| Zero A/B consumption, empty-species boundary signs | Pass |
| Five invalid mode/configuration/grid cases | Rejected as required |
| Nine original coupled trajectories versus separate Radau equations | Pass numerical/accounting checks |
| Capacity recovery at C=1,10,100 | Pass unchanged .001 final active-error target |
| Original rapid-mixing grid 1,10,100,1000 | **Failed** unchanged .005 final-error target |
| Separate additional mixing point 10000 | **Passed** unchanged .005 target and same numerical/accounting checks |

Original primary fast-mixing errors were [0.7065595021890594, 0.3400475891166124, 0.050157394972244876, 0.005213697869783119]. They decreased at every step, but the highest registered rate gave **.005213697869783119 > .005**. The original driver stopped at that assertion. All original code, plan, primary arrays, failure and criterion remain frozen at b8fea77; the separate audit reproduced the failure rather than discarding it. Original qualification remains failed under its original range.

Before integrating the extra rate, MG087a froze kappaA=kappaB=10000, the same fixture/horizon/solvers/scales and unchanged .005 criterion, with no further automatic extension if it failed. New error **0.000524351935444**; separately computed error **0.000524347584884**. Source/plan/spec/hash snapshot published at b8fea77 before follow-up integrations. This is one explicitly disclosed additional qualification point, not a retry of the original trajectories or retrospective selection of a favorable tolerance.

Capacity active errors were [0.0033353285290649505, 0.00033407797556126914, 3.341336519742936e-05], strictly decreasing and final **3.34133651974e-05 < .001**. Stocks/reference denominators increased together while initial activities remained1; Omega was fixed. This is recovery over T=1 for a declared nonuniform field, not long-horizon pattern recovery.

Across all ten coupled trajectories: maximum separately recomputed normalized population discrepancy **2.10711192494e-07 <= 2e-6**; ledger reconstruction error **8.70636895911e-13 <= 1e-8**; normalized closed/import balance residual **8.73170796633e-15 <= 1e-8**. Minimum recorded component **0**, without clipping. Original audit passed all nine raw/source/config/time/ledger checks and the added trajectory passed the same numerical/accounting criteria.

## Scope, execution and retention

Mean-field concentrations use x^2*y; they are not finite-count SSA drift with X(X-1)Y/Omega^2. Active G/X/Y stay on a fixed periodic ring; environmental edges are separate and fixed. Equal-volume symmetric graphs and local contacts remain assumptions. Maintained-reference stocks are uniform; no heterogeneous clamped-bath withdrawals are silently added.

The heterogeneous-stock mixing comparison starts at t=.25 (then .5,1), after the initial mixing layer. Resource states at t=0 are not asserted to equal a spatially mixed bath. A single small nonlinear example is not a general fast-mixing theorem, and no finite rate is universally sufficient. The kappa1000 failure is a useful warning against assuming the linear qualification scale automatically suffices for nonlinear demand.

Original coupled integrations finished before the failed criterion; no original trajectory rerun. Subsequent original separate audit took 4.649994seconds and the one-point extension/audit took 0.685080seconds, each below its prospective300-second ceiling. One local process, numerical-library threads1, zero random draws. Same-author separate equations, directed incidence and Radau solver share fixture constants/NumPy/SciPy; this is not external replication. Deadline guards and failure metadata remain in the tools. Original attempt.json is a retained startup record; failure.json and original_audit.json identify its terminal failed qualification.

All original/follow-up population, reaction, directed-hop and import arrays, independent references, configurations, plans, source hashes and failure records are published under results/MG11G2. raw_manifest.json and retention_verified.json identify and verify the binary decision inputs. No pattern/stochastic run is authorized by this implementation work.

## Recommended next use

Prepare a small deterministic pattern-response study with the active lattice and reaction laws fixed, matched environmental stocks/rates/contact roles, maintained/isolated-environment controls and explicit depletion measurements. Its graph choices, horizon, observables, numerical refinements and measured local budget should be specified and approved before production. Choose mixing scales based on the relevant nonlinear demand, not the G1 probe scale. Do not treat a fixed-template classification loss alone as particle death, given the earlier MG10 background-response findings.
