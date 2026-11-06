# MG10-S1 — Finite-reactant stochastic sampler qualification

MG046, 2026-11-06. **Prospective design; no scientific ensemble has run.**
Richard signed off MG10-C1 and instructed proceeding. This prerequisite to MG10-D
qualifies an explicit fluctuation law, not a full-lattice particle survey. MG10-C1
is CLOSED with its three outcomes; original MG10 remains UNRESOLVED, not rescored.
No external review was requested with this execution instruction.

## Model and bounded approximation

Use the existing seven-species P17 continuous-time Markov jump model. State order
G,X,Y,A,B,Z,W; cumulative reaction order R1e,R2,R2r,R3e,R4,R5e,R1re. R1re is zero
at P17. Active neighbor hops have the old D/h^2 rates. A0/B0 are fixed positive
reference counts; fixtures may start after prescribed depletion. Outputs are kept.
This explicitly specifies classical reaction/hopping shot noise. It is not a
fundamental-tick law, a physical derivation of the reservoirs, or entanglement.

For a leap delta, all hazards use one start-state snapshot. No arriving product
reacts again in that leap. Count-based R4 catalysis follows the old convention:
its hazard is nX*(nX-1)*nY/Omega^2, with no new permanent-token reservation rule.

A conversion is Binomial(nA,Omega*a*delta/A0). G and Y use bounded multinomial
exits for reactions and hops. Joint R3 construction at each site:

    lambda3 = b*nX*nB/B0; m=min(nX,nB)
    K3 ~ Binomial(m,lambda3*delta/m), or zero if m=0
    p3X=lambda3*delta/nX, or zero if nX=0

Draw the remaining X exits from nX-K3 by a multinomial with probabilities
(g,1-g,Dx/h^2,Dx/h^2)*delta/(1-p3X). This guarantees K3<=nB and K3<=nX,
without clipping or rejection. E[K3]=lambda3*delta; because E[nX-K3]=nX*(1-p3X),
each other X exit also has its required Euler conditional mean. The finite-delta
joint distribution is an approximation and must be checked against exact SSA.

Choose delta<=requested maximum and the next checkpoint remainder, and enforce
delta*max(rateA,rateGtotal,rateYtotal,rateXtotal,lambda3/min(nX,nB))<=.2 over sites.
Skip the last term if its denominator is zero. Thus p3X<=.2 and conditional other
X probability sum<=.25. No probability clipping, negative-count repair or ignored
overdraw is permitted. Smaller adaptive steps are recorded, not hidden. Exact SSA
is independently expressed and applies one reaction/hop at a time in continuous time.

## Fixed small-system cases and ensembles

Checkpoints .05,.10,.25; h=1. Species arrays below are per node. Missing output
counts are zero. Diffusion is zero in the one-node cases; (1,1,12) in transport.

| Case | N, Omega | G; X; Y | A; B | Fixed A0; B0 | Z; W |
|---|---|---|---|---|---|
| closed | 1,1 | 16;16;2 | 40;40 | 40;40 | 0;0 |
| scarce_B | 1,1 | 2;12;3 | 8;1 | 8;1 | 0;0 |
| depleted | 1,1 | 2;8;2 | 0;0 | 8;8 | 8;8 |
| transport | 3,2 | [10,14,18];[12,20,28];[1,2,4] | [20,40,10];[1,8,40] | [40,40,40];[40,40,40] | zeros;zeros |

The depleted/transport states are specified predepleted fixtures with fixed
denominators, not reset normalization or newly simulated preparation histories.
Small reference stocks stress the boundary; no full-lattice capacity-limit claim.

For every case: exact SSA and three leap maxima 1/1000,1/2000,1/4000. Sixteen fixed
groups, each R=32768 trajectories in 64 chunks of512; no replacements or extensions.
Warm-up/timing cases use disjoint development keys and do not enter these ensembles.
Use runner.seed_for(10,case_index,method_index,stage,chunk_index), entropy20261027.
Method indices SSA0,coarse1,middle2,fine3; stages4601 development,4602 timing,4603
production. Check all production chunk seeds are distinct before registration;
if not, stop rather than silently selecting new streams. Work scheduling cannot
change stream keys. Replicates consume successive PRNG subsequences inside a fixed
chunk. Inference assumes these pseudorandom draws model independent trajectories.

## Frozen distribution comparison

At each checkpoint retain total X, total B, total W, total R1/R3/R4 counts, and
node0 X/Y. Compare their empirical CDFs with the independently seeded SSA group.
Primary comparison is the finest leap versus SSA: 4 cases*3 times*8 observables
=96 marginal comparisons. Coarse/middle distances are descriptive refinement
evidence, not extra equivalence decisions or a fitted convergence-order claim.

For each comparison let d_hat be the exact supremum empirical CDF difference.
At family error alpha=.01, set e=sqrt(log(4*96/alpha)/(2*32768)). The two-sample
DKW/union-bound band is [max(0,d_hat-2e),min(1,d_hat+2e)] simultaneously for all
96 distances under the stated independent-sampling assumption. Dependence among
observables and re-use of an SSA sample is allowed by the union bound.
The prespecified equivalence margin is .05 CDF units, not a p-value threshold.

- QUALIFIED ON TESTED SMALL SYSTEMS: all required exact-case/ledger checks pass,
  all16 groups complete at R, and every primary upper bound<=.05.
- NOT QUALIFIED ON TESTED CASES: complete data/valid implementation, at least one
  primary lower bound>.05. This is approximation disagreement, not a physical result.
- UNRESOLVED: an interval straddles the margin, missing data, timeout or budget stop.
- QUALIFICATION FAILED: invalid transitions/probabilities/counts, ledger or exact
  oracle failure; stop broader production and retain evidence.

This qualifies only these marginals/cases/steps/horizon. It is not total-variation
equivalence of the full joint law, universal convergence, or a particle claim.

## Qualification, cost and retention

Before the scientific ensembles: exact integer reaction/hop ledgers; zero stock;
scarce B never overdrawn; analytic single-channel A-death SSA/discrete-leap laws;
one-step joint R3 mean/law and all-channel conditional drift; unchanged IDEAL
regression; invalid/missing/interruption/budget labeling and controller checks;
independently recomputed outputs/statistics. Frozen IDEAL code is not changed and
new CHEMOSTAT byte-for-byte sampling equivalence is not claimed.

Lead-allocated local ceilings under the instruction to proceed: qualification/
timing360 worker-seconds; scientific ensembles3600 worker-seconds (1 core-hour),
per job600 seconds, <=8 workers with library threads1. Poll<=.5 seconds including
in-flight elapsed time; terminate only owned jobs. No automatic extension.
Measure all16 classes with512 timing trajectories each, then conservatively
project 2x compute scaling to32768 plus measured startup/checkpoint overhead.
Require total projection<=2700 seconds (25% reserve) before registration. If not,
stop and report a concrete smaller design proposal without silently changing R.

Commit specification, implementation, tests, seed allocation, forecast, durable
write/read check and registration before the16 groups. Retain every trajectory's
three snapshots of seven populations, seven reaction counts and six directed-hop
counts, with exact reconstruction checks, source/ref hashes, environment, all
chunk seeds, costs and failures. Save available chunk checkpoints on interruption.
At N3 each raw group is 32768*3*20*3*8=47185920 bytes before compression, below
100MiB. No full event history is required to recompute these decisions; SSA paths
can be replayed from fixed chunk seeds. Keep raw records in persistent local
results/MG10S1/ and commit/push all decision inputs. A separate analysis must read
raw arrays without importing the new sampler/evaluator, declaring shared helpers.
PI closure is separate. Full-lattice particle runs require a later registered
scope with persistence/identity criteria and a physically motivated fluctuation model.
