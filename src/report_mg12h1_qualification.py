"""MG108 report/figure only; does not run or rescore qualification."""
import json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare_mg12h1 import ROOT,read
BASE=ROOT/'results/MG12H1/qualification'
def main():
    e=read(BASE/'evaluation.json');run=read(BASE/'run_state.json');reg=read(BASE/'registration.json');ledger=read(BASE/'ledger.json');assert e['audit_passed'] and e['complete'];assert len(ledger['records'])==16 and all(x['status']=='COMPLETE' for x in ledger['records']);passed=e['qualified'];maximum=max(x['upper'] for x in e['formal']);observed=max(x['observed_TV'] for x in e['formal']);rows=[]
    labels=['Active token + reactions','Full local chemistry','Separate A/B networks','Feedback catalysts']
    for c,label in enumerate(labels):
        exact=max(x['upper'] for x in e['formal'] if x['cell']==c and x['method']==0);fine=max(x['upper'] for x in e['formal'] if x['cell']==c and x['method']==3);rows.append(f'| {label} | {exact:.9f} | {fine:.9f} | 0.10 |')
    fig,ax=plt.subplots(figsize=(9.5,4.7));pos=np.arange(4);width=.34
    for offset,method,label,color in [(-width/2,0,'Exact reference executor','#617e8b'),(width/2,3,'Strict ticks: maximum 0.005','#157a69')]:
        heights=[max(x['upper'] for x in e['formal'] if x['cell']==c and x['method']==method) for c in range(4)];bars=ax.bar(pos+offset,heights,width,label=label,color=color)
        ax.bar_label(bars,fmt='%.3f',padding=3,fontsize=9)
    ax.axhline(.1,color='#ad543f',linestyle='--',label='Registered limit 0.10');ax.set_xticks(pos,['Active token','Full chemistry','A/B networks','Feedback']);ax.set_ylabel('Largest joint-count TV upper bound\nacross times 0.05 and 0.10');ax.set_ylim(0,max(.125,maximum*1.2));ax.set_title('MG12-H1: registered small-system endpoint qualification');ax.legend(loc='upper left',fontsize=9);ax.grid(axis='y',alpha=.15);fig.tight_layout();fig.savefig(BASE/'MG12H1_qualification_summary.png',dpi=170);plt.close(fig)
    summary=('The explicit hypergraph engine passed its formal small-system qualification. Both its continuous-time executor and the strict-tick executor—with catalysts limited to one participation per tick—reproduced the tested joint species-count endpoint distributions within the registered error bound. All saved event histories and accounting checks passed. This establishes a controlled bridge between our count-based kinetics and explicit relational rewrites on the tested fixed substrates; it supports moving toward more general connectivity. It does not yet establish accurate particle-scale dynamics, evolving geometry, a physical clock, emergent spatial dimension or a quantum vacuum source.' if passed else 'The explicit hypergraph engine completed its formal small-system experiment and its saved accounting/causal audit, but did not meet every registered endpoint-distribution accuracy criterion. The failed bounds are retained without changing the test. This identifies a qualification limit to resolve before extending connectivity or testing particles; it does not by itself identify the source of the discrepancy or invalidate earlier closed stages.')
    text=f'''# MG12-H1 formal qualification results — MG108

{summary}

**Outcome: {e['verdict']}. PI scientific closure pending.**

All16jobs and262144scientific paths completed. No warmups entered the sample, and no failure, timeout, cap stop, retry, replacement or post-outcome criterion change occurred. Saved-data audit passed for all paths.

![Registered joint-count endpoint bounds](MG12H1_qualification_summary.png)

## What was tested

The frozen H0/H1 engine represents tokens and sites through typed incidence relations. Chemical events match local participants and replace their state relations; transport changes attachments. Its continuous-time executor re-enumerates after each event. Its separate strict-tick executor freezes tick-start matches, prioritizes exponential proposals, reserves all token participants including catalysts, rejects conflicts without replacement and commits disjoint rewrites together. New products and arrivals wait until another tick. These are distinct finite-time dynamics, not two names for the same algorithm.

The unchanged four fixtures contain at most6tokens and3sites: one active token with reactions and hopping; six-token full local chemistry; two A plus two B on distinct weighted networks; and two X plus two Y with feedback only. Their finite count generators have6,138,36and3reachable states, respectively. Topology stays fixed. Pairwise transport and a three-participant reaction do not establish any spatial dimension.

Each fixture/method used16384independent fixed seeded paths. Methods were exact continuous time and strict maximum ticks.02,.01,.005, observed at0,.05,.1. The finest.005was already designated in the development design. Only exact and finest were formally scored; coarse/medium results are descriptive and were not used to select a setting.

## Registered scientific result

The observable is the **joint species-count distribution across all sites** at each positive endpoint. It is stronger than the earlier arbitrary-state-index CDF diagnostic, but still does not specify the complete temporal or token-identity law. The categorical reference comes from the separately written count-generator master equation.

| Fixture | Exact: largest upper | Finest tick: largest upper | Registered limit |
|---|---:|---:|---:|
{chr(10).join(rows)}

Across16formal comparisons, largest observed total variation was **{observed:.12f}** and largest upper bound **{maximum:.12f}**. Supported comparisons: **{sum(x['supported'] for x in e['formal'])}/16**. Total variation measures the maximum probability difference for an event defined by the joint count state. The limit0.10 is a deliberately coarse ten-percentage-point endpoint-law bound, not a claim of ten-percent error in every reaction rate or concentration.

For each fixture's Sreachable states and R=16384paths, the registered upper is observedTV+sqrt((S ln2+ln(16/.001))/(2R)). A multinomial L1 concentration bound and union bound give **nominal99.9% simultaneous confidence** for the16specified bounds, conditional on fixed sampling/PRNG assumptions and numerical reference fidelity. Shared trajectories/endpoints do not require independence for the union bound. This confidence is not a probability that Model G physics is correct.

Both sparse and dense matrix-exponential reference calculations agreed within the frozen1e-11tolerance; the feedback-only probabilities also matched their analytic three-state chain within1e-11. Probability/row-sum checks passed. These support the computed reference but are not a rigorous floating-point error certificate. No clipping or normalization was used to conceal errors.

The coarser16comparisons are retained in evaluation.json for description. Their sampling distances need not decrease monotonically with tick refinement. Passing the finest endpoint bounds does not establish convergence at every population, horizon or substrate.

## Accounting and independent saved-data check

The sampling driver checked causal traces on every path. The saved-data checker then reread all262144paths and re-derived rule effects, catalyst reservations, conflicts, participant histories, closed identity/count totals and observation frames from literal saved relations. It replayed **{e['accepted_events_replayed']}accepted events**. All registered seeds, completed denominators, metadata/raw checksums, time grids, fixture definitions and source hashes matched. All path prefixes became complete; no partial jobs were scored as complete.

The checker recomputed categorical frequencies using both bincount and sorted/searchsorted arithmetic, and recomputed concentration margins through two equivalent formulas. It shares the mathematical constants, fixtures and source-hash helper; role replay and count-law analysis are separate from primary incidence matching/execution. This is same-author recomputation, not independent external replication. The master-equation oracle does not call the primary matcher/commit/executors.

## Cost, provenance and retention

Richard's instruction “go ahead with the formal qualification” authorized this bounded qualification. The exact protocol, source hashes, fixtures, all262144seeds and deterministic resolution of{len(reg['seed_collisions_resolved'])}uint32collisions were published before any qualification random draws or reference integrations. Earlier development kernels and data remain unchanged. No external review was requested.

One local worker/library-thread1 completed sampling, reference construction and automatic analysis/audit in **{run['elapsed_seconds']:.6f}seconds ({run['elapsed_seconds']/60:.2f}minutes; {run['elapsed_seconds']/3600:.6f}single-worker hours)**, below the900second cap. The prospective measured forecast was{reg['forecast_seconds']:.2f}seconds; neither figure is an estimate for a full-lattice explicit-token simulation. Finished **{run['finished']}** UTC. Reporting and repository publication are separate from that capped run interval.

All16count archives,16eligible-event trace JSONL files,4reference probability archives, reference checks, full seed/config registration, per-job completion/progress records, budget ledger, analysis and this report/figure are retained publicly in this directory. Refer to raw_manifest.json and retention_verified.json for hashes and post-publication verification. Source text hashes normalize LF where specified; working-file and Git text bytes can differ by line endings. No new run is authorized by this result; PI closure is separate.

## Next implication

The next bounded design can extend fixed connectivity beyond the ring and test explicitly higher-arity substrate relations, while retaining the same reaction/causal contract. Known2D/3Dfixtures can validate transport and dimension measurements; arbitrary relational substrates need not have a single dimension. Evolving-topology rules and substrate-generated supply remain later gates. Before particle-scale work, combinatorial matching costs require an optimized representation proven equivalent to the identity-rich reference, plus a separate large-population and longer-horizon qualification. Keep G5 closed with its unresolved stochastic response; this result does not rescore it.
'''
    replacements={'All16jobs and262144scientific paths':'All 16 jobs and 262,144 scientific paths','at most6tokens and3sites':'at most 6 tokens and 3 sites','have6,138,36and3reachable states':'have 6, 138, 36 and 3 reachable states','used16384independent':'used 16,384 independent','ticks.02,.01,.005':'ticks 0.02, 0.01, 0.005','at0,.05,.1':'at 0, 0.05, 0.1','finest.005was':'finest 0.005 was','Across16formal':'Across 16 formal','limit0.10':'limit 0.10','Sreachable':'S reachable','R=16384paths':'R=16,384 paths','nominal99.9%':'nominal 99.9%','the16specified':'the 16 specified','frozen1e-11tolerance':'frozen 1e-11 tolerance','coarser16comparisons':'coarser 16 comparisons','all262144paths':'all 262,144 paths','121036accepted events':'121,036 accepted events','all262144seeds':'all 262,144 seeds','of5uint32collisions':'of 5 uint32 collisions','worker/library-thread1':'worker with one numerical-library thread','304.760835seconds':'304.760835 seconds','5.08minutes':'5.08 minutes','0.084656single-worker hours':'0.084656 single-worker hours','the900second cap':'the 900-second cap','was416.12seconds':'was 416.12 seconds','All16count archives,16eligible-event trace JSONL files,4reference probability archives':'All 16 count archives, 16 eligible-event trace JSONL files and 4 reference probability archives','Known2D/3Dfixtures':'Known 2D/3D fixtures'}
    for old,new in replacements.items():text=text.replace(old,new)
    (BASE/'MG12H1_qualification_results.md').write_text(text,encoding='utf-8')
if __name__=='__main__':main()
