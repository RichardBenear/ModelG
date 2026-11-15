"""G-first summaries; historical X status never upgrades full GXY persistence."""
import json,math
import numpy as np
from pathlib import Path
from scipy.stats import beta,t
from run_mg185_phase2_G import PUBLIC,RAW,read
from mg10e2c1_pool import atomic_json

def cp(k,n):return [0 if k==0 else float(beta.ppf(.025,k,n-k+1)),1 if k==n else float(beta.ppf(.975,k+1,n-k))]
def run():
 plan=read(PUBLIC/'plan.json');records=read(RAW/'main'/'ledger.json')['records'];rows=[];G=[]
 for omega in (100,300,1000,3000,10000):
  prepared=[r for r in records if r['job']['omega']==omega and r['job']['condition']=='prepared'];complete=[r for r in prepared if r['status']=='COMPLETE'];surv=sum(r['historical_X_persistence'] for r in complete);strict=sum(r['strict_sampled_unique_X_track'] for r in complete);unknown=30-len(complete);rows.append(dict(omega=omega,planned=30,complete=len(complete),unresolved_computational=unknown,historical_X_survived=surv,strict_sampled_unique= strict,survival_CI95=cp(surv,30) if not unknown else None,survival_fraction_identification_bounds=[surv/30,(surv+unknown)/30],GXY_persistence='NOT_INFERRED_FROM_X'))
  for condition in('prepared','homogeneous'):
   summaries=[]
   for r in records:
    if r['job']['omega']==omega and r['job']['condition']==condition and r['status']=='COMPLETE':summaries.append(read(RAW/'main'/r['id']/'summary.json'))
   for block in range(10):
    for which in('fixed_mean','aligned_mean'):
     values=[s['blocks'][block][which][0] for s in summaries if s['blocks'][block][which][0] is not None];n=len(values);mean=float(np.mean(values)) if n else None;half=float(t.ppf(.975,n-1)*np.std(values,ddof=1)/np.sqrt(n)) if n>1 else None;G.append(dict(omega=omega,condition=condition,block=block,observable=which,contributors=n,planned=30 if condition=='prepared' else 6,mean_G_projection=mean,CI95=[mean-half,mean+half] if half is not None else None,conditional_on_unique_original=which=='aligned_mean',valid_snapshot_counts=[s['blocks'][block]['unique_original_snapshots'] for s in summaries],not_independent_time_frames=True))
 atomic_json(PUBLIC/'G_population_summary.json',dict(survival=rows,G_blocks=G,old_verdicts_unchanged=True,full_GXY_claim=False,conditional_G_missing_not_filled=True))
if __name__=='__main__':run()
