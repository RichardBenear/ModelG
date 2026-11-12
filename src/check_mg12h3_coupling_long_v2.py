"""MG157 saved-only count/reservation recomputation and prospective analysis."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import json,math,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from mg12h3_coupling_long import BASE,read,hashes,DELTAS
from mg12h3_finer_precision import setup
from mg10e2c1_pool import atomic_json
def audit():
    start=time.perf_counter();p=read(BASE/'plan.json');assert hashes()==p['source_hashes'];assert read(BASE/'production_attempt.json')['status']=='COMPLETE';initial,contrast=setup();model=contrast.model;data={};endpoints=0;sentinel_events=0
    for j in p['jobs']:
        path=BASE/'production'/j['id'];meta=read(path.with_suffix('.json'));assert meta['status']=='COMPLETE';assert meta['source_hashes']==p['source_hashes'];previous=np.stack([initial,initial]);frames=[previous.copy()];last_rng=None;events=0
        for segment in range(10):
            with np.load(path/f'checkpoint_{segment:02d}.npz') as z:
                np.testing.assert_array_equal(z['initial'],previous);final=z['final'];hist=z['histogram'];events+=int(z['events'])
                if last_rng is not None:np.testing.assert_array_equal(z['rng_before'],last_rng)
                last_rng=z['rng_after'].copy()
                for side in range(2):
                    flat=previous[side].ravel().copy()
                    for i,ch in enumerate(model.channels):
                        for coordinate,v in ch.change:flat[coordinate]+=int(hist[side,i])*v
                    np.testing.assert_array_equal(flat,final[side].ravel());assert flat.min()>=0 and sum(map(int,flat))==sum(map(int,initial.ravel()));endpoints+=1
                if j['sentinel']:
                    with np.load(path/f'trace_{segment:02d}.npz') as archive:
                        t={key:archive[key] for key in archive.files}
                        available=previous.reshape(2,-1).copy();output=available.copy();recorded=np.zeros_like(hist);cursor=0;next_boundary=[round(d*1e9) for d in j['deltas']];steps=next_boundary.copy();last=-1.
                        for b,bt in enumerate(t['boundary_times']):
                            boundary=min(next_boundary);assert float(bt)==boundary/1e9;mask=sum(1<<s for s in range(2) if next_boundary[s]==boundary);assert mask==int(t['boundary_masks'][b])
                            while cursor<len(t['times']) and t['times'][cursor]<bt:
                                mark=float(t['times'][cursor]);assert last<mark;last=mark;index=int(t['channels'][cursor]);ch=model.channels[index];event_mask=int(t['masks'][cursor]);assert event_mask in (1,2,3)
                                for side in range(2):
                                    if event_mask&(1<<side):
                                        assert ch.rate(available[side])>0
                                        for coordinate,v in ch.reserve:assert available[side,coordinate]>=v;available[side,coordinate]-=v
                                        for coordinate,v in ch.change:output[side,coordinate]+=v
                                        recorded[side,index]+=1
                                cursor+=1
                            np.testing.assert_array_equal(available.reshape(2,7,500),t['boundary_available'][b]);np.testing.assert_array_equal(output.reshape(2,7,500),t['boundary_output'][b])
                            for side in range(2):
                                if mask&(1<<side):available[side]=output[side];next_boundary[side]+=steps[side]
                        assert cursor==len(t['times']);np.testing.assert_array_equal(recorded,hist);np.testing.assert_array_equal(output.reshape(2,7,500),final);sentinel_events+=cursor
                previous=final.copy();frames.append(previous.copy())
        assert events==meta['events']
        with np.load(path/'observations.npz') as z:
            np.testing.assert_array_equal(frames,z['counts']);np.testing.assert_array_equal(z['seeds'],j['seeds']);actual=np.asarray([[contrast.weights@x.ravel() for x in frame] for frame in z['counts']]);np.testing.assert_allclose(actual,z['contrasts'],rtol=1e-12,atol=1e-12)
        data[j['id']]=dict(contrasts=actual,seconds=meta['elapsed_seconds'])
    atomic_json(BASE/'raw_audit_v2.json',dict(passed=True,chain_checkpoint_endpoints=endpoints,sentinel_events=sentinel_events,source_hashes=hashes(),elapsed_seconds=time.perf_counter()-start,shared_dependencies='Frozen channel definitions, source/input hashes and fixed projection weights; same-author reconstruction, not external replication.'))
    rng=np.random.default_rng(p['seeds']['bootstrap']);rows=[];indices=[]
    for c in range(2):
        groups={kind:[data[f'c{c}_{kind}_r{r:02d}'] for r in range(8)] for kind in ('paired','independent')};arrays={k:np.asarray([x['contrasts'] for x in v]) for k,v in groups.items()};cost={k:np.asarray([x['seconds'] for x in v]) for k,v in groups.items()};differences={k:a[:,:,0]-a[:,:,1] for k,a in arrays.items()};vp=np.var(differences['paired'],axis=0,ddof=1);vi=np.var(differences['independent'],axis=0,ddof=1);ratio=vp[1:]/vi[1:];cost_ratio=cost['paired'].mean()/cost['independent'].mean();ip=rng.integers(0,8,(10000,8));ii=rng.integers(0,8,(10000,8));indices.append(np.stack([ip,ii]));dp=differences['paired'][:,-1,2];di=differences['independent'][:,-1,2];boots=np.var(dp[ip],axis=1,ddof=1)/np.var(di[ii],axis=1,ddof=1);cb=cost['paired'][ip].mean(axis=1)/cost['independent'][ii].mean(axis=1);adjusted=boots*cb;valid=bool(np.isfinite(adjusted).all());upper=float(np.quantile(adjusted,.975)) if valid else None;point=float(ratio[-1,2]);row=dict(comparison=c,deltas=DELTAS[c],variance_ratio=ratio.tolist(),mean_paired_difference=differences['paired'].mean(axis=0).tolist(),paired_difference_SE=(np.std(differences['paired'],axis=0,ddof=1)/math.sqrt(8)).tolist(),mean_pair_cost_seconds={k:float(v.mean()) for k,v in cost.items()},cost_ratio=float(cost_ratio),primary_Y_variance_ratio=point,primary_Y_cost_adjusted_ratio=point*cost_ratio,primary_Y_bootstrap_upper_97_5=upper,bootstrap_finite=valid,pilot_supported=valid and upper<1 and point<1);rows.append(row)
    np.savez_compressed(BASE/'bootstrap_indices.npz',indices=indices,seed=p['seeds']['bootstrap']);atomic_json(BASE/'analysis_results.json',dict(pilot_benefit_supported=all(r['pilot_supported'] for r in rows),comparisons=rows,bootstrap_replicates=10000,source_hashes=hashes(),limits='Approximate bootstrap coverage with eight pairs; no tick-equivalence or final variance qualification.',analysis_elapsed_seconds=time.perf_counter()-start));print(json.dumps(rows))
if __name__=='__main__':audit()
