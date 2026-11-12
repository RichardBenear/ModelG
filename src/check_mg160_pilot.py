"""MG160-P1 saved-only causal replay and separate FFT pattern recomputation."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import time,json,traceback,hashlib
import numpy as np
from numba import njit
from run_mg160_persistence import BASE,read,hashes
from mg12h3_finer_precision import setup
from mg10e2c1_pool import atomic_json
from mg160_local_archive import local_archive
BASE=local_archive(BASE)

@njit
def replay(times,masks,ids,bt,bm,ba,bo,initial,ri,rn,ci,cv,nchannels):
    available=initial.copy();output=initial.copy();hist=np.zeros((3,nchannels),np.int64);cursor=0;last=0.;previous_boundary=0.;next_tick=np.array([100000,20000],np.int64);steps=next_tick.copy()
    for boundary_index in range(len(bt)):
        ns=min(next_tick[0],next_tick[1]);boundary=ns/1e9;assert bt[boundary_index]==boundary
        expected=0
        for s in range(2):
            if next_tick[s]==ns:expected+=1<<s
        assert bm[boundary_index]==expected
        while cursor<len(times) and times[cursor]<boundary:
            mark=times[cursor];assert mark>last and mark>previous_boundary;channel=ids[cursor];mask=masks[cursor];assert 0<=channel<nchannels and 1<=mask<=3
            branch=0 if mask==3 else (1 if mask==1 else 2);hist[branch,channel]+=1
            for side in range(2):
                if mask&(1<<side):
                    for k in range(ri.shape[1]):
                        i=ri[channel,k]
                        if i>=0:assert available[side,i]>=rn[channel,k];available[side,i]-=rn[channel,k]
                    for k in range(ci.shape[1]):
                        i=ci[channel,k]
                        if i>=0:output[side,i]+=cv[channel,k];assert output[side,i]>=0
            last=mark;cursor+=1
        assert np.array_equal(available,ba[boundary_index]);assert np.array_equal(output,bo[boundary_index])
        for side in range(2):
            if expected&(1<<side):available[side]=output[side].copy();next_tick[side]+=steps[side]
        previous_boundary=boundary
    assert cursor==len(times) and previous_boundary==.001
    return output,hist

def recomputer(reference,canonical,omega):
    n=canonical.shape[1];k=2*np.pi*np.fft.fftfreq(n,d=.2);kernel=np.exp(-.5*(k*.4)**2)
    def smooth(x):return np.fft.ifft(np.fft.fft(x-x.mean())*kernel).real
    template=smooth(reference[1]);norm=float(template@template)
    def raw(pop):
        f=smooth(pop[1]/omega);energy=float(f@f)
        if energy<=1e-28:return dict(amplitude=0.,unaligned_amplitude=0.,correlation=0.,shift_sites=None,phase_ambiguous=False,participation_length=0.,residual_fraction=0.)
        scores=np.fft.ifft(np.fft.fft(f)*np.conj(np.fft.fft(template))).real;best=int(scores.argmax());peak=float(scores[best]);tied=np.flatnonzero(abs(scores-peak)<=1e-10*max(abs(peak),1e-20));distance=np.minimum((tied-best)%n,(best-tied)%n);amplitude=peak/norm
        return dict(amplitude=amplitude,unaligned_amplitude=float(f@template/norm),correlation=peak/np.sqrt(norm*energy),shift_sites=int((best+n//2)%n-n//2),phase_ambiguous=bool(np.any(distance>1)),participation_length=float(.2*energy**2/np.sum(f**4)),residual_fraction=float(np.linalg.norm(f-amplitude*np.roll(template,best))/np.sqrt(energy)))
    baseline=raw(canonical)
    def measure(pop):
        m=raw(pop);m['relative_amplitude']=m['amplitude']/baseline['amplitude'];m['relative_width']=m['participation_length']/baseline['participation_length'];m['coherent']=bool(m['relative_amplitude']>=.5 and m['correlation']>=.8);m['shape_retained']=bool(.5<=m['relative_width']<=2.);m['pattern_retained']=m['coherent'] and m['shape_retained'];return m
    return measure

def outcome(times,metrics,key):
    active=[bool(m[key]) for m in metrics];bad=[i for i,v in enumerate(active) if not v];first_on=[i for i,v in enumerate(active) if v]
    return dict(initial_present=active[0],retained_at_all_frames=bool(all(active)),positive_frames=sum(active),first_loss_interval=[float(times[bad[0]-1]),float(times[bad[0]])] if active[0] and bad else None,right_censored_at=float(times[-1]) if all(active) else None,first_positive_time=float(times[first_on[0]]) if first_on else None)

def check():
    start=time.perf_counter();attempt=BASE/'P1_analysis_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'))
    try:
        plan=read(BASE/'P1_plan.json');assert hashes()==plan['source_hashes'] and read(BASE/'P1_attempt.json')['status']=='COMPLETE';ledger=read(BASE/'P1_production/ledger.json');assert len(ledger['records'])==8 and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records']);canonical,contrast=setup();model=contrast.model;nc=len(model.channels);ri=np.full((nc,2),-1,np.int64);rn=np.zeros_like(ri);ci=np.full((nc,4),-1,np.int64);cv=np.zeros_like(ci)
        for i,ch in enumerate(model.channels):
            assert ch.coefficient>0 and ch.reserve==ch.factors and sum(v for _,v in ch.change)==0
            for k,(coord,num) in enumerate(ch.reserve):ri[i,k]=coord;rn[i,k]=num
            for k,(coord,v) in enumerate(ch.change):ci[i,k]=coord;cv[i,k]=v
        with np.load(BASE/'prepared_inputs.npz') as z:prepared=z['initial'];homogeneous=z['homogeneous'];reference=z['reference'];omega=float(z['omega'])
        measure=recomputer(reference,prepared,omega);rows=[];endpoints=0;sentinel_events=0;inputs={}
        for job in plan['jobs']:
            path=BASE/'P1_production'/job['id'];meta=read(path.with_suffix('.json'));assert meta['source_hashes']==plan['source_hashes'];initial=prepared if job['condition']=='prepared' else homogeneous;previous=np.stack([initial,initial]);frames=[previous.copy()];last_rng=json.dumps(np.random.default_rng(job['seed']).bit_generator.state);total=0;reaction=np.zeros((2,6),np.int64)
            for segment in range(100):
                checkpoint=path/f'checkpoint_{segment:03d}.npz'
                with np.load(checkpoint) as z:record={k:z[k] for k in z.files}
                np.testing.assert_array_equal(record['initial'],previous);assert str(record['rng_before'])==last_rng;last_rng=str(record['rng_after']);final=record['final'];hist=record['histogram'];assert hist.shape==(3,nc) and hist.min()>=0 and int(hist.sum())==int(record['events']);total+=int(record['events'])
                for side in range(2):
                    flat=previous[side].ravel().copy();accepted=hist[0]+hist[side+1]
                    for i,ch in enumerate(model.channels):
                        for coordinate,v in ch.change:flat[coordinate]+=int(accepted[i])*v
                        if ch.rule in ('R1','R2','R2r','R3','R4','R5'):reaction[side,('R1','R2','R2r','R3','R4','R5').index(ch.rule)]+=accepted[i]
                    np.testing.assert_array_equal(flat,final[side].ravel());assert flat.min()>=0 and sum(map(int,flat))==sum(map(int,initial.ravel()));endpoints+=1
                if job['sentinel']:
                    with np.load(path/f'trace_{segment:03d}.npz') as archive:t={k:archive[k] for k in archive.files}
                    out,replayed=replay(np.asarray(t['times'],float),np.asarray(t['masks'],np.int64),np.asarray(t['channels'],np.int64),t['boundary_times'],t['boundary_masks'],t['boundary_available'].reshape(-1,2,3500),t['boundary_output'].reshape(-1,2,3500),previous.reshape(2,-1),ri,rn,ci,cv,nc);np.testing.assert_array_equal(out,final.reshape(2,-1));np.testing.assert_array_equal(replayed,hist);np.testing.assert_array_equal(t['boundary_available'][-1],record['available']);sentinel_events+=len(t['times'])
                previous=final.copy();frames.append(previous.copy())
            assert total==meta['events'];observations=path/'observations.npz';inputs[observations.relative_to(BASE).as_posix()]=hashlib.sha256(observations.read_bytes()).hexdigest()
            with np.load(observations) as z:np.testing.assert_array_equal(z['counts'],frames);assert int(z['seed'])==job['seed'];times=z['times'];np.testing.assert_array_equal(times,np.arange(101)*.001)
            stored=read(path/'metrics.json')['values'];metrics=[[measure(pop) for pop in frame] for frame in frames];trajectories=[]
            for side in range(2):
                trajectory=[f[side] for f in metrics]
                for f,m in enumerate(trajectory):
                    original=stored[f][side]
                    for key in ('amplitude','unaligned_amplitude','correlation','participation_length','residual_fraction','relative_amplitude','relative_width'):np.testing.assert_allclose(m[key],original[key],rtol=1e-9,atol=1e-11)
                    for key in ('phase_ambiguous','coherent','shape_retained','pattern_retained'):assert m[key]==original[key]
                    if not m['phase_ambiguous']:assert m['shift_sites']==original['shift_sites']
                fields=np.array(frames)[:,side];stocks=fields.sum(axis=2);assert np.all(np.diff(stocks[:,3:5],axis=0)<=0) and np.all(np.diff(stocks[:,5:7],axis=0)>=0);delta=stocks[-1]-stocks[0];assert -delta[3]==reaction[side,0] and -delta[4]==delta[5]==reaction[side,3] and delta[6]==reaction[side,5]
                trajectories.append(dict(side=side,delta=plan['deltas'][side],coherence=outcome(times,trajectory,'coherent'),shape=outcome(times,trajectory,'shape_retained'),pattern=outcome(times,trajectory,'pattern_retained'),initial=trajectory[0],final=trajectory[-1],minimum_relative_amplitude=min(m['relative_amplitude'] for m in trajectory),minimum_correlation=min(m['correlation'] for m in trajectory),relative_width_range=[min(m['relative_width'] for m in trajectory),max(m['relative_width'] for m in trajectory)],maximum_abs_wrapped_shift_sites=max(abs(m['shift_sites']) for m in trajectory if m['shift_sites'] is not None),ambiguous_frames=sum(m['phase_ambiguous'] for m in trajectory),species_total_change=delta.tolist()))
            rows.append(dict(id=job['id'],condition=job['condition'],seed=job['seed'],trajectories=trajectories));atomic_json(BASE/'P1_analysis_progress.json',dict(jobs_complete=len(rows),required=8,sentinel_events=sentinel_events))
        atomic_json(BASE/'P1_analysis_results.json',dict(passed=True,jobs=rows,source_hashes=hashes(),input_sha256=inputs,checkpoint_chain_endpoints=endpoints,sentinel_events=sentinel_events,method='Independent FFT smoothing/cross-correlation and extent/outcome recomputation; compiled integer event/reservation replay and separate histogram aggregation. Shared frozen channel/inputs and same author; not external replication.',elapsed_seconds=time.perf_counter()-start));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start));print(json.dumps([dict(id=r['id'],outcomes=[x['pattern'] for x in r['trajectories']]) for r in rows]))
    except BaseException:
        atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':check()
