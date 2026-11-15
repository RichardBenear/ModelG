"""MG184 Phase 2 long-horizon adapter; frozen kernel and detector unchanged."""
import numpy as np
from runner import reference
from lattice import P17,steady_state,tau_leap_euler_1d
from measure import analyze_run,analyze_homogeneous,fourier_shift
H=.2;TAU=6.666666666666667e-4;B=29.;L=100.;OMEGAS=(100,300,1000,3000,10000)
def ref():return reference(B,H,L)
def initial(omega,condition,seed,reference_data=None):
 r=ref() if reference_data is None else reference_data
 if condition=='prepared':return np.rint(omega*(r['S0'][:,None]+r['phi'])).astype(np.int64)
 if condition=='homogeneous':return np.random.default_rng(seed).poisson(omega*r['S0'][:,None],size=(3,len(r['x']))).astype(np.int64)
 raise ValueError('Unknown condition')
def advance(pop,omega,duration,start,seed,tau=TAU):
 r=ref();steps=round(duration/tau);rec=round(.5/tau);assert abs(steps*tau-duration)<1e-9 and abs(rec*tau-.5)<1e-12
 a=[v.copy() for v in pop];out=tau_leap_euler_1d(*a,float(omega),P17['a'],B,P17['g'],P17['p'],P17['q'],P17['dx'],P17['dy'],H,tau,steps,float(start),0.,r['x'],rec,int(seed));fields=np.stack(out[:3],axis=1);times=out[3];return dict(counts=fields,times=times,clips=int(out[4]),final=np.array(a),steps=steps)
def summarize(counts,times,omega,condition,horizon=1000.,reference_data=None):
 r=ref() if reference_data is None else reference_data;counts=np.asarray(counts);times=np.asarray(times);assert counts.shape==(len(times),3,len(r['x'])) and abs(times[0])<1e-9 and abs(times[-1]-horizon)<1e-7 and np.all(np.diff(times)>0)
 if condition=='prepared':
  legacy=analyze_run(counts[:,0],counts[:,1],counts[:,2],times,omega,H,r['S0'],r['A_s'],r['lam'],T=horizon);positions=np.array(legacy['positions']);valid=np.isfinite(positions)&(np.array(legacy['n_dets'])==1);historical=bool(legacy['eligible']);strict=bool(np.all(valid))
 else:
  legacy=analyze_homogeneous(counts[:,0],counts[:,1],counts[:,2],times,omega,H,r['S0'],r['A_s'],r['lam']);positions=np.full(len(times),np.nan);valid=np.zeros(len(times),bool);historical=None;strict=None
 field=counts/omega;template=r['phi']-r['phi'].mean(axis=1,keepdims=True);norm=np.linalg.norm(template,axis=1);unit=template/norm[:,None];centered=field-field.mean(axis=2,keepdims=True);fixed=np.einsum('tsn,sn->ts',centered,unit);aligned=np.full((len(times),3),np.nan);residual=aligned.copy()
 for i in np.where(valid)[0]:
  v=fourier_shift(centered[i],H,positions[i]);q=np.einsum('sn,sn->s',v,unit);aligned[i]=q;residual[i]=np.linalg.norm(v-q[:,None]*unit,axis=1)/norm
 blocks=[]
 for i in range(10):
  mask=(times>=horizon*i/10)&(times<(horizon*(i+1)/10) if i<9 else times<=horizon);keep=mask&valid;blocks.append(dict(left=horizon*i/10,right=horizon*(i+1)/10,snapshots=int(mask.sum()),unique_original_snapshots=int(keep.sum()),aligned_mean=(aligned[keep].mean(axis=0).tolist() if keep.any() else [None]*3),fixed_mean=fixed[mask].mean(axis=0).tolist()))
 last_unmatched=bool(condition=='prepared' and not valid[-1]);return dict(horizon=horizon,condition=condition,historical_X_persistence=historical,strict_sampled_unique_X_track=strict,end_unmatched_or_ambiguous=last_unmatched,legacy=legacy,blocks=blocks,full_GXY_persistence='NOT_SCORED_FROM_X_ALONE'),dict(times=times,fixed_projection=fixed,aligned_projection=aligned,aligned_relative_amplitude=aligned/norm,aligned_shape_residual=residual,unique_original=valid,positions=positions,reference_norms=norm)
