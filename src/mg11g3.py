"""MG092 B-only transport intervention; frozen G2 reaction and A-network laws."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,time,hashlib
import numpy as np
from scipy.integrate import BDF
from scipy.sparse import coo_matrix
from mg11g2 import rhs as old_rhs,initial
from mg11g2_pattern_v2 import environment,sparsity,solve_pattern,SETTINGS,S0,OBS as FULL_OBS,metrics as old_metrics
from mg11g2_pattern_v2_oracle import population_rhs,laplacians
ROOT=Path(__file__).resolve().parents[1];DONORS=ROOT/'results/MG11G2/pattern_response_v2/production';OBS=np.arange(61)/4
ORIGINS=('LOCAL','CROSS_RING');ARMS=('LOCAL','CROSS_RING','OFF');SWITCH_TIME=10.;PARENT_INDEX=40

def jobs():return [dict(id=f'o{o}_s{s}_v{v}_a{a}',origin=o,start=s,setting=v,arm=a) for o in (0,1) for s in (0,1) for v in (0,1) for a in (0,1,2)]
def donor_path(job):return DONORS/f'c{2+job["origin"]}_s{job["start"]}_v{job["setting"]}.npz'
def factory(job):
    if job not in jobs():raise ValueError('Noncanonical job')
    with np.load(donor_path(job)) as raw:
        assert np.array_equal(raw['times'],FULL_OBS);pop=raw['state'][PARENT_INDEX,:3500].reshape(7,500).copy();case=dict(initial=pop,a0=raw['a0'].copy(),b0=raw['b0'].copy(),omega=float(raw['omega']),diff=raw['diff'].copy())
    ga=environment(ORIGINS[job['origin']]);gb=environment(ORIGINS[job['arm']] if job['arm']<2 else ORIGINS[job['origin']]);return case,ga,gb,np.array([1/3,0. if job['arm']==2 else 1/3])

def rhs(t,y,case,ga,gb,k):
    n=ga['n'];e=len(ga['src']);out=old_rhs(t,y,case,ga,[k[0],0.],'FINITE');flow=k[1]*gb['weight']*y[4*n+gb['src']]
    np.add.at(out,4*n+gb['src'],-flow);np.add.at(out,4*n+gb['dst'],flow);out[20*n+e:20*n+2*e]=flow;return out

def pattern(ga,gb):
    base=sparsity(ga);n=ga['n'];e=len(ga['src']);rows=np.r_[4*n+gb['src'],4*n+gb['dst'],20*n+e+np.arange(e)];cols=np.tile(4*n+gb['src'],3)
    return (base+coo_matrix((np.ones(len(rows),dtype=bool),(rows,cols)),shape=base.shape).tocsr()).astype(bool)

def solve(job,times,deadline,checkpoint=None):
    case,ga,gb,k=factory(job);times=np.asarray(times);n=ga['n'];e=len(ga['src']);assert n==gb['n'] and e==len(gb['src'])
    if job['start']==0:
        def remap(y):
            z=y.copy();z[...,20*n+e:20*n+2*e]*=gb['weight']/ga['weight'];return z
        def save(t,y):
            if checkpoint:checkpoint(t,remap(y))
        result=solve_pattern(case,ga,k,'FINITE',times,job['setting'],deadline,checkpoint=save);result['state']=remap(result['state']);return result
    def fun(t,y):
        if time.perf_counter()>deadline:raise TimeoutError('G3 solve deadline')
        return rhs(t,y,case,ga,gb,k)
    solver=BDF(fun,0.,initial(case,ga),float(times[-1]),jac_sparsity=pattern(ga,gb),**SETTINGS[job['setting']]);saved=[initial(case,ga)];index=1
    while solver.status=='running':
        solver.step()
        if solver.status=='failed':raise RuntimeError('G3 BDF failed')
        if index<len(times) and times[index]<=solver.t:
            dense=solver.dense_output()
            while index<len(times) and times[index]<=solver.t:
                saved.append(dense(times[index]));index+=1
                if checkpoint:checkpoint(times[:index],np.array(saved))
    if len(saved)!=len(times) or not np.all(np.isfinite(saved)):raise RuntimeError('Incomplete/nonfinite G3 record')
    return dict(times=times,state=np.array(saved),homogeneous_reduction=False,nfev=solver.nfev,njev=solver.njev)

def reference_rhs(t,y,case,ga,gb,k):
    out=population_rhs(t,y,case,ga,[k[0],0.],'FINITE').reshape(7,ga['n']);L,_=laplacians(gb);out[4]+=k[1]*(L@y.reshape(7,ga['n'])[4]);return out.ravel()

def combined(job,pop,times):
    with np.load(donor_path(job)) as raw:past=raw['state'][:PARENT_INDEX+1,:3500].reshape(PARENT_INDEX+1,7,500)
    return np.concatenate((past,pop[1:])),np.r_[FULL_OBS[:PARENT_INDEX+1],SWITCH_TIME+times[1:]]

def measure(job,pop,times,full):
    z,t=combined(job,pop,times);return old_metrics(z,t,full=full)

def donor_manifest():
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted({donor_path(j) for j in jobs()})}
