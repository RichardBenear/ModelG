"""MG115 incidence-derived population-limit bridge; not a strict-tick sampler."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.integrate import solve_ivp
from mg11g3 import factory as prepared,donor_path
from mg11g2_pattern_v2 import factory as original,metrics
ROOT=Path(__file__).resolve().parents[1]
SETTINGS=(dict(rtol=1e-9,atol=1e-11,max_step=.1),dict(rtol=1e-11,atol=1e-13,max_step=.05))
def jobs():
    finite=[dict(id=f'f_o{o}_s{s}_a{a}_v{v}',kind='FINITE',origin=o,start=s,arm=a,setting=v,duration=15.) for o in range(2) for s in range(2) for a in range(2) for v in range(2)]
    return finite+[dict(id=f'{kind}_s{s}_v{v}',kind=kind,start=s,setting=v,duration=25.) for kind in ('MAINTAINED','ISOLATED_C1') for s in range(2) for v in range(2)]
def factory(job):
    if job not in jobs():raise ValueError('Noncanonical job')
    if job['kind']=='FINITE':
        old=dict(id=f'o{job["origin"]}_s{job["start"]}_v1_a{job["arm"]}',origin=job['origin'],start=job['start'],setting=1,arm=job['arm']);case,ga,gb,k=prepared(old);reference=ROOT/'results/MG11G3/production'/(old['id']+'.npz');parent=donor_path(old);mode='FINITE'
    else:
        case,ga,k,mode=original(job['kind'],job['start']);gb=ga;reference=ROOT/f'results/MG11G2/pattern_response_v2/production/c{0 if mode=="CHEMOSTAT" else 1}_s{job["start"]}_v1.npz';parent=None
    return case,ga,gb,k,mode,reference,parent
def operator(n,cells):
    rows=[];cols=[];values=[]
    for sites,weight in cells:
        if len(set(sites))!=len(sites) or len(sites)<2 or weight<=0:raise ValueError('Invalid cell')
        for source in sites:
            for dest in sites:
                if source!=dest:rows.extend((dest,source));cols.extend((source,source));values.extend((weight/(len(sites)-1),-weight/(len(sites)-1)))
    return coo_matrix((values,(rows,cols)),shape=(n,n)).tocsr()
def graph_cells(g):return [((int(i),int(j)),float(w)) for i,j,w in zip(g['src'],g['dst'],g['weight']) if i<j]
def matrices(case,ga,gb):
    n=ga['n'];active=operator(n,[((i,(i+1)%n),1.) for i in range(n)]);a=operator(n,graph_cells(ga));b=operator(n,graph_cells(gb));return active,a,b
def rhs(y,case,mats,k,mode):
    n=case['initial'].shape[1];p=y[:7*n].reshape(7,n);g,x,yy,a,b=p[:5];source=14*a/case['a0'];conversion=29*x*b/case['b0'];feedback=x*x*yy;out=np.zeros((7,n));out[0]=source-g+.1*x;out[1]=g-.1*x-conversion+feedback-.9*x;out[2]=conversion-feedback;out[3]=-source;out[4]=-conversion;out[5]=conversion;out[6]=.9*x
    for s in range(3):out[s]+=case['diff'][s]*(mats[0]@p[s])
    out[3]+=k[0]*(mats[1]@p[3]);out[4]+=k[1]*(mats[2]@p[4]);imports=np.zeros((2,n))
    if mode=='CHEMOSTAT':imports=-out[3:5].copy();out[3:5]=0
    return np.r_[out.ravel(),imports.ravel()]
def sparsity(n,mats):
    rows=[];cols=[]
    for i in range(n):
        for s in range(9):
            for c in range(5):rows.append(s*n+i);cols.append(c*n+i)
    for species,matrix in [(s,mats[0]) for s in range(3)]+[(3,mats[1]),(4,mats[2]),(7,mats[1]),(8,mats[2])]:
        q=matrix.tocoo();rows.extend((species*n+q.row).tolist());cols.extend(((species if species<7 else species-4)*n+q.col).tolist())
    return coo_matrix((np.ones(len(rows),bool),(rows,cols)),shape=(9*n,9*n)).tocsr()
def solve(job,times,deadline,checkpoint=None):
    case,ga,gb,k,mode,ref,parent=factory(job);n=ga['n'];mats=matrices(case,ga,gb);initial=np.r_[case['initial'].ravel(),np.zeros(2*n)];times=np.asarray(times,float)
    if times[0]!=0 or np.any(np.diff(times)<=0):raise ValueError('Invalid grid')
    def fun(t,y):
        if time.perf_counter()>deadline:raise TimeoutError('Population bridge cap')
        return rhs(y,case,mats,k,mode)
    # Manual BDF stepping retains valid output prefixes if the solver is interrupted.
    from scipy.integrate import BDF
    solver=BDF(fun,0.,initial,float(times[-1]),jac_sparsity=sparsity(n,mats),**SETTINGS[job['setting']]);saved=[initial.copy()];index=1
    if checkpoint:checkpoint(times[:1],np.array(saved))
    while solver.status=='running':
        solver.step()
        if solver.status=='failed':raise RuntimeError('Population BDF failed')
        while index<len(times) and times[index]<=solver.t:
            saved.append(solver.dense_output()(times[index]));index+=1
            if checkpoint:checkpoint(times[:index],np.array(saved))
    if index!=len(times):raise RuntimeError('Population bridge incomplete')
    return np.array(saved)
def measure(job,state,times):
    n=500;pop=state[:,:7*n].reshape(-1,7,n);return metrics(pop,times,full=False)
def check(job,state,times):
    case,ga,gb,k,mode,ref,parent=factory(job)
    n=ga['n'];pop=state[:,:7*n].reshape(-1,7,n);imports=state[:,7*n:].reshape(-1,2,n)
    initial=case['initial'];scale=np.maximum(1,initial.max(axis=1))
    with np.load(ref) as reference:
        old=reference['state'][:len(times),:7*n].reshape(-1,7,n)
        np.testing.assert_array_equal(reference['times'][:len(times)],times)
    error=float(np.max(abs(pop-old)/scale[None,:,None]))
    residual=pop.sum(axis=(1,2))-initial.sum()-imports.sum(axis=(1,2))
    balance=float(np.max(abs(residual)))/max(1,float(initial.sum()))
    minimum=float(state.min());assert np.all(np.isfinite(state))
    if mode=='FINITE':assert not np.any(imports)
    passed=error<=1e-5 and balance<=1e-10 and minimum>=-1e-9
    return dict(passed=bool(passed),reference_error=error,balance=balance,minimum=minimum)
