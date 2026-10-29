"""Independent scoring from primitive run times/outcomes; does not import MG3 evaluation code."""
import json,pickle,collections,math
from pathlib import Path
import numpy as np
from scipy import stats
ROOT=Path(__file__).parent; INPUT=ROOT/'inputs'
class Limited(pickle.Unpickler):
    def find_class(self,m,n):
        allow={('numpy._core.multiarray','_reconstruct'):np._core.multiarray._reconstruct,('numpy.core.multiarray','_reconstruct'):np._core.multiarray._reconstruct,('numpy','ndarray'):np.ndarray,('numpy','dtype'):np.dtype,('numpy._core.multiarray','scalar'):np._core.multiarray.scalar,('numpy.core.multiarray','scalar'):np._core.multiarray.scalar}
        if (m,n) not in allow:raise ValueError((m,n))
        return allow[m,n]
def load(p):
    with p.open('rb') as f:return Limited(f).load()
def rows(cell):return [load(p) for p in sorted((INPUT/'results/MG3/runs'/cell).glob('*.pkl'))]
def prim(r):
    t=None if r['first'] is None else r['first']['t_start']
    startup=t is not None and t<50
    event=t is not None and 50<=t<=480
    E=0 if startup else (t-50 if event else 430)
    E1=0 if startup else min(E,150);E2=E-E1
    return dict(startup=startup,event=event,exposure=E,expo_w1=E1,expo_w2=E2,in_w1=event and t<200)
def rate(rr):
    pp=[prim(r) for r in rr];n=sum(p['event'] for p in pp);E=sum(p['exposure'] for p in pp)
    E1=sum(p['expo_w1'] for p in pp);E2=sum(p['expo_w2'] for p in pp)
    n1=sum(p['event'] and p['in_w1'] for p in pp);n2=n-n1
    lo=stats.chi2.ppf(.025,2*n)/2 if n else 0;hi=stats.chi2.ppf(.975,2*n+2)/2
    p=stats.binomtest(n1,n,E1/E).pvalue if n else None
    return dict(runs=len(rr),events=n,startup=sum(p['startup'] for p in pp),exposure=E,k=n/(100*E) if E else None,ci=[lo/(100*E),hi/(100*E)] if E else None,n1=n1,n2=n2,E1=E1,E2=E2,p=p)
def ratio(a,b):
    n=a['events'];m=b['events'];E=a['exposure'];F=b['exposure']
    lo=stats.beta.ppf(.025,n,m+1);hi=stats.beta.ppf(.975,n+1,m)
    return dict(ratio=(n/E)/(m/F),ci=[lo/(1-lo)*F/E,hi/(1-hi)*F/E])
def wave(x,f):
    q=len(x)//2;y=(f[q:]+f[q::-1][:len(x)-q])/2;r=x[q:]
    j=np.flatnonzero(y[:-1]*y[1:]<0);z=r[j]-y[j]*(r[j+1]-r[j])/(y[j+1]-y[j]);z=z[(z>1)&(z<7)]
    return 2*np.mean(np.diff(z)) if len(z)>1 else None
if __name__=='__main__':
    cells={p.name:rows(p.name) for p in (INPUT/'results/MG3/runs').iterdir() if p.is_dir()}
    mismatch=[];consistency=[]
    for c,rr in cells.items():
        for r in rr:
            for k,v in prim(r).items():
                if r[k]!=v:mismatch.append([c,r['rep'],k,r[k],v])
            t=None if r['first'] is None else r['first']['t_start']
            expected=500 if t is None else min(500,10*math.ceil((t+50)/10))
            if expected!=r['t_end']:consistency.append([c,r['rep'],'stop',expected,r['t_end']])
            for b in r['births']:
                if b['t_conf']!=b['t_start']+20:consistency.append([c,r['rep'],'confirmation_time',b])
    R={c:rate(rr) for c,rr in cells.items()};stage={};c2={};ref=np.load(INPUT/'results/MG003/ref_profiles.npz')['h0.2_b29.0'];ld=wave(ref[0],ref[2])
    for om in [300,250,200,150,120,100]:
        rr=cells[f'A_{om}'];stage[om]=dict(count=sum(r['first'] is not None and 50<=r['first']['t_start']<=150 for r in rr),crowd=np.mean([np.mean(r['n_dets'][:99]) for r in rr]))
    for om in [150,120,100]:
        rr=cells[f'B_{om}'];cnt=collections.Counter(b['cls'] for r in rr for b in r['births'] if not b['shared'])
        eligible=[r for r in rr if r['profile'] is not None and r['profile']['eligible']]
        lam=wave(ref[0],np.mean([r['profile']['full'][1] for r in eligible],0)) if eligible else None
        refs=load(INPUT/f'results/MG3/refs/R_{om}.pkl')
        c2[om]=dict(outcomes=dict(cnt),fraction=cnt['persistent']/sum(cnt.values()),n_nuc=len(eligible),n_ref=sum(r['eligible'] for r in refs['runs']),wavelength=lam,deterministic_wavelength=ld,relative_wavelength_error=lam/ld-1 if lam else None)
    xx=np.array([150,120,100]);ww=np.array([R[f'B_{o}']['events'] for o in xx]);yy=np.log([R[f'B_{o}']['k'] for o in xx])
    A=np.column_stack([np.ones(3),xx]);cov=np.linalg.inv(A.T@(ww[:,None]*A));beta=cov@(A.T@(ww*yy));res=yy-A@beta
    se=np.sqrt(cov[1,1]);sse=np.sum(ww*res**2);r2=1-sse/np.sum(ww*(yy-np.average(yy,weights=ww))**2)
    fit=dict(slope=beta[1],ci=[beta[1]-1.96*se,beta[1]+1.96*se],ci_t=[beta[1]-stats.t.ppf(.975,1)*se*np.sqrt(sse),beta[1]+stats.t.ppf(.975,1)*se*np.sqrt(sse)],r2=r2)
    fish=stats.fisher_exact([[R['B_120']['startup'],16-R['B_120']['startup']],[R['V5_120']['startup'],16-R['V5_120']['startup']]])
    birth_runs=[dict(cell=c,rep=r['rep'],first=r['first']['t_start'] if r['first'] else None,t_end=r['t_end'],births=r['births']) for c,rr in cells.items() for r in rr if r['births']]
    wall=sum(r['wall'] for rr in cells.values() for r in rr)
    for p in list((INPUT/'results/MG3/refs').glob('*.pkl'))+list((INPUT/'results/MG3/pilots').glob('*.pkl')):wall+=load(p)['wall']
    out=dict(commit=json.loads((ROOT/'download_manifest.json').read_text())['commit'],bookkeeping_mismatches=mismatch,stop_or_confirmation_mismatches=consistency,rates=R,stageA=stage,C1={o:R[f'A_{o}']['events']+R[f'B_{o}']['events'] for o in [150,120,100]},C2=c2,C3=fit,V5=ratio(R['V5_120'],R['B_120']),C4_ratio=ratio(R['C4_b29.3_120'],R['B_120']),startup_fisher_p=fish.pvalue,core_hours=wall/3600,birth_runs=birth_runs)
    (ROOT/'independent_rescore.json').write_text(json.dumps(out,indent=2,default=lambda a:a.item() if hasattr(a,'item') else str(a)))
    print(json.dumps({k:v for k,v in out.items() if k not in ['rates','birth_runs']},indent=2,default=str))
    print('Cell rates:',json.dumps(R,indent=2))
