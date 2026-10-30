"""Independent MG3b scoring from first-event times, using no repository evaluation functions."""
from pathlib import Path
import json,pickle,csv,collections,math
import numpy as np
from scipy import stats,optimize,special
ROOT=Path(__file__).parent; INPUT=ROOT/'inputs'; DATA=INPUT/'results/MG3b'
class Limited(pickle.Unpickler):
    def find_class(self,m,n):
        d={('numpy._core.multiarray','_reconstruct'):np._core.multiarray._reconstruct,('numpy.core.multiarray','_reconstruct'):np._core.multiarray._reconstruct,('numpy','ndarray'):np.ndarray,('numpy','dtype'):np.dtype,('numpy._core.multiarray','scalar'):np._core.multiarray.scalar,('numpy.core.multiarray','scalar'):np._core.multiarray.scalar}
        if (m,n) not in d:raise ValueError((m,n))
        return d[m,n]
def load(p):
    with p.open('rb') as f:return Limited(f).load()
def rows(cell):return sorted([load(p) for p in (DATA/'runs'/cell).glob('*.pkl')],key=lambda r:r['rep'])
def first(r):return None if r['first'] is None else float(r['first']['t_start'])
def primitives(r,burn=50):
    t=first(r);excluded=t is not None and t<burn;ev=t is not None and burn<=t<=480
    e=0 if excluded else (t-burn if ev else 480-burn)
    return excluded,ev,e
def ci(n,E):return [0 if n==0 else stats.chi2.ppf(.025,2*n)/(200*E),stats.chi2.ppf(.975,2*n+2)/(200*E)]
def summary(rr,burn=50):
    p=[primitives(r,burn) for r in rr];n=sum(x[1] for x in p);E=sum(x[2] for x in p)
    return dict(runs=len(rr),risk_set=sum(not x[0] for x in p),events=n,excluded=sum(x[0] for x in p),exposure=E,k=n/(100*E),ci=ci(n,E))
def ratio(n1,E1,n2,E2,level=.95):
    a=1-level;lo=stats.beta.ppf(a/2,n2,n1+1) if n2 else 0;hi=stats.beta.ppf(1-a/2,n2+1,n1) if n1 else 1
    return dict(ratio=(n2/E2)/(n1/E1),ci=[lo/(1-lo)*E1/E2,hi/(1-hi)*E1/E2 if hi<1 else float('inf')])
def stationarity(rr):
    p=[primitives(r) for r in rr];E=np.array([x[2] for x in p]);total=E.sum()
    s=optimize.brentq(lambda u:np.minimum(E,u).sum()-total/2,0,430)
    ne=sum(x[1] and first(r)-50<s for x,r in zip(p,rr));nt=sum(x[1] for x in p)
    return dict(split=s+50,early_events=int(ne),late_events=int(nt-ne),early_exposure=float(np.minimum(E,s).sum()),late_exposure=float(np.maximum(E-s,0).sum()),p=stats.binomtest(ne,nt,.5).pvalue)
def glm(omegas,ns,es):
    x=np.array(omegas,float);x-=x.mean();n=np.array(ns,float);E=np.array(es,float);N=n.sum();off=np.log(100*E)
    def loglik(s):return s*(n*x).sum()-N*special.logsumexp(off+s*x)
    f=optimize.minimize_scalar(lambda s:-loglik(s),bounds=(-.3,.1),method='bounded',options={'xatol':1e-14});s=f.x
    # Root of the score supplies the final independent optimum.
    def score(s):
        mu=N*special.softmax(off+s*x)
        return ((n-mu)*x).sum()
    s=optimize.brentq(score,-.3,.1,xtol=1e-14)
    cutoff=stats.chi2.ppf(.95,1)/2
    lo=optimize.brentq(lambda v:loglik(s)-loglik(v)-cutoff,-.5,s)
    hi=optimize.brentq(lambda v:loglik(s)-loglik(v)-cutoff,s,.3)
    mu=N*special.softmax(off+s*x);info=(mu*x*x).sum()-(mu*x).sum()**2/N
    dev=2*np.sum(special.xlogy(n,n/mu)-n+mu)
    return dict(S=-s,profile_ci=[-hi,-lo],wald_ci=[-s-1.95996398454/np.sqrt(info),-s+1.95996398454/np.sqrt(info)],deviance=dev,df=len(n)-2,p=stats.chi2.sf(dev,len(n)-2))
if __name__=='__main__':
    cells={p.name:rows(p.name) for p in (DATA/'runs').iterdir() if p.is_dir()};bad=[];allkeys=[];fields={};stops={}
    for cell,rr in cells.items():
        hist=collections.Counter()
        for r in rr:
            allkeys.append(tuple(r['key']));t=first(r);ex,ev,E=primitives(r)
            if not cell.startswith('TIMING'):
                for k,v in [('startup',ex),('event',ev),('exposure',E),('complete',True)]:
                    if r[k]!=v:bad.append([cell,r['rep'],k,r[k],v])
                end=500 if t is None else min(500,10*math.ceil((t+50)/10))
                if r['t_end']!=end:bad.append([cell,r['rep'],'end',r['t_end'],end])
            for b in r['births']:
                c=r['conf_fields'][b['s_conf']]
                if b['t_conf']!=b['t_start']+20 or b['cls']!=c['cls']:bad.append([cell,r['rep'],'confirmation'])
            if r['first'] is not None:
                c=r['conf_fields'][r['first']['s_conf']]
                if not np.array_equal(c['field'],r['first']['field']):bad.append([cell,r['rep'],'first_field'])
            hist.update(c['cls'] for c in r['conf_fields'].values())
        fields[cell]=dict(hist)
        if cell.startswith('P2_'):
            bcounts=[sum(primitives(r)[1] for r in rr[:i]) for i in range(16,len(rr)+1,16)]
            stops[cell]=dict(cumulative_events=bcounts,valid=all(v<20 for v in bcounts[:-1]) and (bcounts[-1]>=20 or len(rr)==64))
    a=summary(cells['P1_tau']);b=summary(cells['P1_tau2']);r90=ratio(a['events'],a['exposure'],b['events'],b['exposure'],.9);r95=ratio(a['events'],a['exposure'],b['events'],b['exposure'])
    p2={};om=[130,120,110,100]
    for o in om:
        rr=cells['P1_tau' if o==120 else f'P2_{o}'];v=summary(rr);v['landmarks']={str(t):summary(rr,t) for t in [100,150]};v['stationarity']=stationarity(rr)
        # Early/late contributions recomputed from first-event times.
        p=[primitives(r) for r in rr];n1=sum(z[1] and first(r)<100 for z,r in zip(p,rr));E1=sum(min(z[2],50) for z in p);n2=v['events']-n1;E2=v['exposure']-E1
        v['early_late']=dict(early_n=int(n1),late_n=int(n2),early_E=E1,late_E=E2,**ratio(n1,E1,n2,E2));p2[o]=v
    ledger=list(csv.DictReader((DATA/'budget_ledger.csv').open()));tot=collections.defaultdict(float);dup=[];seen=set()
    for row in ledger:
        tot[row['item']]+=float(row['wall_s'])/3600
        if row['kind']=='run':
            if row['key'] in seen:dup.append(row['key'])
            seen.add(row['key'])
    # Reconstruct the original P2 12-hour submission check from completed-run costs.
    oldbudget=[];used=0;walls=[]
    for cell in ['P2_130','P2_110','P2_100']:
        rr=cells[cell]
        for i in range(0,len(rr),16):
            per=np.mean(walls) if walls else 7.2*3600/117
            projected=16*per*1.25
            oldbudget.append(dict(cell=cell,start_rep=i,used_h=used/3600,projected_h=projected/3600,allowed=bool(used+projected<=12*3600)))
            batch=[r['wall'] for r in rr[i:i+16]];walls.extend(batch);used+=sum(batch)
    counts={o:p2[o]['events'] for o in om};expo={o:p2[o]['exposure'] for o in om}
    out=dict(commit=json.loads((ROOT/'download_manifest.json').read_text())['commit'],records=sum(map(len,cells.values())),unique_keys=len(set(allkeys)),bookkeeping_mismatches=bad,part1=dict(tau=a,tau2=b,ratio90=r90,ratio95=r95,equivalent=.5<=r90['ci'][0] and r90['ci'][1]<=2,startup_fisher=stats.fisher_exact([[a['excluded'],a['risk_set']],[b['excluded'],b['risk_set']]]).pvalue),part2=p2,slope=glm(om,[counts[o] for o in om],[expo[o] for o in om]),slope_without100=glm(om[:-1],[counts[o] for o in om[:-1]],[expo[o] for o in om[:-1]]),confirmation_fields=fields,stop20=stops,ledger_by_item=dict(tot),ledger_total=sum(tot.values()),duplicate_run_ledger_keys=dup,original_budget_replay=oldbudget)
    (ROOT/'independent_rescore.json').write_text(json.dumps(out,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
    print(json.dumps(out,indent=2,default=str))
