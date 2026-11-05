"""MG9b independent recomputation from raw run JSONs (does not import mg9b or evaluate_mg9b). Usage: python src/check_mg9b.py"""
import os; os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import json,glob,os,numpy as np
from scipy.optimize import minimize, brentq
from scipy import stats
R='results/MG9b'
cells={}
pb={b:json.load(open(f'{R}/partB_{b}.json')) for b in ('29.20','29.35')}
for b in ('29.20','29.35'):
    used={c['Om']:c['runs'] for c in pb[b]['cells']}
    print(b,'partB cells',[(c['Om'],c['runs'],c.get('capped'),c.get('budget_short'),c.get('unstarted')) for c in pb[b]['cells']])
    cl=[]
    for Om,nu in used.items():
        fs=glob.glob(f'{R}/runs/b{b}_B_{Om:g}/*.json'); rs=sorted((json.load(open(f)) for f in fs),key=lambda r:r['rep'])
        assert len(rs)==nu,(b,Om,len(rs),nu)
        assert all(r['b_det']==r['b']==float(b) for r in rs)
        # independent: event = first persistent birth; exposure = time from 50 to min(t_first,480) for non-startup
        n=sum(1 for r in rs if r['complete'] and r['event']); E=sum(r['exposure'] for r in rs if r['complete'])
        inc=sum(1 for r in rs if not r['complete'])
        # recompute exposure from t_first
        E2=0
        for r in rs:
            if r['startup']: continue
            E2+= (min(r['t_first'],480) if r['event'] else 480)-50
        cl.append((Om,n,E,E2,inc))
    cells[float(b)]=cl
    print(b,cl)
ev=json.load(open('results/MG3b/evaluation.json'))
cells[29.0]=[(float(o),int(c['n']),float(c['exposure']),None,0) for o,c in ev['part2']['cells'].items()]
L=100.
def nll(params,bs,Sfun):
    s=0
    for i,b in enumerate(bs):
        a=params[i]; S=Sfun(b,params)
        for Om,n,E,*_ in cells[b]:
            lm=np.log(L*E)+a-S*Om; s-= n*lm-np.exp(lm)
    return s
bs=[29.0,29.2,29.35]
def fit(Sfun,np_extra,x0e):
    x0=[np.log(1e-4)+0.05*120]*3+x0e
    r=minimize(nll,x0,args=(bs,Sfun),method='BFGS',options=dict(gtol=1e-10));return r
for i,b in enumerate(bs):
    r=minimize(lambda p:nll([p[0]],[b],lambda bb,pp:p[1]),[np.log(1e-4)+6,0.05],method='Nelder-Mead',options=dict(xatol=1e-12,fatol=1e-12,maxiter=20000))
    print('S_b',b,r.x[1])
tr=fit(lambda b,p:p[3]+p[4]*(b-29.175),2,[0.05,0.0])
print('trend S0,S1',tr.x[3],tr.x[4],'ll',-tr.fun)
sep=fit(lambda b,p:p[3+bs.index(b)],3,[0.05,0.05,0.05])
lr=2*(tr.fun-sep.fun);print('LR',lr,stats.chi2.sf(lr,1))
def prof(s1):
    r=minimize(lambda p:nll(p,bs,lambda b,pp:pp[3]+s1*(b-29.175)),tr.x[:4],method='BFGS',options=dict(gtol=1e-10));return r.fun
f=lambda s1:prof(s1)-tr.fun-stats.chi2.ppf(.95,1)/2
print('CI',brentq(f,-1,tr.x[4]),brentq(f,tr.x[4],1))
for b in ('29.20','29.35'):
    for c in pb[b]['cells']:
        rs=sorted((json.load(open(f)) for f in glob.glob(f"{R}/runs/b{b}_B_{c['Om']:g}/*.json")),key=lambda r:r['rep'])
        cum=[sum(r['event'] for r in rs[:k]) for k in range(16,len(rs)+1,16)]
        ok=all(x<20 for x in cum[:-1]) and (cum[-1]>=20 or len(rs)==64)
        print('stop20',b,c['Om'],cum,'OK' if ok else 'VIOLATION')
