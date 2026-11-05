"""Independent review of MG035. Usage: python audit.py CHECKOUT. No ModelG imports."""
import json
import sys
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp, softmax
from scipy.stats import chi2

root=Path(sys.argv[1])/'results'
p=root/'MG9b'
def read(path): return json.loads(path.read_text())
ev=read(p/'evaluation.json')
rows=[read(f) for f in (p/'runs').glob('*/*.json')]
assert len({tuple(r['key']) for r in rows})==len(rows)==280
assert all(r['b']==r['b_det'] and abs(r['tau']-1/1500)<1e-17 for r in rows)
def records(cell): return sorted((r for r in rows if r['cell']==cell),key=lambda r:r['rep'])
def passage(r):
    times=[v['t_start'] for v in r['births'] if v['cls']=='persistent']
    t=min(times,default=None)
    assert t==r['t_first']
    startup=t is not None and t<50
    event=t is not None and 50<=t<=480
    return startup,event,(0 if startup else t-50 if event else 430)
for r in rows:
    assert passage(r)==(r['startup'],r['event'],r['exposure'])
    if r['item'] in ('A','B'): assert r['complete']
    else: assert r['t_end']==(20 if r['item']=='pilot' else 50)
    assert all(c['complete'] for c in r['confirmations'])

history=read(root/'MG3b/evaluation.json')['part2']['cells']
groups={29.:np.array([(float(k),v['n'],v['exposure']) for k,v in history.items()])}
out={'revision':'5e476adaeb3b57b231c88ebb3872c762df10becc','arms':{},'record_completion':dict(Counter(str((r['item'],r['complete'])) for r in rows))}
for tag in ['29.20','29.35']:
    b=float(tag); prefix='b'+tag
    pa,pb=(read(p/f'{part}_{tag}.json') for part in ('partA','partB'))
    ranked=[]; visit=[]
    for om in [250,200,170,150,140,130,120,110,100,90,80]:
        A=records(f'{prefix}_A_{om}'); assert len(A)==4
        clip,=records(f'{prefix}_CLIP_{om}'); assert clip['clip_frac']<=1e-4
        crowd=np.mean([r['mean_dets_t_lt_50'] for r in A]); assert crowd<=2
        count=sum(passage(r)[1] for r in A); visit.append(om)
        if count>=2: ranked.append(om)
        if count>=3: break
    top=max(ranked); assert visit==[c['Om'] for c in pa['cells']] and top==pa['Om_T']==pb['Om_T']
    accepted=[top]
    for om in [top-15,top-30,top+15,top+30]:
        scr=records(f'{prefix}_SCR_{om}'); assert len(scr)==5
        if scr[0]['clip_frac']<=1e-4 and np.mean([r['mean_dets_t_lt_50'] for r in scr[1:]])<=2: accepted.append(om)
        if len(accepted)==3: break
    assert accepted==pb['accepted']
    cells=[]; details=[]; shapes=[]; confirmations=[]
    for om,meta in zip(accepted,pb['cells'],strict=True):
        rr=records(f'{prefix}_B_{om}'); hits=np.array([passage(r)[1] for r in rr])
        assert [r['rep'] for r in rr]==list(range(len(rr)))
        first_batch=next((k for k in range(16,len(rr)+1,16) if hits[:k].sum()>=20 or k==64),None)
        assert first_batch==len(rr)==meta['runs'] and hits.sum()==meta['events']
        n=int(hits.sum()); E=sum(passage(r)[2] for r in rr); starts=sum(passage(r)[0] for r in rr)
        k=n/(100*E); cells.append([om,n,E]); stat=ev['per_b'][tag]['cells'][str(om)]
        assert n==stat['events'] and starts==stat['startup'] and E==stat['exposure'] and k==stat['k']
        details.append(dict(omega=om,runs=len(rr),events=n,exposure=E,startup=starts,rate=k))
        for r in rr:
            if passage(r)[1]:
                match,=[c for c in r['confirmations'] if c['sc']==r['first']['s_conf']]
                assert match['cls']=='persistent'
                shapes.append(match)
    groups[b]=np.array(cells)
    singles=sum(c['n_deep_end']==1 for c in shapes)
    assert singles==ev['per_b'][tag]['event_shape']['single_minimum']
    for om in visit:
        confirmations.extend((r['cell'],r['rep'],c['sc'],c['cls']) for r in records(f'{prefix}_A_{om}') for c in r['confirmations'])
    for om in accepted:
        confirmations.extend((r['cell'],r['rep'],c['sc'],c['cls']) for r in records(f'{prefix}_B_{om}') for c in r['confirmations'])
    v2=read(p/f'vb2_{tag}.json')['pairs']
    assert [(x['cell'],x['rep'],x['sc'],x['std_cls']) for x in v2]==confirmations[:10]
    assert all(x['std_complete'] and x['late_complete'] and x['late_cls']==x['std_cls'] for x in v2)
    v3=read(p/f'vb3_{tag}.json')['pairs']
    assert [(x['rep'],x['sc'],x['std_cls']) for x in v3]==[(rep,sc,cls) for cell,rep,sc,cls in confirmations if cell==f'{prefix}_B_{top}']
    assert v3 and all(x['std_complete'] and x['tight_complete'] and x['std_cls']==x['tight_cls'] for x in v3)
    ref=read(p/f'vb1_{tag}.json')
    final=np.load(p/f'vb1_{tag}_final.npz')['phi'][1]
    original=np.load(root/'MG003/ref_profiles.npz')[f'h0.2_b{b}'][2]
    assert ref['complete'] and ref['one_chain'] and ref['n_minima']==[1]*11 and abs(ref['drift'])<=ref['lam']/4
    assert abs(final.min()/original.min()-1)<.05 and sum((final<np.roll(final,1))&(final<=np.roll(final,-1))&(final<-abs(ref['A'])/2))==1
    out['arms'][tag]=dict(cells=details,single_minimum=singles,shape_hist=dict(Counter(c['n_persistent'] for c in shapes)),vb2_pairs=len(v2),vb3_pairs=len(v3))

# Profile out each b intercept. Solve the moment equations directly, using the
# observed event-weighted Omega mean and its predicted exposure-weighted mean.
def profile(data,s):
    om,n,E=data.T; logw=np.log(E)-s*om
    return float(n@(-s*om)-sum(n)*logsumexp(logw)),float(sum(n)*(softmax(logw)@om)-n@om)
def bounds(fn,x,peak):
    f=lambda s:2*(peak-fn(s))-chi2.ppf(.95,1)
    return [brentq(f,-1,x),brentq(f,x,1)]
separate=0.
for b,data in groups.items():
    s=brentq(lambda s:profile(data,s)[1],-1,1); ll=profile(data,s)[0]; separate+=ll
    ci=bounds(lambda s:profile(data,s)[0],s,ll); tag=f'{b:.2f}'
    np.testing.assert_allclose([s,*ci],[ev['per_b'][tag]['fit']['S'],*ev['per_b'][tag]['fit']['ci']],atol=1e-7)
    out['arms'].setdefault(tag,{}).update(S=s,CI=ci)
def joint(t):
    def scores(a):return [profile(data,a+t*(b-29.175)) for b,data in groups.items()]
    a=brentq(lambda a:sum(g for _,g in scores(a)),-2,2)
    pts=scores(a)
    return sum(l for l,_ in pts),sum((b-29.175)*g for b,(_,g) in zip(groups,pts)),a
t=brentq(lambda t:joint(t)[1],-1,1); ll,_,a=joint(t)
ci=bounds(lambda t:joint(t)[0],t,ll); lr=2*(separate-ll)
np.testing.assert_allclose([a,t,*ci,lr],[ev['trend']['S0'],ev['trend']['S1'],*ev['trend']['ci'],ev['trend']['lr']['stat']],atol=1e-7)
cal=read(p/'calibration/calibration.json')
supported=[q for q in cal['cases'] if q['S1_true'] is not None and q['branches']['B3']['n']>=400]
assert {q['S1_true']==0 for q in supported}=={False,True}
for q in supported:
    v=q['branches']['B3']; assert v['coverage']['k']/v['n']>=.93
    if q['S1_true']==0: assert all(v[sign]['k']/v['n']<=.035 for sign in ('neg','pos'))
lg=read(p/'ledger.json'); assert lg['budget_stop'] is None
entries=lg['entries']; rt={(a['cell'],a['rep']):a for a in entries if a['item']!='vb'}
assert len(rt)==280 and all(rt[(r['cell'],r['rep'])]['wall_s']==r['wall'] for r in rows)
hours=sum(e['wall_s'] for e in entries)/3600; assert abs(hours-ev['ledger_core_h'])<1e-9 and hours<50
assert ci[0]<0<ci[1] and ev['primary']['confirmatory']=='NOT RESOLVED' and ev['qualifiers']==[]
out.update(trend=dict(S0=a,S1=t,CI=ci,lr=lr,lr_p=chi2.sf(lr,1)),hours=hours,largest_jobs=[(e['cell'],e['rep'],e['wall_s']) for e in sorted(entries,key=lambda e:-e['wall_s'])[:2]],branch='B3',outcome='NOT RESOLVED')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
