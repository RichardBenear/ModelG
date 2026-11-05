"""Run with a ModelG checkout path. Reads saved records only."""
from collections import Counter
from math import comb, sqrt
import sys, json
from pathlib import Path
from scipy.stats import binomtest, fisher_exact, norm, chi2

root = Path(sys.argv[1])/'results/MG9b'
reported = json.loads((root/'event_shapes_posthoc.json').read_text())
by_b = {}
pooled_null = Counter({0: 1})
total_expected = total_stat = 0
for b in ('29.20','29.35'):
    cells = []
    for folder in sorted((root/'runs').glob(f'b{b}_B_*'), key=lambda p:int(p.name.split('_')[-1])):
        hist = Counter()
        files = list(folder.glob('*.json'))
        records = sorted((json.loads(f.read_text()) for f in files), key=lambda r:r['rep'])
        cum = 0
        for i,r in enumerate(records,1):
            assert r['rep']==i-1 and r['complete']
            hits = [birth for birth in r['births'] if birth['cls']=='persistent']
            first = min(hits,key=lambda h:h['t_start']) if hits else None
            scored = first is not None and 50<=first['t_start']<=480
            assert scored==r['event']
            if scored:
                c, = [c for c in r['confirmations'] if c['sc']==first['s_conf']]
                assert c['complete'] and c['cls']=='persistent' and c['n_deep_end']>=1
                hist[c['n_deep_end']]+=1
                cum+=1
            if i%16==0: assert (i==len(records))==(cum>=20 or i==64)
        omega=int(folder.name.split('_')[-1]); n=sum(hist.values()); k=hist[1]
        ci=binomtest(k,n).proportion_ci(method='wilson')
        orig, = [c for c in reported[b]['cells'] if c['Om']==omega]
        assert (n,k,n-k)==(orig['events'],orig['single'],orig['multiple'])
        assert hist=={int(k):v for k,v in orig['n_deep_hist'].items()}
        assert max(abs(100*x-y) for x,y in zip(ci,orig['wilson95']))<1e-6
        cells.append(dict(omega=omega,n=n,k=k,wilson=list(ci)))
    n0,n1,n2=[c['n'] for c in cells]
    N=n0+n1+n2; K=sum(c['k'] for c in cells); p=K/N
    score=cells[2]['k']-cells[0]['k']; expected=p*(n2-n0)
    z=(score-expected)/sqrt(p*(1-p)*(n0+n2-(n2-n0)**2/N))
    trend_p=2*norm.sf(abs(z))
    hom_p=chi2.sf(sum((c['k']-p*c['n'])**2/(c['n']*p*(1-p)) for c in cells),2)
    assert abs(trend_p-reported[b]['trend_p'])<1e-10 and abs(hom_p-reported[b]['homog_p'])<1e-10
    by_b[b]=dict(cells=cells,K=K,N=N,trend_p=trend_p,homogeneity_p=hom_p)
    # Additional post-hoc test of a common ordering within b. Fix each row's
    # event counts and each b's single total; scores -1,0,1 per 15 Omega.
    dist=Counter()
    for a in range(n0+1):
        for c in range(n2+1):
            mid=K-a-c
            if 0<=mid<=n1: dist[c-a]+=comb(n0,a)*comb(n1,mid)*comb(n2,c)
    assert sum(dist.values())==comb(N,K)
    joined=Counter()
    for i,v in pooled_null.items():
        for j,w in dist.items(): joined[i+j]+=v*w
    pooled_null=joined;total_expected+=expected;total_stat+=score
tail=sum(v for x,v in pooled_null.items() if abs(x-total_expected)>=abs(total_stat-total_expected)-1e-12)
fisher=fisher_exact([[r['K'],r['N']-r['K']] for r in by_b.values()]).pvalue
out=dict(cells_by_b=by_b,aggregate_fisher_p=fisher,exploratory_stratified_exact_score_p=tail/sum(pooled_null.values()),
         statistic=total_stat,conditional_mean=total_expected,upper_onesided_95_after_zero_of_four=1-.05**.25)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
