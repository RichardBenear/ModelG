"""MG129 bounded actual count-writer timing, no formal production."""
from prepare_mg12h3_bridge import ROOT,fixture,config,instantiate,matches,seed_for,atomic_json,np,json,time,hashlib,os
from mg12h3_exchangeable import simulate,tick
BASE=ROOT/'results/MG12H3/efficiency/timing'
FILES=['src/mg12h3_exchangeable.py','src/check_mg12h3_exchangeable.py','src/prepare_mg12h3_efficiency.py','spec/MG12H3_efficiency_timing.md','src/mg12h1_reference.py','src/mg12h2.py','src/mg12h1.py','src/prepare_mg12h3_bridge.py']
def hashes():return {p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in FILES}
def audit(state,counts,trace):
    output=state.counts();available=output.copy();last=-1;accepted=0
    for row in trace:
        if row['step']!=last:available=output.copy();last=row['step']
        change=np.asarray(row['increment']);reserve=np.asarray(row['reserved']);found=False
        for event in matches(instantiate(state,available)):
            if event.rule!=row['rule']:continue
            d=np.zeros_like(output);r=np.zeros_like(output);before={t:(i,s) for t,i,s in event.before};names=('G','X','Y','A','B','Z','W')
            for t,i,s in event.before:r[names.index(s),i]+=1
            for t,i,s in event.after:
                j,v=before[t];d[names.index(v),j]-=1;d[names.index(s),i]+=1
            if np.array_equal(d,change) and np.array_equal(r,reserve):found=True;break
        assert found and np.all(reserve<=available);available-=reserve;output+=change;assert output.min()>=0;accepted+=1
    np.testing.assert_array_equal(output,counts[-1]);return accepted

def plan():
    if (BASE/'plan.json').exists():raise ValueError('Plan exists')
    seeds=[[seed_for(12,c,0,12900,r) for r in range(32)] for c in range(4)];assert len(set(v for group in seeds for v in group))==128
    atomic_json(BASE/'plan.json',dict(source_hashes=hashes(),seeds=seeds,elapsed_cap=60,paths=128,formal=False))
def run():
    p=json.loads((BASE/'plan.json').read_text());assert p['source_hashes']==hashes() and not (BASE/'attempt.json').exists();start=time.perf_counter();deadline=start+60;atomic_json(BASE/'attempt.json',dict(status='STARTED'));timings=[];accepted=0
    for c in range(4):
        state=fixture(c);arrays=[];clock=time.perf_counter()
        with (BASE/f'c{c}_trace.jsonl').open('x') as f:
            for seed in p['seeds'][c]:
                result=simulate(state.counts(),config(state),[0.,.05,.1],seed,deadline=deadline);accepted+=audit(state,result['counts'],result['trace']);arrays.append(result['counts']);f.write(json.dumps(dict(seed=seed,trace=result['trace']))+'\n')
        np.savez_compressed(BASE/f'c{c}.npz',counts=np.asarray(arrays),seeds=p['seeds'][c]);timings.append(time.perf_counter()-clock)
    scales=[]
    for scale in (10,100):
        state=fixture(3);clock=time.perf_counter();result=tick(state.counts()*scale,config(state),.001,np.random.default_rng(seed_for(12,scale,0,12901,0)),deadline=deadline);scales.append(dict(scale=scale,elapsed=time.perf_counter()-clock,events=len(result['trace']),counts=result['counts'].tolist(),trace=result['trace']))
    assert time.perf_counter()<deadline and p['source_hashes']==hashes();atomic_json(BASE/'results.json',dict(passed=True,paths=128,events_audited=accepted,timing=timings,elapsed_seconds=time.perf_counter()-start,scales=scales,source_hashes=hashes(),forecast_131072_paths_seconds=2*sum(timings)*1024))
if __name__=='__main__':
    import sys
    plan() if '--plan' in sys.argv else run()
