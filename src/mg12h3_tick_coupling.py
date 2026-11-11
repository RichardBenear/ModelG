"""MG156 split-propensity strict-tick coupling, frozen marginals intended."""
from fractions import Fraction
import math,time,json
import numpy as np
from mg12h3_sparse import Tree

def pair(model,left,right,deltas,times,seed,deadline=math.inf,rng=None,sentinel=False,check_every=1000,max_events=1000000):
    if type(seed) is not int or not 0<=seed<2**32:raise ValueError('Invalid seed')
    ds=tuple(Fraction(str(d)) for d in deltas);obs=tuple(Fraction(str(t)) for t in times)
    if len(ds)!=2 or any(d<=0 for d in ds) or len(obs)<2 or obs[0]!=0 or any(b<=a for a,b in zip(obs,obs[1:])) or any((t/d).denominator!=1 for t in obs for d in ds):raise ValueError('Require positive ticks and common-boundary observations')
    rng=np.random.default_rng(seed) if rng is None else rng;work=[model.workspace(left),model.workspace(right)];joint=Tree(np.maximum(work[0].tree.values,work[1].tree.values));next_boundary=list(ds);clock=0.;frame=1;counts=[[work[s].output.reshape(7,model.n).copy()] for s in range(2)];trace=[];boundaries=[];hist=np.zeros((3,len(model.channels)),np.int64);events=0;last_available=None
    def check():
        for s in range(2):
            fresh=np.array([c.rate(work[s].available) for c in model.channels]);np.testing.assert_allclose(work[s].tree.values,fresh,rtol=1e-10,atol=1e-6);assert work[s].available.min()>=0 and work[s].output.min()>=0
        fresh=np.maximum(work[0].tree.values,work[1].tree.values);np.testing.assert_allclose(joint.values,fresh,rtol=1e-10,atol=1e-6);np.testing.assert_allclose(joint.total(),math.fsum(fresh),rtol=1e-10,atol=1e-6)
    totals=[sum(map(int,np.asarray(x).ravel())) for x in (left,right)]
    check()
    while frame<len(obs):
        if time.perf_counter()>deadline:raise TimeoutError('Coupling deadline')
        boundary=min(next_boundary);limit=float(boundary);total=joint.total();active=bool(work[0].positive or work[1].positive)
        if active:
            if not math.isfinite(total) or total<=0:raise ValueError('Invalid joint total')
            wait=float(rng.exponential(1/total))
            if not math.isfinite(wait) or wait<=0:raise ValueError('Invalid mark')
            future=clock+wait
        else:future=math.inf
        if future>=limit:
            clock=limit;which=[s for s in range(2) if next_boundary[s]==boundary];mask=sum(1<<s for s in which);last_available=[w.available.reshape(7,model.n).copy() for w in work]
            if sentinel:boundaries.append(dict(time=clock,mask=mask,available=[a.tolist() for a in last_available],output=[w.output.reshape(7,model.n).tolist() for w in work]))
            for s in which:work[s]=model.workspace(work[s].output.reshape(7,model.n));next_boundary[s]+=ds[s]
            joint=Tree(np.maximum(work[0].tree.values,work[1].tree.values));check()
            if obs[frame]==boundary:
                for s in range(2):counts[s].append(work[s].output.reshape(7,model.n).copy())
                frame+=1
            continue
        if events>=max_events:raise ValueError('Joint event cap')
        index=joint.select(float(rng.random()));a,b=work[0].tree.values[index],work[1].tree.values[index]
        if a==b:mask=3
        elif min(a,b)==0:mask=1 if a>b else 2
        else:mask=3 if float(rng.random())*max(a,b)<min(a,b) else (1 if a>b else 2)
        channel=model.channels[index];touched=set()
        for i,n in channel.reserve:touched.update(model.dependencies[i])
        for s in range(2):
            if mask&(1<<s):work[s].accept(index)
        for i in touched:joint.update(i,max(work[0].tree.values[i],work[1].tree.values[i]))
        if not (work[0].positive or work[1].positive):joint=Tree(np.zeros(len(model.channels)))
        clock=future;hist[{3:0,1:1,2:2}[mask],index]+=1;events+=1
        if sentinel:trace.append(dict(time=clock,mask=mask,channel=index))
        if check_every and events%check_every==0:check()
    check()
    for s in range(2):assert sum(map(int,work[s].output))==totals[s]
    return dict(counts=np.asarray(counts),histogram=hist,trace=trace,boundaries=boundaries,events=events,final_available_before_reset=np.asarray(last_available),rng_state=json.dumps(rng.bit_generator.state) if hasattr(rng,'bit_generator') else None)

def replay(model,left,right,deltas,times,result):
    # Independent chronological reservation replay of retained joint events.
    available=[np.asarray(left).copy().ravel(),np.asarray(right).copy().ravel()];output=[x.copy() for x in available];hist=np.zeros_like(result['histogram']);rows=[(x['time'],0,x) for x in result['trace']]+[(x['time'],1,x) for x in result['boundaries']];rows.sort(key=lambda x:(x[0],x[1]));expected=[Fraction(str(d)) for d in deltas];ds=expected.copy();saved=[[np.asarray(left).copy()],[np.asarray(right).copy()]];frame=1
    for t,kind,row in rows:
        if kind==0:
            channel=model.channels[row['channel']];assert t<float(min(expected));hist[{3:0,1:1,2:2}[row['mask']],row['channel']]+=1
            for s in range(2):
                if row['mask']&(1<<s):
                    assert channel.rate(available[s])>0
                    for i,n in channel.reserve:assert available[s][i]>=n;available[s][i]-=n
                    for i,v in channel.change:output[s][i]+=v
        else:
            boundary=min(expected);assert t==float(boundary);which=[s for s in range(2) if expected[s]==boundary];assert row['mask']==sum(1<<s for s in which)
            for s in range(2):np.testing.assert_array_equal(available[s].reshape(7,model.n),row['available'][s]);np.testing.assert_array_equal(output[s].reshape(7,model.n),row['output'][s])
            for s in which:available[s]=output[s].copy();expected[s]+=ds[s]
            if frame<len(times) and boundary==Fraction(str(times[frame])):
                for s in range(2):saved[s].append(output[s].reshape(7,model.n).copy())
                frame+=1
    np.testing.assert_array_equal(hist,result['histogram']);np.testing.assert_array_equal(saved,result['counts']);return True
