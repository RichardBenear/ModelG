"""MG146 factored compensator for linear mean contrasts; sampler unchanged."""
import math,time,json
import numpy as np

class Contrast:
    def __init__(self,model,weights):
        weights=np.asarray(weights,float)
        if weights.ndim==1:weights=weights[None,:]
        if weights.ndim!=2 or weights.shape[1]!=7*model.n or not np.isfinite(weights).all():raise ValueError('Invalid linear weights')
        self.model=model;self.weights=weights.copy();increments=[];groups={};parts={}
        for c in model.channels:
            change=np.array([math.fsum(float(weights[k,i])*v for i,v in c.change) for k in range(len(weights))]);increments.append(change);key=tuple(sorted(c.factors));parts.setdefault(key,[]).append(c.coefficient*change)
        self.increments=np.asarray(increments).reshape(len(model.channels),len(weights));keys=[];coeff=[]
        for key,values in parts.items():
            row=np.array([math.fsum(float(v[k]) for v in values) for k in range(len(weights))])
            if np.any(row):keys.append(key);coeff.append(row)
        self.factors=tuple(keys);self.coefficients=np.asarray(coeff).reshape(len(keys),len(weights));deps=[[] for _ in range(7*model.n)]
        for g,key in enumerate(keys):
            for i,_ in key:deps[i].append(g)
        self.dependencies=tuple(tuple(x) for x in deps)
    def factor(self,g,available):
        value=1
        for i,order in self.factors[g]:
            n=int(available[i])
            for k in range(order):value*=max(0,n-k)
        return float(value)
    def direct(self,flat):return self.weights@np.asarray(flat).ravel()
    def metadata(self):return dict(factors=self.factors,coefficients=self.coefficients.tolist(),weights=self.weights.tolist())

def tick(contrast,pop,delta,rng,deadline=math.inf,max_events=1000000,sentinel=False,check_every=10000):
    if not math.isfinite(delta) or delta<=0:raise ValueError('Invalid tick')
    model=contrast.model;work=model.workspace(pop);groups=len(contrast.factors);values=np.array([contrast.factor(g,work.available) for g in range(groups)]);exposures=np.zeros(groups);last=np.zeros(groups);hist=np.zeros(len(model.channels),np.uint32);clock=0.;events=0;trace=[];checks=[];before=json.dumps(rng.bit_generator.state);next_mark=None
    def check():
        fresh=np.array([c.rate(work.available) for c in model.channels]);np.testing.assert_allclose(work.tree.values,fresh,rtol=1e-10,atol=1e-6);np.testing.assert_allclose(work.tree.total(),math.fsum(fresh),rtol=1e-10,atol=1e-6);drift=fresh@contrast.increments;grouped=values@contrast.coefficients;np.testing.assert_allclose(grouped,drift,rtol=1e-10,atol=1e-10);checks.append(dict(events=events,max_drift_difference=float(np.max(abs(grouped-drift),initial=0))))
    check()
    while work.positive:
        if time.perf_counter()>deadline:raise TimeoutError('Contrast tick deadline')
        total=work.tree.total()
        if not math.isfinite(total) or total<=0:raise ValueError('Invalid total')
        wait=float(rng.exponential(1/total))
        if not math.isfinite(wait) or wait<=0:raise ValueError('Invalid mark')
        future=clock+wait
        if future>=delta:next_mark=future;break
        if events>=max_events:raise ValueError('Event cap')
        index=work.tree.select(float(rng.random()));channel=model.channels[index];touched=set()
        for i,n in channel.reserve:touched.update(contrast.dependencies[i])
        for g in touched:exposures[g]+=values[g]*(future-last[g]);last[g]=future
        work.accept(index)
        for g in touched:values[g]=contrast.factor(g,work.available)
        clock=future;hist[index]+=1;events+=1
        if sentinel:trace.append(dict(time=clock,channel=index))
        if check_every and events%check_every==0:check()
    exposures+=values*(delta-last);check();integral=exposures@contrast.coefficients
    if not np.isfinite(integral).all() or np.any(exposures<0):raise ValueError('Invalid exposures')
    return dict(counts=work.output.reshape(7,model.n),available=work.available.reshape(7,model.n),histogram=hist,exposures=exposures,integral=integral,direct=contrast.direct(work.output),events=events,trace=trace,checks=checks,rng_before=before,rng_after=json.dumps(rng.bit_generator.state),last_accepted_time=clock,next_mark=next_mark)

def aggregate_audit(contrast,initial,result):
    model=contrast.model;hist=np.asarray(result['histogram']);exposure=np.asarray(result['exposures']);assert hist.shape==(len(model.channels),) and hist.dtype.kind in 'iu' and exposure.shape==(len(contrast.factors),) and np.isfinite(exposure).all() and np.all(exposure>=0);reserve=np.zeros(7*model.n,np.int64);change=reserve.copy()
    for c,n in zip(model.channels,hist):
        if n:
            for i,v in c.reserve:reserve[i]+=int(n)*v
            for i,v in c.change:change[i]+=int(n)*v
    np.testing.assert_array_equal(np.asarray(initial).ravel()-reserve,np.asarray(result['available']).ravel());np.testing.assert_array_equal(np.asarray(initial).ravel()+change,np.asarray(result['counts']).ravel());assert np.min(result['available'])>=0 and np.min(result['counts'])>=0 and np.sum(result['counts'])==np.sum(initial)
    np.testing.assert_allclose(exposure@contrast.coefficients,result['integral'],rtol=1e-12,atol=1e-12)

def trace_audit(contrast,initial,result,delta):
    # Independent full-channel rate integrals at each retained event interval.
    available=np.asarray(initial).copy().ravel();previous=0.;integral=np.zeros(len(contrast.weights));hist=np.zeros(len(contrast.model.channels),np.uint32)
    for row in result['trace']:
        t=float(row['time']);assert previous<t<delta;rates=np.array([c.rate(available) for c in contrast.model.channels]);integral+=(t-previous)*(rates@contrast.increments);c=contrast.model.channels[row['channel']];assert c.rate(available)>0
        for i,v in c.reserve:assert available[i]>=v;available[i]-=v
        hist[row['channel']]+=1;previous=t
    integral+=(delta-previous)*(np.array([c.rate(available) for c in contrast.model.channels])@contrast.increments);np.testing.assert_array_equal(hist,result['histogram']);np.testing.assert_array_equal(available,np.asarray(result['available']).ravel());np.testing.assert_allclose(integral,result['integral'],rtol=1e-10,atol=1e-11)
    return integral

def save_tick(path,contrast,initial,result):
    np.savez_compressed(path,initial=initial,final=result['counts'],available=result['available'],histogram=result['histogram'],exposures=result['exposures'],integral=result['integral'],direct=result['direct'],metadata=json.dumps(contrast.metadata()),rng_before=result['rng_before'],rng_after=result['rng_after'],checks=json.dumps(result['checks']),trace=json.dumps(result['trace']),last_accepted_time=result['last_accepted_time'],next_mark=np.nan if result['next_mark'] is None else result['next_mark'])
