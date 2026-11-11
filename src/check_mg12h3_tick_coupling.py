"""MG156 saved-only recomputation, no random draws or coupled execution."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import math,json
import numpy as np
from scipy.linalg import expm
from run_mg12h3_tick_coupling import BASE,ROOT,read,hashes
from prepare_mg12h3_bridge import fixture,config
from mg12h3_sparse import Model
from mg12h3_tick_coupling import replay
from mg12h3_finer_precision import setup
from mg10e2c1_pool import atomic_json
def check():
    plan=read(BASE/'plan.json');assert hashes()==plan['source_hashes'];a=read(BASE/'stage_A_results.json');marginals=[];sentinels=0;aggregates=0;transport=[]
    for c in range(4):
        state=fixture(c);model=Model(len(state.sites),config(state));left=state.counts();right=left.copy();right[1,0]-=1;right[1,1]+=1
        with np.load(BASE/f'marginal_c{c}.npz') as z:
            counts=z['counts'];hist=z['histograms'];np.testing.assert_array_equal(z['seeds'],plan['seeds']['marginal'][c])
        for s,initial in enumerate((left,right)):
            with np.load(BASE/f'exact_c{c}_s{s}.npz') as z:
                q=np.zeros((len(z['remaining']),)*2);np.add.at(q,(z['q_rows'],z['q_cols']),z['q_rates']);assert abs(q.sum(axis=1)).max()<1e-10;prob=expm(.1*q)[0];np.testing.assert_allclose(prob,z['augmented_probabilities'],rtol=1e-10,atol=1e-12)
                projected={}
                for mass,out in zip(prob,z['outputs']):
                    flat=np.zeros(initial.size,np.int64)
                    for code in out:
                        site,sp=divmod(int(code),7);flat[sp*model.n+site]+=1
                    key=tuple(flat);projected[key]=projected.get(key,0.)+float(mass)
                exact=np.array([projected[tuple(k)] for k in z['support']]);np.testing.assert_allclose(exact,z['probabilities'],rtol=1e-10,atol=1e-12);observed={}
                for pop in counts[:,s]:key=tuple(pop.ravel());observed[key]=observed.get(key,0)+1
                assert set(observed)<=set(projected);tv=sum(abs(observed.get(k,0)/len(counts)-v) for k,v in projected.items())/2;bound=math.sqrt((len(projected)*math.log(2)+math.log(8000))/(2*len(counts)));assert tv<=bound;marginals.append(dict(cell=c,side=s,TV=tv,bound=bound))
            increments=np.zeros((len(counts),initial.size),np.int64)
            for channel,ch in enumerate(model.channels):
                accepted=hist[:,0,channel]+hist[:,s+1,channel]
                for coordinate,change in ch.change:increments[:,coordinate]+=accepted*change
            np.testing.assert_array_equal(initial.ravel()+increments,counts[:,s].reshape(len(counts),-1));aggregates+=len(counts)
        for prefix,init,deltas,times in ((f'marginal_c{c}_r0',(left,right),(.1,.1),[0,.1]),(f'regression_c{c}_r0',(left,left),(.01,.01),[0,.1])):
            with np.load(BASE/(prefix+'.npz')) as z:r={k:z[k] for k in z.files}
            r.update(read(BASE/(prefix+'_trace.json')));replay(model,*init,deltas,times,r);sentinels+=1
    with np.load(BASE/'transport.npz') as z:counts=z['counts'];np.testing.assert_array_equal(z['seeds'],plan['seeds']['transport'])
    for s,d in enumerate((.005,.002)):
        # Separate two-state transition matrix power, instead of closed-form scalar power.
        stay=math.exp(-80*d);matrix=np.array([[stay,1-stay],[1-stay,stay]])
        for frame,T in enumerate((.01,.02),1):
            expected=float(np.linalg.matrix_power(matrix,round(T/d))[0,1]);observed=float(counts[:,s,frame,0,1].sum()/len(counts));bound=math.sqrt(math.log(8000)/(2*len(counts)));assert abs(observed-expected)<=bound;transport.append(dict(side=s,time=T,expected=expected,observed=observed,bound=bound))
    for saved,fresh in zip(a['exact_marginals'],marginals):assert abs(saved['TV']-fresh['TV'])<1e-12
    for saved,fresh in zip(a['transport'],transport):assert abs(saved['expected']-fresh['expected'])<1e-12 and saved['observed']==fresh['observed']
    actual=[]
    if (BASE/'stage_B_results.json').exists():
        b=read(BASE/'stage_B_results.json');initial,contrast=setup();model=contrast.model
        for c,deltas in enumerate(((.0001,.00002),(.00005,.00002))):
            coupled=[];independent=[]
            for r in range(4):
                with np.load(BASE/f'actual_c{c}_r{r}.npz') as z:counts=z['counts'];hist=z['histogram'];record={k:z[k] for k in z.files}
                for s in range(2):
                    reconstructed=initial.ravel().copy();accepted=hist[0]+hist[s+1]
                    for n,ch in enumerate(model.channels):
                        for coordinate,change in ch.change:reconstructed[coordinate]+=int(accepted[n])*change
                    np.testing.assert_array_equal(reconstructed,counts[s,-1].ravel());assert sum(map(int,reconstructed))==sum(map(int,initial.ravel()));aggregates+=1
                if r==0:record.update(read(BASE/f'actual_c{c}_r{r}_trace.json'));replay(model,initial,initial,deltas,[0,.0001],record);sentinels+=1
                coupled.append([contrast.weights@x.ravel() for x in counts[:,-1]]);sides=[]
                for s in range(2):
                    with np.load(BASE/f'actual_single_c{c}_r{r}_s{s}.npz') as z:assert int(z['seed'])==plan['seeds']['actual_single'][c][2*r+s];sides.append(contrast.weights@z['counts'].ravel())
                independent.append(sides)
            coupled=np.asarray(coupled);independent=np.asarray(independent);dc=coupled[:,0]-coupled[:,1];di=independent[:,0]-independent[:,1]
            # Centered sums of squares, separate from driver's np.var computation.
            vc=np.sum((dc-dc.mean(axis=0))**2,axis=0)/3;vi=np.sum((di-di.mean(axis=0))**2,axis=0)/3;ratio=vc/vi
            np.testing.assert_allclose(ratio,b['comparisons'][c]['ratio'],rtol=1e-10,atol=1e-12);actual.append(dict(comparison=c,ratio=ratio.tolist(),mean_difference=dc.mean(axis=0).tolist()))
    result=dict(passed=True,source_hashes=hashes(),exact_marginals=marginals,transport=transport,aggregate_chain_endpoints=aggregates,sentinel_replays=sentinels,actual=actual,limits='Same author; shared frozen rule/contrast definitions and retained traces. Dense saved-generator exponential and separate count aggregation, not external replication.');atomic_json(BASE/'separate_recomputation.json',result);print(json.dumps(result))
if __name__=='__main__':check()
