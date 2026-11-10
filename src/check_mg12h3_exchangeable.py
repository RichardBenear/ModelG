"""MG128 deterministic saved-support reservation generator check, no RNG."""
from prepare_mg12h3_bridge import ROOT,fixture,instantiate,config,matches,atomic_json,np
from mg12h3_exchangeable import channels,apply,tick,simulate
from collections import defaultdict

def key(rule,change,reserved):return rule,tuple(change.ravel()),tuple(reserved.ravel())
def run():
    checked=0;total=0
    for c in range(4):
        state=fixture(c)
        with np.load(ROOT/f'results/MG12H3/bridge/c{c}_reference.npz') as z:states=z['states']
        for flat in states:
            a=flat.reshape(7,len(state.sites));explicit=defaultdict(float);grouped=defaultdict(float)
            for event in matches(instantiate(state,a)):
                change=np.zeros_like(a);reserved=np.zeros_like(a)
                for token,site,sp in event.before:
                    s=('G','X','Y','A','B','Z','W').index(sp);reserved[s,site]+=1
                before={token:(site,sp) for token,site,sp in event.before}
                for token,site,sp in event.after:
                    oldsite,oldsp=before[token];change[('G','X','Y','A','B','Z','W').index(oldsp),oldsite]-=1;change[('G','X','Y','A','B','Z','W').index(sp),site]+=1
                explicit[key(event.rule,change,reserved)]+=event.rate
            for channel in channels(a,config(state)):
                rule,change,reserved,rate=channel;grouped[key(rule,change,reserved)]+=rate;remaining,result=apply(a,a,channel);np.testing.assert_array_equal(remaining,a-reserved);np.testing.assert_array_equal(result,a+change);total+=1
            assert explicit.keys()==grouped.keys()
            for k in explicit:np.testing.assert_allclose(explicit[k],grouped[k],rtol=1e-12,atol=1e-12)
            checked+=1
    class Fixed:
        def exponential(self,scale):return .01
        def random(self):return 0.
    # Catalysts must be unavailable after one event; products stay unavailable.
    a=np.zeros((7,1),np.int64);a[1,0]=2;a[2,0]=2
    cfg=(2.,[(10.,4.)],[],(.2,.3,.4),(.4,.7),(False,False,False,False,True,False))
    out=tick(a,cfg,.1,Fixed());assert len(out['trace'])==1 and out['counts'][1,0]==3 and out['available'][1,0]==0
    # A->G followed by G->X is forbidden within the same tick.
    a=np.zeros((7,1),np.int64);a[3,0]=1
    cfg=(2.,[(10.,4.)],[],(.2,.3,.4),(.4,.7),(True,True,False,False,False,False))
    out=tick(a,cfg,.1,Fixed());assert len(out['trace'])==1 and out['counts'][0,0]==1 and out['counts'][1,0]==0
    assert not tick(a,cfg,.01,Fixed())['trace']
    # Transport arrival cannot move again until the next tick.
    a=np.zeros((7,2),np.int64);a[0,0]=1
    cfg=(2.,[(10.,4.)]*2,[('ActiveLink',0,1,1.),('ActiveLink',1,0,1.)],(1.,1.,1.),(1.,1.),(False,)*6)
    out=tick(a,cfg,.1,Fixed());assert len(out['trace'])==1 and out['counts'][0,1]==1
    out=simulate(a,cfg,[0.,.025,.05],0,delta=.1,rng=Fixed());assert out['steps']==2 and len(out['trace'])==2 and out['counts'][-1,0,0]==1
    zero=np.zeros((7,2),np.int64);assert not simulate(zero,cfg,[0.,.025],0,rng=Fixed())['trace']
    try:list(channels(-np.ones((7,2),np.int64),cfg));raise AssertionError('Negative accepted')
    except ValueError:pass
    atomic_json(ROOT/'results/MG12H3/efficiency/algebra_check.json',dict(passed=True,states=checked,channels=total,random_draws=0,integrations=0,scripted_tick_checks=7,scope='Joint count increment and reservation rates, not stochastic executor qualification'))
if __name__=='__main__':run()
