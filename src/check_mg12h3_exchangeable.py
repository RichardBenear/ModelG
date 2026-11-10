"""MG128 deterministic saved-support reservation generator check, no RNG."""
from prepare_mg12h3_bridge import ROOT,fixture,instantiate,config,matches,atomic_json,np
from mg12h3_exchangeable import channels,apply
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
    atomic_json(ROOT/'results/MG12H3/efficiency/algebra_check.json',dict(passed=True,states=checked,channels=total,random_draws=0,integrations=0,scope='Joint count increment and reservation rates, not stochastic executor qualification'))
if __name__=='__main__':run()
