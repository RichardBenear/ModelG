"""MG135 additional no-RNG cross-engine tick and ring inventory checks."""
from check_mg12h3_sparse import *
from mg12h3_exchangeable import tick as dense_tick
import sys,hashlib
class Fixed:
    def exponential(self,scale):return .2*scale
    def random(self):return .37123456

def deep_size(obj,seen=None):
    seen=set() if seen is None else seen
    if id(obj) in seen:return 0
    seen.add(id(obj));size=sys.getsizeof(obj)
    if isinstance(obj,dict):size+=sum(deep_size(k,seen)+deep_size(v,seen) for k,v in obj.items())
    elif isinstance(obj,(tuple,list,set)):size+=sum(deep_size(v,seen) for v in obj)
    elif hasattr(obj,'__dict__'):size+=deep_size(vars(obj),seen)
    return size

def run():
    ticks=[]
    for c in range(4):
        state=fixture(c);cfg=config(state);model=Model(len(state.sites),cfg);pop=state.counts()
        for step in range(5):
            a=tick(model,pop,.1,Fixed());b=dense_tick(pop,cfg,.1,Fixed());np.testing.assert_array_equal(a['counts'],b['counts']);np.testing.assert_array_equal(a['available'],b['available']);assert a['events']==len(b['trace']);pop=a['counts']
        ticks.append(dict(fixture=c,ticks=5,passed=True))
    inventory=json.loads((ROOT/'results/MG12H3/ring_feasibility/inventory.json').read_text());memory=[]
    for row in inventory['records']:
        with np.load(ROOT/f'results/MG12H3/ring_feasibility/{row["id"]}.npz') as z:pop=z['counts'];cfg=(float(z['omega']),list(zip(z['a0'],z['b0'])),json.loads(str(z['links'])),tuple(z['diff']),tuple(z['coupling']),(True,)*6)
        model=Model(500,cfg);work=model.workspace(pop);np.testing.assert_allclose(math.fsum(work.tree.values),row['initial_rate'],rtol=1e-12);memory.append(dict(case=row['id'],estimated_python_object_bytes=deep_size((model,work)),initial_rate_verified=True))
    sources={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in ['src/mg12h3_sparse.py','src/check_mg12h3_sparse.py','src/check_mg12h3_sparse_post.py','spec/MG12H3_sparse_development.md']}
    (BASE/'postcheck.json').write_text(json.dumps(dict(passed=True,scripted_cross_engine_ticks=ticks,memory_estimates=memory,source_hashes=sources,random_draws=0,limits='Python object size traversal is an estimate, not process peak RSS. Deterministic injected clocks, not sampled trajectories or statistical qualification.'),indent=2)+'\n')
if __name__=='__main__':run()
