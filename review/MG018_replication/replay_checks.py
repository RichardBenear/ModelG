"""Regenerate registered event lists, then independently audit tick counts and
sparse node lifetimes/positions. Uses the pinned generator for canonical replay;
does not integrate fields, build evolving meshes, or save large event lists.
"""
from pathlib import Path
import json, sys, hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'inputs/src'))
from evolve import gen_events


def main():
    out={};data=ROOT/'inputs/results/MG5'
    for f in sorted((data/'runs').glob('*.json')):
        rec=json.loads(f.read_text())
        if rec['nu']==0:continue
        key=rec['key'];rp=rec['replay'];sk=rec['spec']['seed']
        spawn=tuple(sk[n] for n in ['gate','b_idx','om_idx','stage','replicate'])+(0,)
        seed=int(np.random.SeedSequence(20261027,spawn_key=spawn).generate_state(1,dtype=np.uint32)[0])
        assert seed==rp['seed']
        ev=gen_events(np.arange(18901),18901,rec['nu'],rec['t_end'],rec['Lx'],rec['Ly'],seed)
        h=hashlib.sha256()
        for name,dt in [('t','<f8'),('type','i1'),('id','<i8'),('x','<f8'),('y','<f8')]:
            h.update(np.asarray(ev[name],dtype=dt).tobytes())
        digest=h.hexdigest()
        assert digest==rp['sha256']==rp['fresh_process_sha256'] and len(ev)==rp['n_events']
        born=ev['type']==1;dead=~born;bid=ev['id'][born];did=ev['id'][dead]
        assert np.array_equal(bid,np.arange(18901,18901+len(bid)))
        assert len(np.unique(did))==len(did) and np.all(np.diff(ev['t'])>0)
        limit=18901+len(bid)
        bt=np.full(limit,-np.inf);dt=np.full(limit,np.inf)
        bt[bid]=ev['t'][born];dt[did]=ev['t'][dead]
        assert np.all(dt>=bt)
        with np.load(data/'checks'/f'{key}.npz') as z:d={k:z[k] for k in z.files}
        delta=rec['delta'];K=round(rec['t_end']/delta)
        # Tick k contains (k*delta,(k+1)*delta]; searchsorted-left assigns
        # exact boundary events to the preceding tick, as registered.
        tick=np.searchsorted(np.arange(1,K+1)*delta,ev['t'],side='left')
        pb=np.bincount(tick[born],minlength=K);pd=np.bincount(tick[dead],minlength=K)
        birth_tick=np.full(limit,-1,dtype=int);birth_tick[bid]=tick[born]
        eph=np.bincount(tick[dead][birth_tick[did]==tick[dead]],minlength=K)
        for col,v in [('n_prop_births',pb),('n_prop_deaths',pd),('n_ephemeral',eph),
                      ('n_eff_births',pb-eph),('n_eff_deaths',pd-eph)]:
            assert np.array_equal(v[:-1],d['tick_'+col]),(key,col)
        max_newborn_position_error=0.
        for i,t in enumerate(d['t']):
            a,b=map(int,d['off'][i:i+2]);ids=d['ids'][a:b]
            boundary=0. if t==0 else (round(t/delta)-1)*delta
            assert np.all((ids>=0)&(ids<limit))
            assert np.all(bt[ids]<=boundary) and np.all(dt[ids]>boundary)
            is_new=ids>=18901
            if np.any(is_new):
                index=ids[is_new]-18901
                actual=np.column_stack([ev['x'][born][index],ev['y'][born][index]])
                err=float(np.max(np.abs(actual-d['pos'][a:b][is_new])))
                assert err==0.;max_newborn_position_error=max(max_newborn_position_error,err)
        out[key]=dict(events=len(ev),sha256=digest,seed=seed,hash_matches=True,
            all_tick_counts_match=True,all_sparse_node_lifetimes_valid=True,
            max_newborn_position_error=max_newborn_position_error,
            final_interval_events_deliberately_unapplied=int(pb[-1]+pd[-1]))
        print(key,len(ev),'hash / ticks / lifetimes / positions agree',flush=True)
    assert len(out)==16
    (ROOT/'replay_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('All 16 event lists reproduce. Total events:',sum(r['events'] for r in out.values()))


if __name__=='__main__':main()
