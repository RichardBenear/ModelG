"""Replay a short run with two excursion labels; retain their fields for independent continuation."""
import sys,json
from pathlib import Path
import numpy as np
from rescore import load
from independent_confirmation import confirm
ROOT=Path(__file__).parent;INPUT=ROOT/'inputs'
sys.path.insert(0,str(INPUT/'src'))
import mg3
if __name__=='__main__':
    r=load(INPUT/'results/MG3/runs/B_100/13.pkl')
    st=mg3.sim_run(mg3.new_state(r['b'],r['Om'],r['tau'],r['key'],r['b_det']))
    out=dict(cell='B_100',rep=13,claimed_births=r['births'],replayed_births=st['births'],first_field_equal=np.array_equal(r['first']['field'],st['first']['field']),t_end=st['t'])
    own=[]
    for b in st['births']:
        phi=st['phi'][b['s_conf']].astype(float)
        own.append(dict(t_start=b['t_start'],claimed=b['cls'],**confirm(phi,r['b'])))
        np.savez(ROOT/f"replay_birth_{b['t_start']:g}.npz",phi=phi)
    out['independent_outcomes']=own
    (ROOT/'replay.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
