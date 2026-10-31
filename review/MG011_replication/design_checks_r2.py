"""MG4c draft 2: detector and accounting counterexamples, not physics runs.
Only the Python standard library is required. No MG4/MG4c source is imported.
"""
from pathlib import Path
import json,math

POS=[(x,y) for y in range(9) for x in range(9)]

def distance(a,b):
    return math.hypot(*[min(abs(x-y),9-abs(x-y)) for x,y in zip(a,b)])

NEIGHBOURS=[[j for j,q in enumerate(POS) if distance(p,q)<=1] for p in POS]
SHAPE=[-math.exp(-distance(p,(4,4))**2/2) for p in POS]

def candidates(field,theta):
    return [i for i,nb in enumerate(NEIGHBOURS) if field[i]<theta and
            all(field[i]<field[j] or (field[i]==field[j] and i<=j) for j in nb)]

def decaying_well(final_depth):
    checks=[[final_depth*math.exp(.1*(300-t))*v for v in SHAPE]
            for t in range(200,301,10)]
    theta=.5*min(checks[-1])
    ids=[candidates(f,theta) for f in checks]
    # One stationary candidate at every time has no possible merge ambiguity;
    # its radius-.6 centroid uses just that node (grid spacing is 1), so links hold.
    single=all(len(x)==1 for x in ids) and len({x[0] for x in ids})==1
    fixed_ids=[candidates(f,-5.117354) for f in checks]
    return dict(final_depth=final_depth,theta=theta,counts=[len(x) for x in ids],
                class_with_self_normalized_threshold='single' if single else 'other',
                counts_with_fixed_MG4_threshold=[len(x) for x in fixed_ids],
                max_amplitude_at_first_check=abs(min(checks[0])))

def category(rows):
    if any(all(x=='single' for x in row) for row in rows):return 'W1'
    if any(x=='single' for row in rows for x in row):return 'W3'
    return 'W2'

def main():
    out=dict(reviewed_commit='877e9a6858aa69ff6678c96af1cc542bf8a3208a',
             proposal_blob='5ca1f4b30dffe88a2d769675b57256e3f3c9aa44',
             run_cap=8*5+2*6*5+2)
    out['threshold_scale_examples']=[decaying_well(e) for e in [1.0,1e-12,1e-30]]
    assert all(x['class_with_self_normalized_threshold']=='single'
               for x in out['threshold_scale_examples'])
    assert all(n==0 for x in out['threshold_scale_examples'][1:]
               for n in x['counts_with_fixed_MG4_threshold'])
    out['refined_four_single_now_covered']=category([['multiplied']*4,['decayed']*4,['single']*4])
    assert out['run_cap']==102 and out['refined_four_single_now_covered']=='W1'
    out['O2_counterexample']=dict(original='multiplied',weighted='decayed',
                                 same_outcome=False,single_recovered=False)
    out['launch_cutoff_counterexample']=dict(completed_at_launch=19.9,in_flight=9,
                                           hours_per_job=1.0,final_hours=19.9+9*1.0)
    assert out['launch_cutoff_counterexample']['final_hours']>22.5
    path=Path(__file__).with_suffix('.json')
    path.write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
