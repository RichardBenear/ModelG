"""Verify that draft 3 rejects the round-2 vanishing-amplitude examples.
Constructed fields only; no Model G dynamics. Python standard library only.
"""
from pathlib import Path
import json,math
from design_checks_r2 import SHAPE,candidates

def main():
    out=dict(reviewed_commit='e9c5e81cbc825512fe9cb1887d24929a862c5c38',
             proposal_blob='cee43bed8f7cbae7ad728e9d5dcc473e4227c309',T_eligibility=[])
    for epsilon in [1e-30,1e-12,1.0,10.0]:
        checks=[[epsilon*math.exp(.1*(300-t))*v for v in SHAPE] for t in range(200,301,10)]
        ids=[candidates(f,-1.0) for f in checks]
        fixed_single=all(len(x)==1 for x in ids) and len({x[0] for x in ids if x})==1
        floor=min(checks[-1])<=-2.0
        out['T_eligibility'].append(dict(final_depth=epsilon,
            fixed_threshold_single=fixed_single,depth_floor_pass=floor,
            eligible=fixed_single and floor))
    assert [x['eligible'] for x in out['T_eligibility']]==[False,False,False,True]
    out['run_cap']=40+60+2
    out['recommendation']='Support registration of Draft 3; no outstanding design blocker'
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
