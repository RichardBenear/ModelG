"""MG5 draft 1: event-work arithmetic, remap example, and static baseline RHS.
No MG5 trajectories or rewriting dynamics are run.
"""
from pathlib import Path
import math,json,sys
ROOT=Path(__file__).resolve().parent
IN=ROOT/'inputs'

def main():
    n=18901
    out=dict(commit='3179342a9b559af06beabdfea98890c627a97e12',proposal_blob='9d6dd2b6800ffcbcef61923a232180c0f341c944',rates=[])
    for nu in [.002,.02,.2]:
        epoch=min(1,.05/nu);horizon=500
        expected_each=n*nu*horizon
        out['rates'].append(dict(nu=nu,epoch=epoch,epochs=horizon/epoch,
            expected_births=expected_each,expected_deaths=expected_each,
            expected_events_per_epoch=2*n*nu*epoch,
            expected_events_per_run=2*expected_each,
            full_rebuild_hours_if_one_second_per_event=2*expected_each/3600))
    out['birth_death_count_sd_limit']=math.sqrt(n)
    out['halved_fast_epoch_steps_at_dt005']=(.25/2)/.05
    out['quadratic_interpolation_example']=dict(triangle_vertices=[[-1,0],[1,0],[0,1]],
        birth=[0,1/3],vertex_values=[1,1,1],interpolated_value=1,
        actual_x2_plus_y2=1/9,bias=8/9)
    out['PartB_background_shift']=dict(X0=14/.9,Y0_28=28/(14/.9),Y0_29=29/(14/.9),
        phiY_adjustment_to_preserve_absolute_Y=-1/(14/.9))
    assert math.isclose(out['rates'][-1]['expected_events_per_run'],3780200)
    assert out['halved_fast_epoch_steps_at_dt005']==2.5
    if IN.exists():
        import numpy as np
        sys.path.insert(0,str(IN/'src'))
        import run_mg4c as rm,graph_ode as go
        dec=json.loads((IN/'results/MG4/part0/decisions.json').read_text())
        out['stored_initial_state_stationarity']=[]
        for slot,pl in [(0,0),(0,1),(1,0),(1,1)]:
            j=rm.job(28,'D',slot,pl);g,h,c=rm.graph_for(j,dec)
            y=np.load(IN/f"results/MG4c/fields/{j['key']}.npz")['final'].ravel()
            P={**go.P17,'b':28}
            old=go.make_system(g,h,P)[0](0,y)
            changed=go.make_system(g,.2,P)[0](0,y)
            out['stored_initial_state_stationarity'].append(dict(key=j['key'],MG4c_h=h,
                max_rhs_at_MG4c_h=float(abs(old).max()),max_rhs_at_MG5_fixed_h=float(abs(changed).max())))
    (ROOT/'design_checks.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
