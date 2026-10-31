"""Finite decision examples and budget arithmetic for MG4c draft 1. No model dynamics."""
from pathlib import Path
import json,math

S,M,D='single','multiplied','decayed'

def draft_categories(grid,refined):
    grid_all=any(all(v==S for v in row) for row in grid)
    refined_all=any(all(v==S for v in row) for row in refined)
    any_single=any(v==S for row in grid+refined for v in row)
    return dict(W1=grid_all,W2=not any_single,W3=any_single and not(grid_all or refined_all))

def main():
    out={}
    # T is required at every new b, including refinement, but not counted in its cap.
    out['runs']=dict(coarse=8*(1+4),refinement_D=2*6*4,
                     refinement_T=2*6,weighted=2,draft_total=90,corrected_max=40+48+12+2)
    out['refined_four_of_four_case']=draft_categories([[M]*4,[D]*4],[[S]*4])
    # D0 remains single at both endpoints of the proposed panel-edge bracket.
    out['invalid_D0_bisection_bracket']={'lower_b':28.5,'lower_classes':[S]*4,
        'upper_b':28.75,'upper_classes':[S,M,S,S],
        'D0_changes_predicate':False,'panel_changes_predicate':True}
    out['lowest_grid_single']={'b':27.0,'lower_neighbour_in_registered_grid':False}
    # Six binary decisions localize a monotone boundary only if it is bracketed.
    out['bracket_width_after_six_steps']=.25/2**6
    # A nonzero single-particle interval can lie between sampled dyadic points.
    lo,hi=28.0,28.25;window=(28.123,28.124);samples=[]
    for _ in range(6):
        b=(lo+hi)/2
        cls=D if b<window[0] else M if b>window[1] else S
        samples.append(dict(b=b,cls=cls))
        if cls==D:lo=b
        elif cls==M:hi=b
        else:break
    out['narrow_window_missed']={'actual_window':window,'samples':samples,
        'single_observed':any(x['cls']==S for x in samples)}
    # Verify the proposal's triangular normalization independently.
    a=.2/math.sqrt(1.5)
    area=math.sqrt(3)/2*a*a;weight=1/math.sqrt(3)
    out['triangular_weight_over_area']=weight/area
    out['MG4_triangular_rate']=1/.2**2
    assert out['runs']['corrected_max']==102
    assert not any(out['refined_four_of_four_case'].values())
    assert not out['narrow_window_missed']['single_observed']
    assert math.isclose(weight/area,1/.2**2,rel_tol=1e-14)
    (Path(__file__).parent/'design_checks.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
