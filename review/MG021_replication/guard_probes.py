"""Control-flow probes of MG021 authorization, with all physics mocked.
No dynamics, production files or upstream source are changed. Reports whether
main_amended would authorize a launch in deliberately invalid fixture cases.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import contextlib,io,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'inputs/src'))
import run_mg6_part0 as mod


def probe(case):
    old=json.loads((ROOT/'inputs/results/MG6/part0/decisions.json').read_text())
    final_state=dict(matches_manifest=None if case=='missing_manifest_entry' else True)
    checks=dict(bound_margin=dict(ok=True),transient=dict(window=[100.,200.]),final_state=final_state)
    timing_status='FAILED-NUMERICAL' if case=='failed_timing_run' else 'COMPLETE'
    wall=2000. if case=='projection_over_80' else 20.
    maximum=.09 if case=='timing_margin_exceeded' else .05
    timing_res=dict(status=timing_status,t_reached=1. if case=='failed_timing_run' else 20.,wall_s=wall,max_phiG=maximum)
    settle_status='FAILED-NUMERICAL' if case=='failed_settle' else 'COMPLETE'
    calls=[]
    class Pool:
        def __init__(self,*a,**kw):pass
        def __enter__(self):return self
        def __exit__(self,*a):pass
        def map(self,fn,args):
            args=list(args);calls.extend(args)
            return [dict(res=timing_res.copy()) for _ in args]
    scratch_parent=(ROOT/'inputs').resolve()
    assert scratch_parent.is_relative_to(ROOT.resolve())
    with TemporaryDirectory(prefix='_guard_probe_',dir=scratch_parent) as td:
        folder=Path(td)
        assert folder.resolve().is_relative_to(scratch_parent)
        (folder/'decisions.json').write_text(json.dumps(old))
        np.savez(folder/f'final_{mod.AMEND_ACC}.npz',pos=np.zeros((1,2)),F=np.zeros((3,1)))
        # These mocks replace physics / hardware timing only. The authorization
        # logic inside the actual pinned main_amended executes unchanged.
        with patch.object(mod,'P0',str(folder)),patch.object(mod.go,'continuum2d',lambda:None),\
             patch.object(mod,'amended_checks',lambda *a,**kw:checks.copy()),\
             patch.object(mod,'graph',lambda:(dict(N=old['N'],Lx=old['Lx'],Ly=old['Ly']),old['h'])),\
             patch.object(mod,'check_quadrature',lambda *a:dict(ok=True,refine_ok=case!='refined_grid_not_timed')),\
             patch.object(mod,'settle',lambda *a:(dict(status=settle_status,wall_s=1.),dict(pos=np.zeros((1,2)),F=np.zeros((3,1))))),\
             patch.object(mod,'ProcessPoolExecutor',Pool),patch.object(mod,'finish_amended',lambda d:d),\
             contextlib.redirect_stdout(io.StringIO()):
            result=mod.main_amended()
        return dict(case=case,proceed_written=result['proceed'],projection=result['timing']['projected_core_h_20_runs'],
            manifest_status=result['final_state']['matches_manifest'],timing_status=timing_status,
            timing_maximum=maximum,settle_status=settle_status,selected_grid=[result['nx'],result['ny']],
            timing_arguments_include_grid=any('nx' in a or 'ny' in a or 'grid' in a for a in calls))


def main():
    cases=['missing_manifest_entry','failed_timing_run','timing_margin_exceeded','failed_settle',
           'projection_over_80','refined_grid_not_timed']
    results=[probe(c) for c in cases]
    # This verifies the identified holes in the reviewed revision. These are
    # deliberately failing conditions, not successful MG6 calibration results.
    assert all(r['proceed_written'] for r in results)
    (ROOT/'guard_probes.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
