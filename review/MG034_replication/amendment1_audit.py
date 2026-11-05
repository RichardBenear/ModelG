"""Static/mocked MG9b Amendment 1 audit: python amendment1_audit.py CHECKOUT.
No production integration. The spawn probe only reads module attributes.
"""
import contextlib, io, json, multiprocessing, sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

def child_probe(_):
    import run_mg9b as d
    return {'simulation_function':d.simulate.__name__, 'output_folder':Path(d.OUT).name}

def audit(root):
    sys.path.insert(0,str(root/'src'))
    import run_mg9b as d, evaluate_mg9b as ev
    out={}
    with TemporaryDirectory() as tmp, contextlib.ExitStack() as stack:
        def replace(k,v): stack.enter_context(patch.object(d,k,v))
        replace('OUT',tmp); Path(tmp,'GO').write_text('synthetic test only')
        replace('Pool',lambda *_:contextlib.nullcontext(None))
        replace('stage_vb1',lambda:{29.35:{'pass_':True},29.20:{'pass_':False}})
        replace('part_a',lambda *_:{'Om_T':150}); replace('part_b',lambda *_:{'cells':'all complete'})
        def stop(*a): d.check_budget(1,'V-b3 29.35',10**8)
        replace('vb3',stop)
        with contextlib.redirect_stdout(io.StringIO()): d.stage_runs()
        bs=d.ledger()['budget_stop']; out['later_reference_already_invalid']={k:bs[k] for k in ('at','rate_done','rate_complete','budget_truncated')}
        out['later_reference_already_invalid']['expected_budget_truncated']=False
    recs=[dict(cell=cell,rep=0,first={'s_conf':80},event=True,confirmations=[dict(sc=80,n_persistent=n,n_deep_end=n)])
          for cell,n in [('b29.35_B_150',1),('b29.35_B_135',2)]]
    out['event_shapes']={'expected_single':1,'forward':ev.event_shapes(recs),'reversed':ev.event_shapes(recs[::-1])}
    with TemporaryDirectory() as tmp,patch.object(d,'OUT',tmp),contextlib.redirect_stdout(io.StringIO()):
        r=d.vb3(None,29.35,150);out['empty_tolerance_control']={'n':r['n'],'result':r['result']}
    with patch.object(d,'OUT','mock_test_output'),patch.object(d,'simulate',lambda *_:None):
        with multiprocessing.get_context('spawn').Pool(1) as pool: out['fresh_worker']=pool.map(child_probe,[None])[0]
    c=json.loads((root/'results/MG9b/calibration/calibration.json').read_text()); labels={}
    for br in ('B3','B2-20','B2-35'):
        support=[(x['S1_true'],x['branches'][br]) for x in c['cases'] if x['S1_true'] is not None and x['branches'][br]['n']>=400]
        eligible={s==0 for s,_ in support}=={True,False}
        passed=eligible and all(b['coverage']['k']/b['n']>=.93 and (s!=0 or all(b[k]['k']/b['n']<=.035 for k in ('neg','pos'))) for s,b in support)
        labels[br]='calibrated' if passed else 'not calibrated'
        assert labels[br]==c['labels'][br]['label']
    assert all(x['branch_freq']['B2-other']==0 for x in c['cases'])
    out['calibration_labels_from_summary_counts']=dict(labels,**{'B2-new':'not calibrated (unsupported)'})
    return out

if __name__=='__main__':
    data=audit(Path(sys.argv[1]).resolve()); s=json.dumps(data,indent=2);print(s); Path(__file__).with_suffix('.json').write_text(s+'\n')
