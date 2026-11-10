"""MG120 no-integration completion preparation using held diagnostics."""
from mg12h3_completion import *
from check_mg12h3_completion import evaluate
def plan():
    if (BASE/'preparation.json').exists():raise ValueError('Already prepared')
    entries=matrix();verify_reuse(entries);ledger=read(DIAG/'ledger.json');timings={r['method']['id']:r['elapsed_seconds'] for r in ledger['records']};forecast=2*(timings['bdf_tight']+timings['bdf_tighter'])+30;assert forecast<=225
    atomic_json(BASE/'plan.json',dict(stage='MG120b',scope='planning/saved-evidence preparation only; fresh integrations await PI',source_original_hashes=completion_hashes(),entries=entries,forecast_seconds=forecast,elapsed_cap=300,workers=1,threads_per_library=1,original_failed_incomplete_preserved=True))
def prepare():
    plan=read(BASE/'plan.json');assert plan['source_original_hashes']==completion_hashes();verify_reuse(plan['entries']);start=time.perf_counter();held=[]
    for e in plan['entries']:
        copy=dict(e)
        if e['origin']=='fresh_corrected_replay':copy['path']=e['diagnostic_reference']
        held.append(copy)
    result=evaluate(held,BASE/'held_pipeline_evaluation.json',held=True);assert result['complete'] and result['qualified']
    # Missing second corrected slot must suppress matrix completion, not become zero.
    missing=[dict(e) for e in held];missing[-1]=dict(missing[-1],path=(BASE/'intentionally_absent.npz').relative_to(ROOT).as_posix());incomplete=evaluate(missing,BASE/'held_missing_evaluation.json',held=True);assert not incomplete['complete'] and not incomplete['qualified']
    # Altering a reused digest must be rejected before reading outcomes.
    corrupt=[dict(e) for e in held];corrupt[0]['source_sha256']='0'*64
    try:verify_reuse(corrupt)
    except AssertionError:pass
    else:raise AssertionError('Changed reused digest accepted')
    bad_method=[dict(e) for e in held];last=dict(bad_method[-1]);last['method']=dict(last['method'],rtol=1e-4);bad_method[-1]=last
    try:evaluate(bad_method,BASE/'must_not_pass.json',held=True)
    except AssertionError:pass
    else:raise AssertionError('Incorrect method provenance accepted')
    assert completion_hashes()==plan['source_original_hashes'];atomic_json(BASE/'preparation.json',dict(passed=True,plan_sha256=hashlib.sha256((BASE/'plan.json').read_bytes()).hexdigest(),elapsed_seconds=time.perf_counter()-start,scope='held old/diagnostic inputs only; no new integrations',source_original_hashes=completion_hashes(),complete_matrix_pass=True,missing_suppressed=True,reused_digest_tamper_rejected=True,method_tamper_rejected=True,forecast_seconds=plan['forecast_seconds']))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plan',action='store_true');p.add_argument('--prepare',action='store_true');a=p.parse_args()
    if a.plan:plan()
    elif a.prepare:prepare()
    else:p.error('Choose action')
