"""MG120 explicit provenance for22validated reuse slots and2fresh replays."""
from prepare_mg12h3_population import *
from run_mg12h3_tolerance import METHODS,diagnostic_hashes
BASE=ROOT/'results/MG12H3/completion';OLD=ROOT/'results/MG12H3/population';DIAG=ROOT/'results/MG12H3/tolerance_check'
FILES=('src/mg12h3_completion.py','src/prepare_mg12h3_completion.py','src/run_mg12h3_completion.py','src/check_mg12h3_completion.py','spec/MG12H3_completion_plan.md')
def matrix():
    complete={r['job']['id']:r for r in read(OLD/'ledger.json')['records']};entries=[]
    for job in jobs():
        fresh=job['kind']=='ISOLATED_C1' and job['start']==1
        if fresh:
            method=METHODS[job['setting']];path=BASE/'production'/(job['id']+'.npz');diagnostic=DIAG/(method['id']+'.npz')
        else:
            assert job['id'] in complete and complete[job['id']]['status']=='COMPLETE';method=dict(id='original_base' if job['setting']==0 else 'original_fine',solver='BDF',**SETTINGS[job['setting']]);path=OLD/'production'/(job['id']+'.npz');diagnostic=None
        entries.append(dict(job=job,origin='fresh_corrected_replay' if fresh else 'reused_validated_MG116',path=path.relative_to(ROOT).as_posix(),method=method,diagnostic_reference=diagnostic.relative_to(ROOT).as_posix() if diagnostic else None,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest() if not fresh else None))
    assert sum(e['origin']=='reused_validated_MG116' for e in entries)==22;return entries
def completion_hashes():
    result=dict(hashes(),**diagnostic_hashes(),**{name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in FILES})
    diagnostic_inputs=[DIAG/(m['id']+suffix) for m in METHODS for suffix in ('.npz','.json')]+[DIAG/name for name in ('registration.json','ledger.json','evaluation.json','raw_manifest.json')]
    for path in diagnostic_inputs:
        result[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes() if path.suffix=='.npz' else path.read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()
    return result
def verify_reuse(entries):
    validated={r['job']['id']:r for r in read(OLD/'ledger.json')['records']}
    for e in entries:
        if e['origin']!='reused_validated_MG116':continue
        p=ROOT/e['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==e['source_sha256'];assert validated[e['job']['id']]['status']=='COMPLETE' and validated[e['job']['id']]['check']['passed']
        expected=dict(id='original_base' if e['job']['setting']==0 else 'original_fine',solver='BDF',**SETTINGS[e['job']['setting']]);assert e['method']==expected
        with np.load(p) as z:assert json.loads(str(z['job_json']))==e['job']
    return True
def save_new(path,entry,times,state):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,state=state,times=times,job_json=json.dumps(entry['job']),method_json=json.dumps(entry['method']),origin='fresh_corrected_replay')
    os.replace(tmp,path)
