"""Run upstream MG4c tests with workspace-local temporary files; report Windows skips explicitly."""
from pathlib import Path
import os,sys,tempfile,importlib.util,json,time,traceback
ROOT=Path(__file__).resolve().parent
temp=ROOT/'tmp';temp.mkdir(exist_ok=True)
tempfile.tempdir=str(temp)
sys.path.insert(0,str(ROOT/'inputs/src'))
result=[]
for filename in ['test_mg4c.py']:
    spec=importlib.util.spec_from_file_location(filename[:-3],ROOT/'inputs/tests'/filename)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    for name,test in list(vars(module).items()):
        if not name.startswith('test_') or not callable(test):continue
        if os.name=='nt' and name in ['test_pipeline_synthetic_case_A','test_pipeline_synthetic_case_B_and_W2']:
            row=dict(test=name,status='skipped',reason='upstream fork-dependent test skips Windows')
        else:
            start=time.perf_counter()
            try:test();row=dict(test=name,status='passed',seconds=time.perf_counter()-start)
            except Exception:row=dict(test=name,status='failed',traceback=traceback.format_exc())
        result.append(row);print(row,flush=True)
        (ROOT/'test_results.json').write_text(json.dumps(result,indent=2))
(ROOT/'test_results.json').write_text(json.dumps(result,indent=2))
if any(row['status']=='failed' for row in result):raise SystemExit(1)
