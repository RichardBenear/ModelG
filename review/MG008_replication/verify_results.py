"""Compare independent arithmetic to saved evaluation, and reverify input provenance."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).parent
if __name__=='__main__':
    own=json.loads((ROOT/'independent_rescore.json').read_text());ref=json.loads((ROOT/'inputs/results/MG3b/evaluation.json').read_text())
    assert own['records']==258 and own['unique_keys']==258 and not own['bookkeeping_mismatches']
    assert not own['duplicate_run_ledger_keys'] and all(x['valid'] for x in own['stop20'].values())
    assert all(x['allowed'] for x in own['original_budget_replay'])
    for om,v in own['part2'].items():
        r=ref['part2']['cells'][om]
        assert v['events']==r['n'] and v['runs']==r['n_runs'] and v['excluded']==r['n_startup']
        assert abs(v['exposure']-r['exposure'])<1e-12 and abs(v['k']-r['k'])<1e-14
        assert np.allclose(v['ci'],[r['k_lo'],r['k_hi']],rtol=1e-12,atol=0)
        assert abs(v['stationarity']['p']-r['stationary']['p'])<1e-10
    for name in ['slope','slope_without100']:
        r=ref['part2']['slope' if name=='slope' else 'slope_without_100'];v=own[name]
        assert abs(v['S']-r['S'])<1e-10
        assert np.allclose(v['profile_ci'],[r['S_profile_lo'],r['S_profile_hi']],rtol=1e-9)
    assert np.allclose(own['part1']['ratio90']['ci'],ref['part1']['ci90'],rtol=1e-12)
    assert own['part1']['equivalent'] and ref['part1']['decision']=='EQUIVALENT'
    assert abs(own['ledger_total']-ref['ledger_core_h']['total'])<1e-10
    manifest=json.loads((ROOT/'download_manifest.json').read_text())
    for row in manifest['files']:
        data=(ROOT/'inputs'/row['path']).read_bytes()
        assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['sha'],row['path']
    out=dict(commit=manifest['commit'],input_hashes_verified=len(manifest['files']),independent_arithmetic_agrees=True)
    p=ROOT/'independent_confirmation.json'
    if p.exists():
        rows=json.loads(p.read_text());out.update(confirmation_fields=len(rows),confirmation_agreement=sum(r['agrees'] for r in rows))
    (ROOT/'verification.json').write_text(json.dumps(out,indent=2));print(out)
