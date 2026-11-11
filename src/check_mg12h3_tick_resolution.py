"""MG144 independent column-kernel and exported spectral-error recomputation."""
from pathlib import Path
import numpy as np,json,hashlib,math
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/tick_resolution'
def run():
 q=json.loads((BASE/'results.json').read_text());rows=[]
 for rel,digest in q['source_hashes'].items():
  p=ROOT/rel;raw=p.read_bytes() if p.suffix=='.npz' else p.read_text(encoding='utf-8').replace('\r\n','\n').encode();assert hashlib.sha256(raw).hexdigest()==digest
 for path in sorted(BASE.glob('*.npz')):
  with np.load(path) as z:L=z['L'];exits=z['exits'];values=z['eigenvalues'];profiles=z['profiles'];continuous=z['continuous'];discrete=z['discrete'];deltas=z['deltas'];horizons=z['horizons'];saved=z['transition']
  case,species=path.stem.rsplit('_',1);rates=L+np.diag(exits)
  for d,delta in enumerate(deltas):
   # Columnwise heterogeneous formula, independent of regular-graph shortcut.
   factors=np.divide(-np.expm1(-exits*delta),exits,out=np.full_like(exits,delta),where=exits!=0);P=rates*factors[None,:]+np.diag(np.exp(-exits*delta));np.testing.assert_allclose(P,saved[d],rtol=1e-12,atol=1e-12)
   factor=float(factors[0]);steps=np.rint(horizons/delta).astype(int)
   for t,horizon in enumerate(horizons):
    base=1+factor*values;assert np.all(base>0);mult=np.exp(steps[t]*np.log1p(factor*values));reference=np.exp(values*horizon);operator_error=float(np.max(abs(mult-reference)));old=next(x for x in q['records'] if x['case']==case and x['species']==species and x['delta']==delta and x['horizon']==horizon);assert abs(operator_error-old['operator_L2_error'])<1e-10;errors=np.sqrt(np.einsum('ij,ij->j',discrete[d,t]-continuous[t],discrete[d,t]-continuous[t]));np.testing.assert_allclose(errors,old['profile_L2_error'],rtol=1e-10,atol=1e-12);assert max(errors)<=operator_error+1e-9;rows.append(dict(case=case,species=species,delta=float(delta),horizon=float(horizon),operator_L2_error=operator_error,actual_profile_L2_error=float(errors[-1]),passed=True))
 summary=[]
 for d in (.001,.0005,.0001,.00002):
  current=[r for r in rows if r['delta']==d];summary.append(dict(delta=d,worst_operator_L2_error=max(r['operator_L2_error'] for r in current),worst_actual_profile_L2_error=max(r['actual_profile_L2_error'] for r in current),screen_1percent_operator=all(r['operator_L2_error']<=.01 for r in current)))
 (BASE/'separate_recomputation.json').write_text(json.dumps(dict(passed=True,records=rows,summary=summary,source_hashes_unchanged=True,shared='Saved generator eigenvalues/means, independent columnwise transition construction and log1p exponential powers; same author, not external replication.'),indent=2)+'\n');print('elapsed',q['elapsed_seconds']);print(json.dumps(summary));print('files',len(list(BASE.glob('*.npz'))))
if __name__=='__main__':run()
