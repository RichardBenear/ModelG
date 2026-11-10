"""Saved inventory recomputation using vectorized edge sums."""
from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/ring_feasibility'
def run():
 q=json.loads((BASE/'inventory.json').read_text());rows=[]
 for row in q['records']:
  with np.load(BASE/(row['id']+'.npz')) as z:p=z['counts'];omega=float(z['omega']);a0=z['a0'];b0=z['b0'];links=json.loads(str(z['links']));diff=z['diff'];k=z['coupling']
  totals=[]
  for s in range(5):
   kind='ActiveLink' if s<3 else ('ALink' if s==3 else 'BLink');edges=[x for x in links if x[0]==kind];src=np.array([x[1] for x in edges]);dst=np.array([x[2] for x in edges]);w=np.array([x[3] for x in edges]);flow=p[s,src]*w*(diff[s] if s<3 else k[s-3]);totals.append(float(flow.sum()));assert abs(np.bincount(dst,weights=flow,minlength=500).sum()-np.bincount(src,weights=flow,minlength=500).sum())<1e-6*max(1,flow.sum())
  g,x,y,a,b,_,_=p.astype(float);reaction=float((14*omega*a/a0+g+.1*x+29*b*x/b0+x*(x-1)*y/omega**2+.9*x).sum());np.testing.assert_allclose(totals,row['transport_rates'],rtol=1e-12,atol=1e-6);np.testing.assert_allclose(reaction+sum(totals),row['initial_rate'],rtol=1e-12);assert p.sum()==row['total_tokens'];rows.append(dict(id=row['id'],passed=True,initial_rate=reaction+sum(totals)))
 (BASE/'separate_recomputation.json').write_text(json.dumps(dict(passed=True,records=rows,shared='Saved inputs and scalar Model G formulas; vectorized edge arithmetic separate from inventory loop. Same author, not external replication.'),indent=2)+'\n')
if __name__=='__main__':run()
