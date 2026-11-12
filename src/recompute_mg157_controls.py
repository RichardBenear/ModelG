import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,time,json,math
sys.path.insert(0,'src')
import numpy as np
from scipy.linalg import expm
from collections import defaultdict
from mg12h3_coupling_long import BASE,read,state
from mg10e2c1_pool import atomic_json
start=time.perf_counter();kernels={};count=0
for path in sorted(BASE.glob('control_kernel_*.npz')):
 with np.load(path) as z:
  q=np.zeros((len(z['remaining']),)*2);np.add.at(q,(z['q_rows'],z['q_cols']),z['q_rates']);prob=expm(float(z['delta'])*q)[0];np.testing.assert_allclose(prob,z['augmented_probabilities'],atol=1e-12,rtol=1e-10);law=defaultdict(float)
  for mass,out in zip(prob,z['outputs']):
   flat=np.zeros(7,np.int64)
   for code in out:site,sp=divmod(int(code),7);assert site==0;flat[sp]+=1
   law[tuple(flat)]+=float(mass)
  kernels[float(z['delta']),tuple(z['initial'])]=dict(law);count+=1
with np.load(BASE/'control_samples.npz') as z:samples=z['counts'];np.testing.assert_array_equal(z['seeds'],read(BASE/'plan.json')['seeds']['control'])
rows=[]
for side,d in enumerate((.01,.004)):
 law={tuple(state().counts().ravel()):1.}
 for step in range(round(.04/d)):
  nxt=defaultdict(float)
  for initial,mass in law.items():
   for final,p in kernels[d,initial].items():nxt[final]+=mass*p
  law=dict(nxt)
  if step+1 in (round(.02/d),round(.04/d)):
   f=1 if step+1==round(.02/d) else 2
   with np.load(BASE/f'control_law_s{side}_f{f}.npz') as z:np.testing.assert_allclose([law[tuple(k)] for k in z['support']],z['probabilities'],atol=1e-12,rtol=1e-10)
   freq=defaultdict(int)
   for pop in samples[:,side,f]:freq[tuple(pop.ravel())]+=1
   tv=sum(abs(freq.get(k,0)/len(samples)-v) for k,v in law.items())/2;bound=math.sqrt((len(law)*math.log(2)+math.log(4000))/(2*len(samples)));assert tv<=bound;rows.append(dict(side=side,frame=f,TV=tv,bound=bound))
atomic_json(BASE/'control_separate_recomputation.json',dict(passed=True,kernels=count,rows=rows,elapsed_seconds=time.perf_counter()-start,method='Dense saved-generator exponential, separate kernel composition and sample aggregation; same author, shared saved generator.'));print(json.dumps(rows))
