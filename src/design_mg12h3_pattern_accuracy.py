"""MG145 delegated pattern-accuracy choice, saved-data planning only."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import numpy as np,json,math,hashlib
from scipy.stats import t,norm
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/finer_tick_design'
def run():
 p=ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz'
 with np.load(p) as z:counts=z['counts'];reference=z['concentrations']
 scales=np.array([140/9,140/9,29/(140/9)]);records=[]
 for s,name in enumerate(('G','X','Y')):
  template=reference[s]/scales[s];template=template-template.mean();rms=float(np.sqrt(np.mean(template**2)));assert rms>0;template/=rms;initial=float(np.mean(template*counts[s]/10/scales[s]));assert initial>0;margin=.05*initial;weight=template/(500*10*scales[s]);rows=[]
  with np.load(ROOT/f'results/MG12H3/tick_resolution/o0_a0_Omega1_{name}.npz') as z:mu=z['eigenvalues'];U=z['eigenvectors'];lam=float(z['exits'][0])
  for horizon in (.005,.01):
   variances={}
   for delta in (.0001,.00005,.00002):
    factor=-math.expm1(-lam*delta)/lam;K=(U*((1+factor*mu)**round(horizon/delta)))@U.T;n=counts[s].astype(float);variances[delta]=float(n@(K.T@(weight**2)-(K.T@weight)**2))
   for delta in (.0001,.00005):
    v1,v2=variances[delta],variances[.00002];df=47*(v1+v2)**2/(v1*v1+v2*v2);half=float(t.ppf(1-.05/48,df)*math.sqrt((v1+v2)/48));approx_r=(norm.ppf(1-.05/48)**2)*(v1+v2)/margin**2;rows.append(dict(horizon=horizon,delta=delta,transport_only_R48_halfwidth=half,margin=margin,halfwidth_to_margin=half/margin,large_R_width_match_approximation=float(approx_r)))
  records.append(dict(species=name,template_basis='Unrounded frozen prepared concentration, centered and unit spatial RMS; avoids using rounding artifacts as the intended pattern',initial_integer_state_contrast=initial,margin_5percent=margin,planning=rows))
 out=dict(storage_approved=True,accuracy_authorization='PI delegated accuracy choice; lead selects 5% of initial prepared pattern contrast, not 1% of background concentration',records=records,scope='Transport-only precision planning; not coupled variance or a production sample-size determination',random_draws=0,source_hash=hashlib.sha256(p.read_bytes()).hexdigest());(BASE/'pattern_accuracy_choice.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(records))
if __name__=='__main__':run()
