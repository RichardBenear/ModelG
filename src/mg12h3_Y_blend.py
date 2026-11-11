"""MG151/152 fixed Y blend and predeclared fresh validation statistics."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import numpy as np,json,hashlib,time,math,sys,traceback
import mg12h3_finer_precision as core
from mg10e2c1_pool import atomic_json
ROOT=core.ROOT;BASE=ROOT/'results/MG12H3/Y_blend_validation';WEIGHTS=ROOT/'spec/MG12H3_Y_blend_weights.json';ALPHA=np.asarray(json.loads(WEIGHTS.read_text())['alpha_direct']);assert ALPHA.tolist()==[.31103637554093005,.620838372539336]
FILES=['src/mg12h3_Y_blend.py','src/prepare_mg12h3_Y_blend.py','src/run_mg12h3_Y_blend.py','spec/MG12H3_Y_blend_validation.md','spec/MG12H3_Y_blend_weights.json','results/MG12H3/finer_tick_precision_pilot/separate_recomputation.json','results/MG12H3/finer_tick_precision_pilot/production/ledger.json','results/MG12H3/contrast_estimator/single_channel.npz']+[f'results/MG12H3/contrast_estimator/mixed_c{c}.npz' for c in range(4)]
def hashes():
 h=core.hashes()
 for name in FILES:
  p=ROOT/name;raw=p.read_bytes() if p.suffix=='.npz' else p.read_text(encoding='utf-8').replace('\r\n','\n').encode();h[name]=hashlib.sha256(raw).hexdigest()
 return h

def blend(direct,integrated,time_index):
 a=ALPHA[time_index];return a*np.asarray(direct)+(1-a)*np.asarray(integrated)

def pooled_variance(array):return np.mean(np.var(array,axis=-1,ddof=1),axis=-1)

def validate_variance(direct,integrated,indices,time_index=1):
 direct=np.asarray(direct,float);integrated=np.asarray(integrated,float);indices=np.asarray(indices)
 if direct.shape!=(3,8) or integrated.shape!=(3,8) or not np.isfinite(direct).all() or not np.isfinite(integrated).all() or indices.ndim!=3 or indices.shape[1:]!=(3,8) or indices.dtype.kind not in 'iu' or indices.min()<0 or indices.max()>=8:raise ValueError('Incomplete/invalid paired input')
 mixed=blend(direct,integrated,time_index);vd=pooled_variance(direct);vh=pooled_variance(mixed);point=vh/vd if vd>0 else math.nan;per_tick=np.divide(mixed.var(1,ddof=1),direct.var(1,ddof=1),out=np.full(3,np.nan),where=direct.var(1,ddof=1)>0)
 a=direct[np.arange(3)[None,:,None],indices];b=mixed[np.arange(3)[None,:,None],indices];den=pooled_variance(a);numerator=pooled_variance(b);ratios=np.divide(numerator,den,out=np.full(len(indices),np.nan),where=den>0);finite=bool(np.isfinite(ratios).all() and np.isfinite(point) and np.isfinite(per_tick).all())
 if finite:
  lo,hi=np.percentile(ratios,[2.5,97.5]);slack=1e-12
  if hi<1-slack and np.all(per_tick<1-slack):verdict='VARIANCE BENEFIT SUPPORTED IN THIS VALIDATION'
  elif lo>1+slack:verdict='EVIDENCE OF HIGHER POOLED VARIANCE'
  elif hi<1-slack:verdict='POOLED BENEFIT; TICK-SPECIFIC VALIDATION UNRESOLVED'
  else:verdict='UNRESOLVED VARIANCE BENEFIT'
 else:lo=hi=math.nan;verdict='UNRESOLVED: NONFINITE OR ZERO VARIANCE'
 def number(x):return float(x) if np.isfinite(x) else None
 result=dict(time=(.005,.01)[time_index],alpha_direct=float(ALPHA[time_index]),pooled_direct_variance=float(vd),pooled_blend_variance=float(vh),blend_over_direct=number(point),per_tick_blend_over_direct=[number(x) for x in per_tick],bootstrap_95_interval=[number(lo),number(hi)],bootstrap_replicates=len(indices),nonfinite_bootstrap_count=int(np.count_nonzero(~np.isfinite(ratios))),verdict=verdict,coverage='Approximate stratified paired percentile bootstrap; small-sample coverage not guaranteed')
 return result,ratios
