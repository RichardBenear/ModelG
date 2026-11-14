"""Independent-reference scalar noise calibration and bounded acquisition."""
import numpy as np
from scipy.stats import chi2,beta,f
from mg174_modes import ModeTracker
BANDS=np.array_split(np.arange(1,251),8)
def estimate(model,L,fields,seed):
 base=ModeTracker(model,L,0,seed);ys=np.fft.rfft(fields-fields.mean(axis=2)[:,:,None],axis=2);q=base.inner(ys,ys);total=float(q.sum());df=16*998;bands=[];okay=True
 for ix in BANDS:
  mask=np.zeros(251);mask[ix]=1;qb=float(np.sum(np.einsum('rsf,fst,rtf,f->r',np.conj(ys),base.P,ys,base.freqweights*mask).real/500));db=int(16*2*base.freqweights[ix].sum());lo,hi=beta.ppf([.001/16,1-.001/16],db/2,(df-db)/2);passed=bool(lo<=qb/total<=hi);okay &= passed;bands.append(dict(q=qb,df=db,lo=float(lo),hi=float(hi),passed=passed))
 scale=total/df;ci=total/chi2.ppf([.975,.025],df);return dict(scale=scale,df=df,q=q.tolist(),bands=bands,adequate=bool(okay),ci95=ci.tolist())
def acquire(records):
 for t in range(1,6):
  if records[t-1]['emitted']=='single' and records[t]['emitted']=='single':return t
 return None
def outcome(records,truth,positive,absent):
 acq=acquire(records)
 if positive:
  tail=records[acq+1:] if acq is not None else [];tr=truth[acq+1:] if acq is not None else [];n=len(tail);single=sum(r['emitted']=='single' for r in tail);accurate=sum(r['emitted']=='single' and min((r['center']-g[0])%500,(g[0]-r['center'])%500)<=2 and abs(r['width']-g[1])<=.12500001 for r,g in zip(tail,tr));success=acq is not None and single>=np.ceil(.9*n) and accurate>=np.ceil(.9*n);return dict(acquisition=acq,post_frames=n,single=single,accurate=accurate,success=bool(success))
 false_acq=any(records[t-1]['emitted']=='single' and records[t]['emitted']=='single' for t in range(1,len(records))) if absent else None;false=records[15]['emitted']=='single' if not absent else None;return dict(acquisition=acq,false_acquisition=false_acq,false_single_deadline=false,success=bool(not false_acq if absent else acq is not None and not false))
def tracker(model,L,rho,seed,cal):
 obj=ModeTracker(model,L*np.sqrt(cal['scale']),rho,seed);obj.qcut=float(998*f.ppf(.999,998,cal['df']));return obj
