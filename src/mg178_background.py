"""Independent spatial reference and separate temporal reference calibration."""
import numpy as np
from scipy.stats import beta,f,norm
from mg176_acquisition import estimate
from mg174_modes import ModeTracker
from mg178_persistent import PersistentTracker
ALPHA=.001/12

def calibrate(model,L,spatial,temporal,seed):
 cal=estimate(model,L,spatial,seed);base=ModeTracker(model,L,0,seed);ys=np.fft.rfft(spatial-spatial.mean(axis=2)[:,:,None],axis=2);yt=np.fft.rfft(temporal-temporal.mean(axis=2)[:,:,None],axis=2);den=float(base.inner(yt[:-1],yt[:-1]).sum());rho=float(base.inner(yt[1:],yt[:-1]).sum()/den);valid=abs(rho)<.99;used=float(np.clip(rho,-.99,.99));variance=1-used**2;innov=(yt[1:]-used*yt[:-1])/np.sqrt(variance);qtime=float(base.inner(innov,innov).sum());df=63*998-1;scale=cal['scale'];se=float(np.sqrt(qtime*variance/df/den));cal['rho']=used;cal['raw_rho']=rho;cal['rho_ci95']=[rho-1.96*se,rho+1.96*se];ev,V=np.linalg.eigh(base.P);W=np.einsum('fij,fj,fkj->fik',V,np.sqrt(np.maximum(ev,0)),np.conj(V));white=np.einsum('fab,rbf->raf',W,ys);energy=np.sum(base.freqweights[None,None,:]*abs(white)**2,axis=(0,2))/500;fraction=float(energy[0]/energy.sum());pc=float(min(1.,2*min(beta.cdf(fraction,3992,3992),beta.sf(fraction,3992,3992))));mean=ys.mean(axis=0);qm=float(16*base.inner(mean,mean));qc=float(base.inner(ys-mean,ys-mean).sum());pm=float(f.sf(15*qm/qc,998,14970));ratio=qtime/df/scale;pv=float(min(1.,2*min(f.cdf(ratio,df,15968),f.sf(ratio,df,15968))));lag=float(base.inner(innov[2:],innov[:-2]).sum()/np.sqrt(base.inner(innov[2:],innov[2:]).sum()*base.inner(innov[:-2],innov[:-2]).sum()));pl=float(2*norm.sf(abs(lag)*np.sqrt(61*998)));bands=[]
 for b in cal['bands']:
  fraction=b['q']/sum(cal['q']);lo,hi=beta.ppf([ALPHA/2,1-ALPHA/2],b['df']/2,(15968-b['df'])/2);bands.append(bool(lo<=fraction<=hi))
 cal['checks']=dict(rho_valid=valid,bands=bands,channel_p=pc,mean_p=pm,temporal_variance_p=pv,lag2_p=pl,temporal_q=qtime,temporal_df=df,lag2=lag);cal['adequate']=bool(valid and all(bands) and min(pc,pm,pv,pl)>=ALPHA);return cal

def make_tracker(model,L,seed,cal,keep=8):
 obj=PersistentTracker(model,L*np.sqrt(cal['scale']),cal['rho'],seed,secondary_keep=keep);obj.qcut=float(998*f.ppf(.999,998,15968));return obj
