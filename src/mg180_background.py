"""Joint-profile reference variance test integrated with alignment tracker."""
import numpy as np
from scipy.stats import f
from mg178_background import calibrate as old_calibrate,ALPHA
from mg179_profile import fit_reference
from mg174_modes import ModeTracker
from mg180_alignment import AlignedTracker

def calibrate(model,L,spatial,temporal,seed):
 cal=old_calibrate(model,L,spatial,temporal,seed);base=ModeTracker(model,L,0,seed);yt=np.fft.rfft(temporal-temporal.mean(axis=2)[:,:,None],axis=2);stats=[sum(cal['q']),float(base.inner(yt[1:],yt[1:]).sum()),float(base.inner(yt[1:],yt[:-1]).sum()),float(base.inner(yt[:-1],yt[:-1]).sum()),float(base.inner(yt[0],yt[0]))];profile=fit_reference(*stats);cal['profile']=profile;cal['old_temporal_variance_p']=cal['checks']['temporal_variance_p'];cal['checks']['temporal_variance_p']=profile['p'];cal['checks']['rho_valid']=cal['checks']['rho_valid'] and profile['valid'];cal['adequate']=bool(cal['checks']['rho_valid'] and all(cal['checks']['bands']) and min(cal['checks'][k] for k in ['channel_p','mean_p','temporal_variance_p','lag2_p'])>=ALPHA);cal['conditional_rho']=cal['rho'];cal['conditional_scale']=cal['scale'];return cal

def make_tracker(model,L,seed,cal,keep=8):
 obj=AlignedTracker(model,L*np.sqrt(cal['scale']),cal['rho'],seed,secondary_keep=keep);obj.qcut=float(998*f.ppf(.999,998,15968));return obj
