"""MG031 Draft 1 static checks only. No new integration or production imports.
Usage: python design_checks.py CHECKOUT
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

root=Path(sys.argv[1])/'results/MG9c'; z=np.load(root/'snapshots.npz')
state=z['front_D02_29.30000_t30000'].reshape(3,-1); n=500; h=.2
k=2*np.pi*np.fft.rfftfreq(n,d=h)
def shift(f,s):
    return np.fft.irfft(np.fft.rfft(f,axis=1)*np.exp(-1j*k*s),n=n,axis=1)
out={'reviewed_revision':'02e3a37ee4328ef9900673f1cec9ccd91010db50',
     'one_node_error':float(abs(shift(state,h)-np.roll(state,1,axis=1)).max()),'round_trip':{}}
for s in [h/4,h/2,3*h/4]:
    err=abs(shift(shift(state,s),-s)-state).max(axis=1)
    out['round_trip'][str(s)]=dict(max_error_per_field_G_X_Y=err.tolist(),max_error=float(err.max()))
tr=json.loads((root/'runs/front_D02_29.30000.json').read_text())['track']
t=np.array([r['t'] for r in tr]); c=np.array([(r['xL']+r['xR'])/2 for r in tr])
mask=t>=12000; u=t[mask]-12000; values=c[mask]-c[mask][0]
# Fit the existing decelerating center as a counterexample to C0, not a forecast gate.
def model(u,p):return p[0]*(-np.expm1(-u/p[1]))
fit=least_squares(lambda p:model(u,p)-values,[-.05,25000],bounds=([-10,100],[0,1e7]),
                  xtol=1e-13,ftol=1e-13,gtol=1e-13)
def slope(a,b):
    ts=np.arange(a,b+1,50); return float(np.polyfit(ts,model(ts-12000,fit.x),1)[0])
old=slope(25000,30000); new=slope(32000,40000)
actualold=float(np.polyfit(t[t>=25000],c[t>=25000],1)[0])
out['restart_counterexample']=dict(description='Smooth continuation of a fit to existing tracks; no restart artifact by construction.',
    fitted_tau=float(fit.x[1]),fit_rms=float(np.sqrt(np.mean((model(u,fit.x)-values)**2))),
    old_window_speed=old,new_window_speed=new,actual_old_midpoint_speed=actualold,
    relative_speed_difference=abs(new-old)/abs(old),exceeds_25_percent_plus_floor=abs(new-old)>.25*abs(old)+5e-9)
out['near_zero_counterexamples']={
    'rho':{'early_speed':0,'late_speed':0,'ratio':'undefined'},
    'fine_sign_reversal':{'baseline':1e-10,'fine':-1e-10,'floor':5e-9,
        'SIMILAR':True,'REDUCED_by_opposite_sign_clause':True},
    'C1_ratio_only':{'baseline_speed':-4e-7,'tight_speed':-4e-7,'baseline_rho':.8,'tight_rho':1.,
        'speed_agrees':True,'ratio_difference':.2,'draft_NUMERICS_SENSITIVE_condition':'not explicit'}}
Path(__file__).with_name('design_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
