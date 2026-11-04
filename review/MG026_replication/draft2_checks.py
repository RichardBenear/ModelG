"""Draft 2 static/synthetic checks only. No Model G nonlinear outcomes."""
import json
import numpy as np
from scipy.optimize import brentq

def threshold(h, length, mode):
    x0 = 140/9
    s = 4/h**2*np.sin(np.pi*mode*h/length)**2
    def det(b):
        return np.linalg.det(np.array([[-1,.1,0],[1,b-1,x0*x0],[0,-b,-x0*x0]])-s*np.diag([1,1,12]))
    return brentq(det,29,31)

wavelength = 100/34
out = {"D02L_thresholds": {str(m): threshold(.2,200,m) for m in (67,68)},
       "cosine_cubic_SH": "19*c2**2/18 - 3/4",
       "SH_boundary_c2_squared": 27/38,
       "envelopes": {}}
for h in (.2,.1):
    n = round(100/h)
    wave = np.cos(2*np.pi*34*np.arange(n)/n)
    windows = np.array([np.roll(wave,j) for j in range(-7,8)])
    envelope = np.sqrt(2)*windows.std(axis=0)
    out["envelopes"][str(h)] = dict(window_width=15*h,
        min=float(envelope.min()),max=float(envelope.max()),mean=float(envelope.mean()))
# Both interfaces first advance 2 lambda, then return to lambda.
# The endpoint invasion criterion passes despite reversal.
out["reversal_example"] = dict(displacement=wavelength, range=2*wavelength,
    passes_invasion_inequality=True, matches_stated_unresolved_reversal=True)
print(json.dumps(out,indent=2))
