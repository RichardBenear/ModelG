"""Static MG9 design checks only. No nonlinear integration or P1/P2 outcome computed."""
import json
import numpy as np
from scipy.optimize import brentq, minimize_scalar

D = np.diag([1., 1., 12.])
X0 = 140/9

def jac(b):
    return np.array([[-1., .1, 0.], [1., b-1, X0**2], [0., -b, -X0**2]])

def threshold(s):
    return brentq(lambda b: np.linalg.det(jac(b)-s*D), 0, 10000)

opt = minimize_scalar(threshold, bounds=(.1,20), method="bounded",
                      options={"xatol":1e-13})
s = opt.x
out = dict(continuous_k=dict(bc=threshold(s), effective_k_squared=s,
                             continuum_k=float(np.sqrt(s))))
for h in (.2,.1):
    n = round(100/h)
    modes = np.arange(1,n//2+1)
    lam = 4/h**2*np.sin(np.pi*modes/n)**2
    bs = np.array([threshold(z) for z in lam])
    m = int(np.argmin(bs))
    bc = bs[m]
    values, vectors = np.linalg.eig(jac(bc)-lam[m]*D)
    i = np.argmin(abs(values))
    v = vectors[:,i].real
    v /= v[1]
    vl, lefts = np.linalg.eig((jac(bc)-lam[m]*D).T)
    left = lefts[:,np.argmin(abs(vl))].real
    left /= np.dot(left,v)
    jp = np.array([[0,0,0],[0,1,0],[0,-1,0]])
    out[str(h)] = dict(continuous_k=float(2/h*np.arcsin(h*np.sqrt(s)/2)),
                       finite_L100_bc=bc, critical_mode=int(modes[m]),
                       finite_k=float(2*np.pi*modes[m]/100), alpha=float(left@jp@v))
out["listed_integrations_minimum"] = 80+33+5*14
out["declared_cap"] = 120
print(json.dumps(out,indent=2))
