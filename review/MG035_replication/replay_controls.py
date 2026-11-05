"""Replay archived late-control fields with an independent sparse lattice implementation."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp

def replay(field, b):
    n = field.shape[1]
    eye = sparse.eye(n, format='csr')
    lap = sparse.diags([np.ones(n-1), -2*np.ones(n), np.ones(n-1)], [-1,0,1], format='lil')
    lap[0,-1] = lap[-1,0] = 1
    lap = lap.tocsr() / .2**2
    x0 = 140/9
    def rhs(t, u):
        g,x,y = u.reshape(3,n)
        q = b*x + x0*x0*y + b/x0*x*x + 2*x0*x*y + x*x*y
        return np.concatenate([lap@g-g+.1*x, lap@x+g-x+q, 12*(lap@y)-q])
    def jac(t, u):
        _,x,y = u.reshape(3,n)
        a = sparse.diags(b + 2*b/x0*x + 2*(x0+x)*y)
        c = sparse.diags((x0+x)**2)
        return sparse.bmat([[lap-eye,.1*eye,None],[eye,lap-eye+a,c],[None,-a,12*lap-c]],format='csc')
    times = np.r_[200:301:10,1900:2001:10]
    sol = solve_ivp(rhs,(0,2000),field.astype(float).ravel(),method='BDF',jac=jac,
                    atol=1e-10,rtol=1e-8,t_eval=times)
    assert sol.success and np.array_equal(sol.t,times)
    return sol.y.T.reshape(-1,3,n)

def classify(seq, depth, wavelength):
    paths = None
    for f in seq[:,1,:]:
        peaks=np.flatnonzero((f < np.roll(f,1)) & (f <= np.roll(f,-1)) & (f < -abs(depth)/2))
        if paths is None:
            paths=set(peaks)
        else:
            updated=set()
            for prev in paths:
                if len(peaks):
                    dist=abs((peaks-prev+len(f)//2)%len(f)-len(f)//2)*.2
                    if min(dist)<=wavelength: updated.add(peaks[np.argmin(dist)])
            paths=updated
    return 'persistent' if paths else ('excursion' if abs(seq[-1]).max()<.001 else 'unresolved')

if __name__=='__main__':
    p=Path(sys.argv[1])/'results/MG9b'
    out=[]
    for tag in ['29.20','29.35']:
        ref=json.loads((p/f'vb1_{tag}.json').read_text())
        records=json.loads((p/f'vb2_{tag}.json').read_text())['pairs']
        with np.load(p/f'vb2_{tag}_fields.npz') as fields:
            assert len(fields.files)==len(records)==10
            for r in records:
                key=f"{r['cell']}__{r['rep']}__{r['sc']}"
                seq=replay(fields[key],float(tag))
                labels=[classify(s,ref['A'],ref['lam']) for s in (seq[:11],seq[11:])]
                assert labels==[r['std_cls'],r['late_cls']], (key,labels,r)
                out.append(dict(b=tag,key=key,standard=labels[0],late=labels[1]))
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(len(out),'archived fields matched at both windows')
