"""Audit the fixed-tau Omega->infinity map against the proposed ODE reference."""
import json
from pathlib import Path
import numpy as np
import numba as nb
from lattice_snapshot import P17, lattice_ode_1d

@nb.njit
def advance(c,tau,T,h,b):
    c=c.copy();N=c.shape[0]
    lossG=-np.expm1(-tau);lossX=-np.expm1(-(1+b)*tau)
    for k in range(round(T/tau)):
        v=np.empty_like(c)
        for i in range(N):
            G,X,Y=c[i]
            r4=min(tau*X*X*Y,Y)
            v[i,0]=G+14*tau-lossG*G+.1/(1+b)*lossX*X
            v[i,1]=X+lossG*G-lossX*X+r4
            v[i,2]=Y+b/(1+b)*lossX*X-r4
        for i in range(N):
            for s in range(3):
                D=12 if s==2 else 1
                c[i,s]=v[i,s]+D*tau/h**2*(v[(i-1)%N,s]-2*v[i,s]+v[(i+1)%N,s])
    return c

if __name__=='__main__':
    results=[];h=.2;N=500;back=np.array([140/9,140/9,261/140])
    for b in [29.,29.16]:
        P={**P17,'b':b};back[2]=b/back[1]
        x,G,X,Y,res=lattice_ode_1d(P,N,h,400,seed_scale=1.5)
        phi=np.column_stack([G,X,Y]);c=phi+back
        assert X.min() < -4, 'ODE reference did not nucleate'
        for factor in [1,2]:
            tau=.4*h*h/12/factor
            final=advance(c,tau,100,h,b)-back
            mask=abs(x)<10
            errs=np.linalg.norm((final-phi)[mask],axis=0)/np.linalg.norm(phi[mask],axis=0)
            result=dict(b=b,tau=tau,reference_residual=res,
                        reference_center=phi[N//2].tolist(),map_center=final[N//2].tolist(),
                        E_per_field=errs.tolist(),E=float(max(errs)))
            print(json.dumps(result),flush=True);results.append(result)
    Path(__file__).with_name('limit_profile_audit.json').write_text(json.dumps(results,indent=2))
