"""Independent deterministic audit of the proposed Phase 2 integrator; no gate runs."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_continuous_lyapunov

def audit(tau,h=.2,b=29.,N=500):
    a,g,p,q=14.,.1,1.,1.
    lossG=-np.expm1(-q*tau);lossX=-np.expm1(-(1+b)*tau)
    X=a*tau*(1+b)/((1-g)*lossX)
    G=(a*tau+g/(1+b)*lossX*X)/lossG
    Y=b*lossX/((1+b)*tau*X)
    Jstep=np.array([[1-lossG,g*lossX/(1+b),0],
                    [lossG,1-lossX+2*tau*X*Y,tau*X*X],
                    [0,b*lossX/(1+b)-2*tau*X*Y,1-tau*X*X]])
    lam=-4*np.sin(np.pi*np.arange(N)/N)**2/h**2
    amps=(1+tau*lam[:,None]*np.array([1,1,12]))[:,:,None]*Jstep
    ev=np.linalg.eigvals(amps)
    maxrate=float(np.log(np.abs(ev).max())/tau)
    x0=140/9;y0=261/140
    odeJ=np.array([[-1,.1,0],[1,b-1,x0*x0],[0,-b,-x0*x0]])
    var=np.zeros(3)
    ch=[([1,0,0],a),([-1,1,0],x0),([1,-1,0],g*x0),([0,-1,1],b*x0),([0,-1,0],(1-g)*x0),([0,1,-1],x0*x0*y0)]
    noise=sum(np.outer(v,v)*rate for v,rate in ch)
    for ell in lam:
        drift=odeJ+ell*np.diag([1,1,12])
        B=noise+np.diag(-2*ell*np.array([1,1,12])*np.array([x0,x0,y0]))
        var+=np.diag(solve_continuous_lyapunov(drift,-B))/N
    return dict(tau=tau,limit_fixed_point=dict(G=G,X=X,Y=Y),
                relative_shift_percent=dict(G=100*(G/x0-1),X=100*(X/x0-1),Y=100*(Y/y0-1)),
                max_step_map_growth_rate=maxrate,continuous_lna_variance=var.tolist())

if __name__=='__main__':
    result=[audit(.4*.2**2/12/factor) for factor in [1,2,4,8]]
    Path(__file__).with_name('deterministic_limit_audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
