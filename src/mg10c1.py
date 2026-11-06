"""MG044 follow-up numerical model and prospective resolved-minimum diagnostic."""
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp, BDF
from scipy.signal import find_peaks
import mg10_reservoir as old
from measure import smooth_periodic, circ

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/MG10C1"
SCALE=np.array([140/9,140/9,29/(140/9),1.,1.])
B0=29*(140/9)*300


def background_rhs(t,z):
    g,x,y,ua,ub=z
    source=14*(1-ua); gx=g; xg=.1*x
    xy=29*x*(1-ub); yx=x*x*y; sink=.9*x
    return np.array([source-gx+xg,gx-xg-xy+yx-sink,xy-yx,source/4200,xy/B0])


def background_ref_rhs(t,z):
    g,x,y,u,v=z
    return np.array([-g+.1*x+14-14*u,
                     g-(30-29*v)*x+x*x*y,
                     (29-29*v)*x-x*x*y,
                     (1-u)/300,
                     x*(1-v)/((140/9)*300)])


def reduced(job):
    primary=job["kind"]=="reduced_primary"
    rhs=background_rhs if primary else background_ref_rhs
    horizon=job.get("T",300)
    sol=solve_ivp(rhs,(0,horizon),[140/9,140/9,29/(140/9),0,0],
                  method="DOP853" if primary else "Radau",rtol=1e-11,
                  atol=[1e-12]*3+[1e-14]*2,max_step=.1,dense_output=True)
    if not sol.success: raise RuntimeError(sol.message)
    # Serialize the actual adaptive dense polynomials, not a coarser resampling.
    # SciPy 1.16.1 DOP853 F and Radau Q layouts are verified by qualification.
    knots=sol.sol.ts
    coefficients=[]
    for segment in sol.sol.interpolants:
        if primary:
            c=np.zeros((1,5))
            for i,f in enumerate(reversed(segment.F)):
                c[0]+=f
                new=np.zeros((len(c)+1,5))
                if i%2==0: new[1:]=c
                else: new[:-1]+=c;new[1:]-=c
                c=new
            c[0]+=segment.y_old
        else:
            c=np.vstack([segment.y_old,segment.Q.T])
        padded=np.zeros((8,5));padded[:len(c)]=c
        coefficients.append(padded)
    coefficients=np.array(coefficients)
    spline=polynomial_solution(knots,coefficients)
    midpoint=(knots[1:]+knots[:-1])/2
    interpolation_error=float(np.max(abs(np.array([spline(t) for t in midpoint])-sol.sol(midpoint).T)/SCALE))
    if interpolation_error>1e-10: raise RuntimeError(f"Background interpolation {interpolation_error}")
    times=np.arange(horizon+1,dtype=float)
    z=sol.sol(times).T
    if not np.all(np.isfinite(z)) or np.any(z<0) or np.any(z[:,3:]>1): raise FloatingPointError("Reduced state")
    return dict(times=times,state=z,knots=knots,coefficients=coefficients,
                interpolation_error=np.array(interpolation_error))


def polynomial_solution(knots,coefficients):
    def value(t):
        i=int(np.clip(np.searchsorted(knots,t,side="left")-1,0,len(knots)-2))
        u=(t-knots[i])/(knots[i+1]-knots[i])
        result=coefficients[i,-1].copy()
        for row in coefficients[i,-2::-1]: result=result*u+row
        return result
    return value


def base_spline(file):
    with np.load(file) as raw:
        return polynomial_solution(raw["knots"],raw["coefficients"])


def mode_matrix(background,modes):
    g,x,y,ua,ub=background
    modes=np.asarray(modes)
    lap=4*np.sin(np.pi*modes/500)**2/.2**2
    mat=np.zeros((len(modes),4,4))
    mat[:,0,0]=-1-lap; mat[:,0,1]=.1
    mat[:,1,0]=1; mat[:,1,1]=-1-29*(1-ub)+2*x*y-lap
    mat[:,1,2]=x*x; mat[:,1,3]=29*x
    mat[:,2,1]=29*(1-ub)-2*x*y; mat[:,2,2]=-x*x-12*lap
    mat[:,2,3]=-29*x
    mat[:,3,1]=29*(1-ub)/B0; mat[:,3,3]=-29*x/B0
    scales=SCALE[[0,1,2,4]]
    return mat*scales[None,None,:]/scales[None,:,None]


def linear(job,directory):
    primary=job["kind"]=="linear_primary"
    background=base_spline(Path(directory)/("reduced_primary.npz" if primary else "reduced_reference.npz"))
    modes=np.array(job.get("modes",list(range(1,251))),dtype=int)
    n=len(modes)
    def rhs(t,flat):
        return np.einsum("mij,jm->im",mode_matrix(background(t),modes),flat.reshape(4,n)).ravel()
    def jac(t,flat):
        matrices=mode_matrix(background(t),modes)
        return sp.bmat([[sp.diags(matrices[:,r,c]) for c in range(4)] for r in range(4)],format="csc")
    initial=np.zeros((4,n));initial[1]=1
    horizon=job.get("T",300)
    sol=solve_ivp(rhs,(0,horizon),initial.ravel(),method="BDF" if primary else "Radau",
                  jac=jac,rtol=1e-9,atol=1e-12,max_step=.2,t_eval=np.arange(horizon+1,dtype=float))
    if not sol.success or not np.all(np.isfinite(sol.y)): raise RuntimeError("Linear integration incomplete/nonfinite")
    return dict(times=sol.t,response=sol.y.T.reshape(-1,4,n),modes=modes)


def full(job,checkpoint=None):
    horizon=job.get("T",300)
    n=500
    state=old.initial("H")
    state[1]+=SCALE[1]*job["epsilon"]*np.cos(2*np.pi*job["mode"]*np.arange(n)/n)
    saved=[state.copy()];times=[0.]
    if job["method"]=="E2":
        for t in range(1,horizon+1):
            old.advance_euler(state,1,True,6000,6000)
            if not old.valid_state(state): raise FloatingPointError(f"Euler t={t}")
            saved.append(state.copy());times.append(float(t))
            if checkpoint and (t%25==0 or t==horizon): checkpoint(dict(times=np.array(times),state=np.array(saved)))
    elif job["method"]=="Oplus":
        indices=np.concatenate([np.arange(i*n,(i+1)*n) for i in range(15) if i!=11])
        def unpack(v):
            a=np.zeros(15*n);a[indices]=v;return a
        atol=np.repeat([1e-12]*3+[1e-13]*2+[1e-9]*10,n)[indices]
        solver=BDF(lambda t,z:old.ode_rhs(t,unpack(z),1,True)[indices],0.,state.ravel()[indices],float(horizon),
                   rtol=1e-10,atol=atol,max_step=.25,
                   jac=lambda t,z:old.ode_jac(t,unpack(z),1,True)[indices,:][:,indices])
        next_time=1
        while solver.status=="running":
            solver.step()
            if solver.status=="failed": raise RuntimeError("BDF failed")
            if next_time<=solver.t:
                dense=solver.dense_output()
                while next_time<=solver.t and next_time<=horizon:
                    snapshot=unpack(dense(next_time)).reshape(15,n)
                    if not old.valid_state(snapshot): raise FloatingPointError(f"ODE t={next_time}")
                    saved.append(snapshot.copy());times.append(float(next_time))
                    if checkpoint and (next_time%25==0 or next_time==horizon): checkpoint(dict(times=np.array(times),state=np.array(saved)))
                    next_time+=1
    else: raise ValueError(job["method"])
    return dict(times=np.array(times),state=np.array(saved))


def resolved_minima(x_field,ref):
    """Prominence-qualified old depth rule; counts resolved minima, not particles."""
    phi=np.asarray(x_field,dtype=float)-SCALE[1]
    mean=phi.mean()
    smooth=smooth_periodic(phi-mean,.2)+mean
    eta=1e-8*max(1.,SCALE[1],abs(ref["A_s"]))
    n=len(phi)
    peaks,_=find_peaks(-np.tile(smooth,3),prominence=eta)
    indices=peaks[(peaks>=n)&(peaks<2*n)]-n
    indices=indices[smooth[indices]<ref["thr"]]
    kept=[]
    for i in indices[np.argsort(smooth[indices])]:
        if all(circ(ref["x"][i],ref["x"][j],100)>=ref["lam"]/2 for j in kept): kept.append(i)
    return np.array(sorted(kept),dtype=int)


def fixtures(ref):
    x=ref["x"];rows=[]
    eta=1e-8*max(1.,SCALE[1],abs(ref["A_s"]))
    for level in (0.,.1,1.,5.,10.,SCALE[1],30.):
        rows.append((f"flat_{level}",np.full(500,level),0))
        rows.append((f"roundoff_{level}",np.full(500,level)+1e-12*np.cos(2*np.pi*16*np.arange(500)/500),0))
    rows.append(("reference",SCALE[1]+ref["phi"][1],1))
    for count in (1,2):
        f=np.full(500,5.)
        for center in (0.,25.)[:count]:
            d=circ(x,center,100)
            f-=np.exp(-.5*(d/.8)**2)
        for shift in (0,137,250): rows.append((f"gaussian_{count}_shift{shift}",np.roll(f,shift),count))
    phase=2*np.pi*16*np.arange(500)/500
    rows.append(("subfloor_ripple",5+eta/100*np.cos(phase),0))
    rows.append(("resolved_ripple",5+eta*100*np.cos(phase),16))
    return rows
