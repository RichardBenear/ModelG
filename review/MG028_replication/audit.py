"""Independent MG028 algebra and records audit; no production-module imports or ODE reruns."""
import json
import hashlib
import sys
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.optimize import brentq, minimize_scalar

root=Path(sys.argv[1])/"results/MG9"
read=lambda p:json.loads((root/p).read_text())
p1=read("p1.json")
ev=read("evaluation.json")
x0=140/9
D=np.diag([1,1,12])
def J(b):
    return np.array([[-1,.1,0],[1,b-1,x0*x0],[0,-b,-x0*x0]])
def bc(s):
    return brentq(lambda b:np.linalg.det(J(b)-s*D),29,40,xtol=1e-13)
sopt=minimize_scalar(bc,bounds=(3,6),method="bounded",options={"xatol":1e-12}).x
phase=np.linspace(0,2*np.pi,128,endpoint=False)
cos=np.cos(phase)
coeff={}
for name,h,s in [("i",.2,4/.2**2*np.sin(np.pi*34*.2/100)**2),
                 ("ii",.2,sopt),("iii",None,sopt),
                 ("iv",.1,4/.1**2*np.sin(np.pi*34*.1/100)**2)]:
    b=bc(s)
    matrix=J(b)-s*D
    _,_,vh=np.linalg.svd(matrix)
    v=vh[-1];v=v/v[1]
    _,_,vh=np.linalg.svd(matrix.T)
    left=vh[-1];left=left/(left@v)
    # Obtain mean and second harmonic forcing by phase quadrature.
    leading=v[:,None]*cos
    quadratic=(b/x0)*leading[1]**2+2*x0*leading[1]*leading[2]
    forcing=np.array([0,1,-1])[:,None]*quadratic
    w0=np.linalg.solve(J(b),-forcing.mean(axis=1))
    s2=4*s if h is None else 4*s*(1-h*h*s/4)
    w2=np.linalg.solve(J(b)-s2*D,-2*(forcing*np.cos(2*phase)).mean(axis=1))
    correction=w0[:,None]+w2[:,None]*np.cos(2*phase)
    cubic=(2*(b/x0)*leading[1]*correction[1]+
           2*x0*(leading[1]*correction[2]+leading[2]*correction[1])+
           leading[1]**2*leading[2])
    coefficient=float((left@np.array([0,1,-1]))*2*np.mean(cubic*cos))
    np.testing.assert_allclose(coefficient,p1["L3"][name]["L3"],rtol=2e-6)
    coeff[name]=coefficient

fits={}
for path in sorted((root/"runs").glob("V2_*.json")):
    r=json.loads(path.read_text())
    assert r["complete"] and r["stop_t"]>=40
    t,a=np.array(r["t"]),np.array(r["A"])
    win=t>=20
    tau=t[win]-20
    sig=p1["sigma34_at_b_i"]
    basis=np.expm1(-2*sig*tau)/sig
    design=np.column_stack([np.exp(-2*sig*tau),basis])
    intercept,l3=np.linalg.lstsq(design,1/a[win]**2,rcond=None)[0]
    change=abs(a[win][0]**2/a[win][-1]**2-1)
    fits[path.stem]=dict(l3=float(l3),inverse_square_change=float(change))
    assert abs(l3-coeff["i"])/coeff["i"]<.001
assert fits["V2_A0.05"]["inverse_square_change"]<.25
assert all(fits[k]["inverse_square_change"]>=.25 for k in ["V2_A0.1","V2_A0.2"])
assert abs(fits["V2_A0.1"]["l3"]-fits["V2_A0.1_tight"]["l3"])/fits["V2_A0.1"]["l3"]<.02

lam=100/34
counts=Counter()
slow={}
for path in sorted((root/"runs").glob("front_grid_*.json")):
    r=json.loads(path.read_text())
    if "track" not in r:
        counts[r["label"]]+=1
        continue
    assert r["complete"] and r["success"] and r["t_reached"]==r["T"]
    tr=r["track"]
    status=[v["status"] for v in tr]
    first=next((i for i,s in enumerate(status) if s in ("fill","loss")),None)
    if first is not None:
        assert len(set(status[first:]))==1
        label="INVADES" if status[first]=="fill" else "RETREATS"
    else:
        assert set(status)=={"ok"}
        tr=[v for v in tr if v["t"]>=100]
        assert not any(v["breakup"] or v["remote"] for v in tr)
        t=np.array([v["t"] for v in tr])
        disps=[-np.array([v["xL"] for v in tr]),np.array([v["xR"] for v in tr])]
        disps=[v-v[0] for v in disps]
        end=[float(v[-1]) for v in disps]
        vel=[float(np.polyfit(t[t>=t[-1]-200],v[t>=t[-1]-200],1)[0]) for v in disps]
        # None of these recorded tracks reverses by half a wavelength.
        for v in disps:
            assert np.max(np.maximum.accumulate(v)-v)<lam/2
        if min(end)>=lam:
            label="INVADES"
        else:
            assert r["T"]==1800 and all(lam/2<v<lam for v in end)
            assert r["residual"]>1e-5
            label="UNRESOLVED-SLOW"
            slow[str(r["b"])]=dict(displacements=end,late_velocities=vel,residual=r["residual"])
    assert label==r["label"]
    counts[label]+=1
assert sum(counts.values())==31 and counts["UNRESOLVED-SLOW"]==3 and counts["INVADES"]==13
assert ev["P1"]=="SUBCRITICAL" and ev["P2"]=="UNRESOLVED"

# Independently reconstruct stationary RHS from every saved preparation state.
z=np.load(root/"states.npz")
residual_error=0.
hashes=0
for key in z.files:
    path=root/"runs"/(key+".json")
    if not path.exists():
        continue
    r=json.loads(path.read_text())
    f=z[key]
    if "sha256" in r:
        assert hashlib.sha256(np.asarray(f,dtype=np.float64).tobytes()).hexdigest()==r["sha256"]
        hashes+=1
    if "b" not in r or key.startswith("bench"):
        continue
    g,x,y=f.reshape(3,-1)
    lap=lambda v:(np.roll(v,1)+np.roll(v,-1)-2*v)/.2**2
    nonlinear=(x0+x)**2*(r["b"]/x0+y)-x0*r["b"]
    rhs=np.array([lap(g)-g+.1*x,lap(x)+g-(1+r["b"])*x+nonlinear,
                  12*lap(y)+r["b"]*x-nonlinear])
    residual_error=max(residual_error,abs(float(abs(rhs).max())-r["residual"]))
assert residual_error<1e-9
ledger=read("ledger.json")["integrations"]
assert len(ledger)==150 and len({r["key"] for r in ledger})==150
print(json.dumps(dict(cubic_coefficients=coeff,growth_fits=fits,grid_counts=dict(counts),slow_fronts=slow,
                      state_hashes_checked=hashes,maximum_residual_recompute_error=residual_error,
                      integrations=len(ledger),core_hours=sum(r["wall_s"] for r in ledger)/3600),indent=2))
