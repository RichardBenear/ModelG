"""MG192 post-hoc diagnostic (not a registered test): continuation of the uniform soliton in a, steps 2e-5 relative, until Newton fails."""
import sys, json; sys.path.insert(0,'src')
import numpy as np, mg192_mobility as m
d=np.load('results/MG192_gravity_mobility/soliton_adjoint_L100_N1024.npz');U0=d['U']
g=m.Grid(100.,1024);D2=g.D2();E,rows=m.even_maps(g.N)
def solve(u,a,b=29.0):
    for it in range(25):
        F=m.rhs(g,u,a=a,b=b).ravel();r=np.max(np.abs(F))
        if r<1e-10: return u,r,True
        J=m.jac(g,u,D2,b=b);Jr=J[rows]@E
        u=u+(E@np.linalg.solve(Jr,-F[rows])).reshape(3,g.N)
        if np.any(u<=0) or not np.all(np.isfinite(u)): return u,np.inf,False
    return u,r,False
out={}
for direction in (+1,-1):
    u=U0.copy();rows_out=[]
    for i in range(1,200):
        dl=direction*i*2e-5
        a=14*(1+dl);u2,r,ok=solve(u.copy(),a)
        depth=float(u2[1].min()-m.hom(a=a)[1]) if ok else None
        if ok:
            J=m.jac(g,u2,D2);ev=np.linalg.eigvals(J[rows]@E);mx=float(np.max(ev.real))
        else: mx=None
        rows_out.append((dl,ok,depth,mx))
        print(direction,dl,ok,depth,mx,flush=True)
        if not ok or depth>-1.0: break
        u=u2
    out[str(direction)]=rows_out
json.dump(out,open('results/MG192_gravity_mobility/posthoc_fold.json','w'))
