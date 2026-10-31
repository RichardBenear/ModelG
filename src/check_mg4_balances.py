"""MG4 balance-check investigation (GATES.md MG4, Part A result): recomputes the node-sum balance residuals of saved
final states in float64 (the registered check) and in extended precision (np.longdouble, 80-bit on x86 Linux).
    python src/check_mg4_balances.py            # every run whose registered balance check failed
    python src/check_mg4_balances.py R0_p4_sol  # named runs
"""
import sys,json,math,numpy as np
sys.path.insert(0,'src'); import run_mg4_partA as pa, graph_ode as go, graphs as gr
from lattice import steady_state,P17
dec=json.load(open(pa.DEC))
def rhs_ld(g,h,y):
    P=P17; G0,X0,Y0=[np.longdouble(v) for v in steady_state(P)]; b,gg,p,q=[np.longdouble(P[k]) for k in ("b","g","p","q")]
    N=g["N"]; y=y.astype(np.longdouble); I,J=g["I"],g["J"]; hh=np.longdouble(h)**2
    def L(f):
        out=np.zeros(N,np.longdouble); d=(f[J]-f[I]); np.add.at(out,I,d); np.add.at(out,J,-d); return out/hh
    pG,pX,pY=y[:N],y[N:2*N],y[2*N:]; X,Y=X0+pX,Y0+pY; r4=X*X*Y-X0*X0*Y0
    return np.concatenate([L(pG)-q*pG+gg*pX, np.longdouble(P["dx"])*L(pX)+p*pG-(1+b)*pX+r4, np.longdouble(P["dy"])*L(pY)+b*pX-r4]), y
import glob, os
keys = sys.argv[1:] or sorted(os.path.basename(f)[:-5] for f in glob.glob("results/MG4/partA/runs/*.json")
                              if not json.load(open(f))["balances"]["ok"])
print("run | at run time (float64, registered check): r2, tolerance | saved float32 state recomputed: float64 r1, r2 | "
      "extended precision r1, r2 | sum|rhs|")
for key in keys:
    r=json.load(open(f"results/MG4/partA/runs/{key}.json")); g,h=pa.get_graph(dict(fam=r["fam"],slot=r["slot"],torus=r["torus"]),dec); N=g["N"]
    y=np.load(pa.field_file(key))["final"].astype(float).ravel()
    rhs,_,_=go.make_system(g,h); d=rhs(300,y); b64=go.balances(y,d,N)
    dl,yl=rhs_ld(g,h,y); I=[yl[k*N:(k+1)*N].sum() for k in range(3)]; dI=[dl[k*N:(k+1)*N].sum() for k in range(3)]
    r2=sum(dI)-(-(1-1)*I[0]-(1-0.1)*I[1]); r1=dI[0]-(-I[0]+np.longdouble(0.1)*I[1])
    scale=float(np.abs(d).sum())
    rb = r["balances"]
    print(key, "| run r2 %.2e tol %.2e |" % (rb["r2"], 1e-10 + 1e-9 * rb["M2"]), "f64 r1 %.2e r2 %.2e | longdouble r1 %.2e r2 %.2e | sum|rhs| %.2e (eps*that %.1e)"%(b64["r1"],b64["r2"],float(r1),float(r2),scale,scale*2.2e-16))
