"""MG195 (confirmatory): inverse-square attraction in the conserved-G test model (EX008), 3D.

Single-particle M and mu (radial, l = 1 linear response) at P1-P3 and C2; direct axisymmetric co-moving pair solves
(EX009 solver) at d = 18, 22, 27, 33, 40; controls C1 (k = 0) and C2 (ether G exchange lambda (G_inf - G),
lambda = 0.01); numerical checks at P1, d = 22. Spec: spec/MG195_inverse_square.md.
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import ex009_pair2d as p2  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "MG195_inverse_square")
PARTIAL = os.path.join(OUT, "partial")
POINTS = {"P1": dict(cx=1, cy=0, k=0.2, lam=0.0), "P2": dict(cx=0, cy=1, k=0.2, lam=0.0),
          "P3": dict(cx=1, cy=0, k=0.05, lam=0.0), "C1": dict(cx=1, cy=0, k=0.0, lam=0.0),
          "C2": dict(cx=1, cy=0, k=0.2, lam=0.01)}
DS = (18.0, 22.0, 27.0, 33.0, 40.0)


class RadialL(e.Model):
    def __init__(self, R, dr, k, cx, cy, lam):
        super().__init__(R, dr, k, cx, cy)
        self.lam = lam

    def F(self, u):
        f = super().F(u)
        f[3] += self.lam * (e.GINF - u[3])
        return f

    def J(self, u, l=0):
        n = self.n
        return (super().J(u, l) - sp.block_diag([sp.csr_matrix((3 * n, 3 * n)), self.lam * sp.identity(n)])).tocsc()


class AxiL(p2.Axi):
    lam = 0.0

    def F(self, u):
        f = super().F(u)
        f[3] += self.lam * (e.GINF - u[3])
        return f

    def J(self, u):
        n = self.n
        return (super().J(u) - sp.block_diag([sp.csr_matrix((3 * n, 3 * n)), self.lam * sp.identity(n)])).tocsc()


def radial(name, R):
    c = POINTS[name]
    prev = np.load(os.path.join(e.OUT, "particle_GX_k0.05_R60.npz"))
    M = RadialL(R, 0.05, 0.0, c["cx"], c["cy"], c["lam"])
    u = np.array([np.interp(M.r, prev["r"], prev["u"][s]) for s in range(4)])
    u[3] = e.GINF
    u, res, ok = M.newton(u)
    for kk in np.linspace(0, c["k"], 6)[1:]:
        M = RadialL(R, 0.05, kk, c["cx"], c["cy"], c["lam"])
        u, res, ok = M.newton(u)
    out = dict(point=name, R=R, ok=bool(ok and u[2, 0] > 1), residual=res, M=M.integral(c["k"] * u[3] * u[2]))
    if c["k"] > 0:
        out["mu"] = M.drift(u)["mu"]
    out["eigs_l0"] = e.leading_eigs(M, u, 0)
    out["eigs_l1"] = e.leading_eigs(M, u, 1)
    np.savez_compressed(os.path.join(OUT, f"radial_{name}_R{R:g}.npz"), r=M.r, u=u)
    return out


def predict(name, d, rad):
    c = POINTS[name]
    if c["k"] == 0:
        return 0.0
    kap = np.sqrt(c["lam"])
    return rad["mu"] * rad["M"] / (4 * np.pi) * np.exp(-kap * d) * (1 / d ** 2 + kap / d)


def pair(args):
    name, d, h, far = args
    c = POINTS[name]
    p2.K, p2.CX, p2.CY = c["k"], c["cx"], c["cy"]
    AxiL.lam = c["lam"]
    prof = np.load(os.path.join(OUT, f"radial_{name}_R120.npz"))
    Mrad = float(np.sum(0)) if c["k"] == 0 else json.load(open(os.path.join(OUT, "radial.json")))[name]["M"]
    a = AxiL(d, h=h, far=far)
    t0 = time.time()
    u, v, res, ok, it = a.solve(a.guess(prof["r"], prof["u"], Mrad))
    tag = f"{name}_d{d:g}_h{h:g}_far{far:g}"
    os.makedirs(PARTIAL, exist_ok=True)
    np.savez_compressed(os.path.join(PARTIAL, f"pair_{tag}.npz"), u=u, v=v, d=d, h=h, far=far)
    return dict(point=name, d=d, h=h, far=far, ok=bool(ok), residual=res, iters=it, v=float(v),
                M_particle=a.integral(c["k"] * u[3] * u[2]), wall=time.time() - t0)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rad = {}
    rads = []
    for name in POINTS:
        for R in (60.0, 120.0):
            r_ = radial(name, R)
            print("radial", r_, flush=True)
            rads.append(r_)
        rad[name] = r_                       # R = 120 values are used for predictions
    json.dump(rad, open(os.path.join(OUT, "radial.json"), "w"), indent=1, default=float)
    jobs = [(n, d, 0.1, 150.0) for n in POINTS for d in DS] + [("P1", 22.0, 0.075, 150.0), ("P1", 22.0, 0.1, 300.0)]
    res = []
    with Pool(int(os.environ.get("MODELG_WORKERS", "4"))) as pool:
        for r_ in pool.imap_unordered(pair, jobs):
            r_["prediction"] = predict(r_["point"], r_["d"], rad[r_["point"]])
            print("pair", r_, flush=True)
            res.append(r_)
            json.dump(dict(radial_all=rads, radial=rad, pairs=res), open(os.path.join(OUT, "results.json"), "w"),
                      indent=1, default=float)
    out = dict(radial_all=rads, radial=rad, pairs=res, wall_seconds=time.time() - t0)
    out["verdicts"] = score(out)
    json.dump(out, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print(json.dumps(out["verdicts"], indent=1, default=float))


def score(out):
    P = {(r["point"], r["d"], r["h"], r["far"]): r for r in out["pairs"]}
    base = {k: v for k, v in P.items() if k[2] == 0.1 and k[3] == 150.0}
    c1 = max(abs(P[("C1", d, 0.1, 150.0)]["v"]) for d in DS)
    v22 = P[("P1", 22.0, 0.1, 150.0)]["v"]
    g = abs(P[("P1", 22.0, 0.075, 150.0)]["v"] - v22) / abs(v22)       # relative numerical changes
    f = abs(P[("P1", 22.0, 0.1, 300.0)]["v"] - v22) / abs(v22)
    floor = 10 * c1                                                     # absolute null floor
    num_ok = max(g, f) <= 0.02
    ver = dict(floor=floor, c1_max=c1, grid_change_rel=g, far_change_rel=f, numerics_ok=num_ok, points={})
    for n in ("P1", "P2", "P3"):
        rows = [base[(n, d, 0.1, 150.0)] for d in DS]
        conv = all(r["ok"] and r["residual"] <= 1e-8 for r in rows) and num_ok
        neg = all(r["v"] < 0 and abs(r["v"]) > floor for r in rows)
        ratios = {r["d"]: r["v"] / r["prediction"] for r in rows}
        within = all(abs(ratios[d] - 1) <= 0.05 for d in DS if d >= 22)
        dd = np.array([d for d in DS if d >= 22])
        vv = np.array([abs(base[(n, d, 0.1, 150.0)]["v"]) for d in dd])
        expo = float(-np.polyfit(np.log(dd), np.log(vv), 1)[0])
        if min(abs(r["prediction"]) for r in rows) <= floor:
            conv = False
        if not conv:
            verdict = "UNRESOLVED"
        elif neg and within and abs(expo - 2) <= 0.10:
            verdict = "SUPPORTED"
        else:
            verdict = "NOT SUPPORTED"
        ver["points"][n] = dict(converged=conv, all_attractive_above_floor=neg, ratios=ratios, exponent=expo,
                                verdict=verdict)
    c2 = base[("C2", 40.0, 0.1, 150.0)]
    inv = P[("P1", 40.0, 0.1, 150.0)]["prediction"]
    ver["C2"] = dict(v40=c2["v"], yukawa_pred=c2["prediction"], inverse_square_pred_P1=inv,
                     within_factor2=bool(0.5 <= c2["v"] / c2["prediction"] <= 2.0) if c2["prediction"] else False,
                     suppressed=bool(abs(c2["v"]) < 0.5 * abs(inv)))
    ver["study"] = ("INVERSE-SQUARE CONFIRMED" if all(ver["points"][n]["verdict"] == "SUPPORTED" for n in ("P1", "P2", "P3"))
                    and ver["C2"]["within_factor2"] and ver["C2"]["suppressed"] else "per point")
    return ver


if __name__ == "__main__":
    main()
