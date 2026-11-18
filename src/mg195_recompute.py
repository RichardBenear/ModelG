"""MG195 separate recomputation path (spec section 5).

From the saved co-moving pair states (results/MG195_inverse_square/partial/pair_*.npz) and the radial profiles:
  1. v_ls: least-squares speed  v = -<dz u, F(u)> / <dz u, dz u>, with dz u from numpy.gradient on the cell
     centres (not the solver's Dz operator); also the residual |F(u) + v_ls dz u|_max.
  2. M_flux: G consumption of the particle at z = d/2 from the inward G flux through a cylinder box
     rho < 6, |z - d/2| < 6 (uniform h-grid region), with numpy.gradient derivatives (not the reaction integral).
  3. M_radial_flux: the same from the R = 120 radial profile, flux 4 pi r^2 dG/dr at r = 10.
  4. Exponent refit over d = 22..40 from v_ls; prediction ratio with mu from radial.json (shared) and M_radial_flux.
Shared dependencies: saved states, the F() of src/ex009_pair2d.py / src/mg195_inverse_square.py (for item 1's
F), radial.json mu, NumPy. Items 2-4 do not use the reaction integral or the solver's v.
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import mg195_inverse_square as m  # noqa: E402
import ex009_pair2d as p2  # noqa: E402

OUT = m.OUT


def box_flux(a, G, zc, L=6.0):
    nr, nz = a.nr, a.nz
    rc = a.R[:, 0]
    zz = a.Z[0, :]
    Gm = G.reshape(nr, nz)
    dGdr = np.gradient(Gm, rc, axis=0)
    dGdz = np.gradient(Gm, zz, axis=1)
    ir = np.argmin(abs(rc - L))
    jt, jb = np.argmin(abs(zz - (zc + L))), np.argmin(abs(zz - (zc - L)))
    rr, Rs = rc[:ir], rc[ir]
    top = np.trapezoid(2 * np.pi * rr * dGdz[:ir, jt], rr)
    bot = -np.trapezoid(2 * np.pi * rr * dGdz[:ir, jb], rr)
    zs = zz[jb:jt + 1]
    side = 2 * np.pi * Rs * np.trapezoid(dGdr[ir, jb:jt + 1], zs)
    return float(top + bot + side)          # net inward-negative flux = consumption inside the box


def main():
    rad = json.load(open(os.path.join(OUT, "radial.json")))
    res = json.load(open(os.path.join(OUT, "results.json")))
    solver = {(r["point"], r["d"], r["h"], r["far"]): r for r in res["pairs"]}
    Mrf = {}
    for n in m.POINTS:
        pr = np.load(os.path.join(OUT, f"radial_{n}_R120.npz"))
        r, G = pr["r"], pr["u"][3]
        Mrf[n] = float(4 * np.pi * 10.0 ** 2 * np.interp(10.0, r, np.gradient(G, r)))
    rows = []
    for f in sorted(glob.glob(os.path.join(m.PARTIAL, "pair_*.npz"))):
        s = np.load(f)
        name = os.path.basename(f).split("_")[1]
        d, h, far = float(s["d"]), float(s["h"]), float(s["far"])
        c = m.POINTS[name]
        p2.K, p2.CX, p2.CY = c["k"], c["cx"], c["cy"]
        m.AxiL.lam = c["lam"]
        a = m.AxiL(d, h=h, far=far)
        u = s["u"]
        Fu = a.F(u)
        zz = a.Z[0, :]
        dzu = np.array([np.gradient(u[i].reshape(a.nr, a.nz), zz, axis=1).ravel() for i in range(4)])
        w = np.tile(a.V, 4)
        v_ls = float(-np.sum(w * dzu.ravel() * Fu.ravel()) / np.sum(w * dzu.ravel() ** 2))
        resid = float(np.max(np.abs(Fu + v_ls * dzu)))
        Mf = box_flux(a, u[3], d / 2) if c["k"] > 0 else 0.0
        sv = solver[(name, d, h, far)]
        pred = 0.0
        if c["k"] > 0:
            kap = np.sqrt(c["lam"])
            pred = rad[name]["mu"] * Mrf[name] / (4 * np.pi) * np.exp(-kap * d) * (1 / d ** 2 + kap / d)
        rows.append(dict(point=name, d=d, h=h, far=far, v_solver=sv["v"], v_ls=v_ls, rel_diff=(v_ls - sv["v"]) / sv["v"]
                         if sv["v"] else None, resid_ls=resid, M_flux_box=Mf, M_integral=sv["M_particle"],
                         prediction_recomputed=pred, ratio_recomputed=v_ls / pred if pred else None))
    expo = {}
    for n in ("P1", "P2", "P3"):
        dd = np.array([22.0, 27.0, 33.0, 40.0])
        vv = np.array([abs(next(r["v_ls"] for r in rows if r["point"] == n and r["d"] == d and r["h"] == 0.1
                                and r["far"] == 150.0)) for d in dd])
        expo[n] = float(-np.polyfit(np.log(dd), np.log(vv), 1)[0])
    out = dict(M_radial_flux_r10=Mrf, M_radial_integral={n: rad[n]["M"] for n in rad}, pairs=rows, exponent_ls=expo)
    json.dump(out, open(os.path.join(OUT, "recompute.json"), "w"), indent=1, default=float)
    for r in rows:
        print(r["point"], r["d"], r["h"], r["far"], "v %.6e ls %.6e rel %s ratio %s  M box %.4f int %.4f" % (
            r["v_solver"], r["v_ls"], r["rel_diff"], r["ratio_recomputed"], r["M_flux_box"], r["M_integral"]))
    print("exponents", expo, "\nM radial flux", Mrf)


if __name__ == "__main__":
    main()
