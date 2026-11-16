"""EX003 (exploratory): screening range of the particle's G field as q -> g p, and the q = 0 case.

Family: p = 1, g = 0.1, q = 0.1 + 0.9 eps, a = 14 eps (keeps G0 = X0 = 140/9), b = (29/30.04408) * b_c(q).
eps = 1 is the published point. The 3D radial solver is src/ex002_far_field_3d.py with its module constants set.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex002_far_field_3d as e  # noqa: E402

ROOT = e.ROOT
OUT = os.path.join(ROOT, "results", "EX003_screening")
X0 = 140 / 9
D = np.diag([1.0, 1.0, 12.0])


def jac_h(q, b, g=0.1, p=1.0):
    Y0 = b / X0
    return np.array([[-q, g, 0], [p, -(1 + b) + 2 * X0 * Y0, X0 ** 2], [0, b - 2 * X0 * Y0, -X0 ** 2]])


def b_c(q):
    ks = np.linspace(0.01, 4, 800)

    def grow(b):
        return max(np.max(np.linalg.eigvals(jac_h(q, b) - D * k * k).real) for k in ks)
    lo, hi = 20.0, 40.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if grow(mid) > 0 else (mid, hi)
    return lo


def rates(q, b):
    ss = np.arange(4.0)
    c = np.polyfit(ss, [np.linalg.det(jac_h(q, b) + D * s) for s in ss], 3)
    out = []
    for s in np.roots(c):
        k = np.sqrt(complex(s))
        out.append(dict(decay=float(abs(k.real)), wavenumber=float(abs(k.imag))))
    return sorted(out, key=lambda d: d["decay"])


def params(eps):
    q = 0.1 + 0.9 * eps
    bc = b_c(q)
    return q, 14.0 * eps, 29.0 / 30.04408 * bc, bc


def set_params(q, a):
    e.Q_, e.A = q, a


def slow_tail(r, phi, kap_s, lo, hi):
    """Fit r*phi = C e^{-kap r} on [lo, hi] with kap free; returns C, kap."""
    m = (r >= lo) & (r <= hi)
    y = r[m] * phi[m]
    if np.all(np.abs(y) < 1e-300):
        return dict(C=0.0, kappa=None)
    rr = least_squares(lambda p: (p[0] * np.exp(-p[1] * r[m]) - y) / (np.abs(y) + 1e-300),
                       [y[0] * np.exp(kap_s * r[m][0]), kap_s])
    return dict(C=float(rr.x[0]), kappa=float(rr.x[1]))


def main(q0_only=False):
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    if q0_only:
        rec = json.load(open(os.path.join(OUT, "results.json")))
        return q0_part(rec, t0)
    rec = dict(family="p=1, g=0.1, q=0.1+0.9 eps, a=14 eps, b=(29/30.04408) b_c(q)", linear=[], particles=[])
    for eps in (1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001):
        q, a, b, bc = params(eps)
        rt = rates(q, b)
        rec["linear"].append(dict(eps=eps, q=q, a=a, b=b, b_c=bc, rates=rt))
        print(f"eps {eps}: q {q:.4f} a {a:.4f} b {b:.4f} rates {[(round(x['decay'], 5), round(x['wavenumber'], 5)) for x in rt]}",
              flush=True)
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)

    # 3D particles by continuation in eps
    d = np.load(os.path.join(ROOT, "results", "EX002_far_field_3d", "profile_R60_neumann_std.npz"))
    r0, u0 = d["r"], d["u"]
    targets = {1.0, 0.3, 0.1, 0.03, 0.01}
    path = sorted(set(np.round(np.geomspace(1.0, 0.01, 25), 12).tolist()) | targets, reverse=True)
    prev_r, prev_u = r0, u0
    for eps in path:
        q, a, b, bc = params(eps)
        set_params(q, a)
        kap_s = rates(q, b)[0]["decay"]
        R = float(min(200.0, max(60.0, 12.0 / kap_s)))
        rad = e.Radial(R, 0.02, "neumann", b)
        ui = np.array([np.interp(rad.r, prev_r, prev_u[s] - prev_u[s][-1]) + rad.h0[s] for s in range(3)])
        ui, res, it, ok = rad.newton(ui, tol=1e-9, maxit=40)
        if not ok:
            ui = rad.settle(ui, dt=0.2, T=100.0)
            ui, res, it, ok = rad.newton(ui, tol=1e-9, maxit=40)
        phi = ui - rad.h0[:, None]
        alive = ok and phi[1, 0] < -1.0
        print(f"eps {eps:.4f}: R {R:.0f} ok {ok} res {res:.1e} core phiX {phi[1, 0]:.3f}", flush=True)
        if not alive:
            rec["particles"].append(dict(eps=float(eps), ok=bool(ok), alive=False, core_phiX=float(phi[1, 0])))
            break
        prev_r, prev_u = rad.r, ui
        if any(abs(eps - t) < 1e-9 for t in targets):
            r = rad.r
            lo = max(30.0, 0.0)   # oscillating tail (rate ~0.516) is below 1e-7 of the core by r = 30
            hi = min(R - 10.0, lo + 6.0 / kap_s)
            st = slow_tail(r, phi[0], kap_s, lo, hi) if hi > lo + 5 else dict(C=None, kappa=None)
            stX = slow_tail(r, phi[1], kap_s, lo, hi) if hi > lo + 5 else dict(C=None, kappa=None)
            p = dict(eps=float(eps), q=q, a=a, b=b, R=R, ok=bool(ok), alive=True, residual=res,
                     core_phiX=float(phi[1, 0]), core_phiG=float(phi[0, 0]),
                     int_G=rad.integral(phi[0]), int_X=rad.integral(phi[1]),
                     int_G_core=rad.integral(phi[0] * (r < 3)),
                     kappa_slow_linear=kap_s, slow_tail_G=st, slow_tail_X=stX,
                     rphiG_at_1_over_kappa=float(np.interp(1 / kap_s, r, r * phi[0])),
                     rphiG_at_30=float(np.interp(30.0, r, r * phi[0])))
            rec["particles"].append(p)
            np.savez_compressed(os.path.join(OUT, f"profile_eps{eps:.3g}.npz"), r=r, u=ui)
            print("   ", {k: v for k, v in p.items() if k not in ("q", "a", "b")}, flush=True)
            json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)

    return q0_part(rec, t0)


def q0_part(rec, t0):
    # q = 0: homogeneous ODE from the published ether state, with the -qG term removed (p G term kept)
    def rhs(t, y, q):
        G, X, Y = y
        b = 29.0
        return [-q * G + 0.1 * X + 14.0, G - (1 + b) * X + X * X * Y, b * X - X * X * Y]
    y0 = [X0, X0, 29.0 / X0]
    sol = solve_ivp(rhs, (0, 200), y0, args=(0.0,), method="Radau", rtol=1e-8, atol=1e-8, t_eval=[0, 1, 10, 50, 100, 200])
    rec["q0_homogeneous"] = dict(t=sol.t.tolist(), G=sol.y[0].tolist(), X=sol.y[1].tolist(), Y=sol.y[2].tolist(),
                                 note="published a, b, g, p; -qG removed; G0 formula a/(q - g p) = -140 (no positive "
                                      "steady state)")
    print("q=0 homogeneous:", {k: np.round(v, 3).tolist() for k, v in rec["q0_homogeneous"].items() if k != "note"},
          flush=True)
    rec["wall_seconds"] = time.time() - t0
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print("done", rec["wall_seconds"], flush=True)


if __name__ == "__main__":
    main(q0_only="--q0" in sys.argv)
