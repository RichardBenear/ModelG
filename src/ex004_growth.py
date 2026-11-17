"""EX004 (exploratory): does growth give a Model G structure a net G deficit, and does it reach far near q -> g p?

Budget (any dimension, no boundary flux): d/dt I_G = -q I_G + g I_X,  d/dt S = p I_G - I_X + W,
with I_f = integral of phi_f, S = I_X + I_Y, W = integral of an added X source w(r).
Steady with a source: I_G = g W / (q - g p).  Slow growth (dI_G/dt ~ 0, W = 0): I_G = -g (dS/dt) / (q - g p).
A sink (W < 0) is the steady stand-in for growth; Arm C checks the solver against it, Arm G uses the model's own growth.
Family (EX003): q = 0.1 + 0.9 eps, a = 14 eps, b = (29/30.04408) b_c(q) unless stated.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex002_far_field_3d as e  # noqa: E402
import ex003_screening as x3  # noqa: E402

ROOT = e.ROOT
OUT = os.path.join(ROOT, "results", "EX004_growth")
G, P = 0.1, 1.0


class Src(e.Radial):
    """Radial solver with an extra X source w(r) (per unit volume)."""

    def __init__(self, *a, w=None, **k):
        super().__init__(*a, **k)
        self.w = np.zeros(self.n) if w is None else w

    def F(self, u):
        f = super().F(u)
        f[1] += self.w
        return f


def load_profile(eps, rad):
    d = np.load(os.path.join(ROOT, "results", "EX003_screening", f"profile_eps{eps:.3g}.npz"))
    r0, u0 = d["r"], d["u"]
    return np.array([np.interp(rad.r, r0, u0[s] - u0[s][-1]) + rad.h0[s] for s in range(3)])


def budget(rad, u):
    phi = u - rad.h0[:, None]
    IG, IX, IY = (rad.integral(phi[s]) for s in range(3))
    return phi, IG, IX, IY


def front(rad, phi, thr=1.0):
    idx = np.nonzero(np.abs(phi[1]) > thr)[0]
    return float(rad.r[idx[-1]]) if len(idx) else 0.0


def arm_c(eps, W=-1.0, sig=2.0):
    q, a, b, bc = x3.params(eps)
    x3.set_params(q, a)
    kap_s = x3.rates(q, b)[0]["decay"]
    R = float(max(60.0, 15.0 / kap_s))
    rad = Src(R, 0.05, "neumann", b)
    u0, res0, _, ok0 = rad.newton(load_profile(eps, rad), tol=1e-9, maxit=40)
    w = np.exp(-rad.r ** 2 / (2 * sig ** 2))
    rad.w = W * w / rad.integral(w)
    u1, res1, _, ok1 = rad.newton(u0.copy(), tol=1e-9, maxit=40)
    _, IG0, IX0, _ = budget(rad, u0)
    phi1, IG1, IX1, IY1 = budget(rad, u1)
    dG = u1[0] - u0[0]
    lo = 20.0
    hi = min(R - 10.0, lo + 6.0 / kap_s)
    tail = x3.slow_tail(rad.r, dG, kap_s, lo, hi)
    pred = G * W / (q - G * P)
    m = rad.r > 1.0 / kap_s
    out = dict(eps=eps, q=q, a=a, b=b, R=R, W=W, sigma=sig, ok=bool(ok0 and ok1), residuals=[res0, res1],
               IG_base=IG0, IG_with_sink=IG1, IX_with_sink=IX1, IG_pred=pred, rel_err=(IG1 - pred) / pred,
               kappa_slow_linear=kap_s, tail_dG=tail, frac_dG_beyond_1_over_kappa=rad.integral(dG * m) / rad.integral(dG),
               core_phiX=float(phi1[1, 0]))
    print("arm C", {k: v for k, v in out.items() if k not in ("q", "a")}, flush=True)
    return out


def fold(eps, R=60.0):
    q, a, b0, bc = x3.params(eps)
    x3.set_params(q, a)
    rad = e.Radial(R, 0.05, "neumann", b0)
    u, res, _, ok = rad.newton(load_profile(eps, rad), tol=1e-9, maxit=40)
    assert ok
    b, path = b0, []
    db = 0.1 * (bc - b0) / (30.04409 - 29.0)
    while db > 2e-4 * (bc - b0) and b < bc:
        r2 = e.Radial(R, 0.05, "neumann", b + db)
        ut, res, _, ok = r2.newton(u - rad.h0[:, None] + r2.h0[:, None], tol=1e-9, maxit=25)
        if ok and ut[1, 0] - r2.h0[1] < -1:
            b, u, rad = b + db, ut, r2
            path.append(dict(b=b, S=rad.integral(u[1] + u[2] - rad.h0[1] - rad.h0[2]),
                             front=front(rad, u - rad.h0[:, None])))
        else:
            db /= 2
    print(f"fold eps {eps}: b0 {b0:.5f} b_fold {b:.5f} b_c {bc:.5f}", flush=True)
    return dict(b0=b0, b_fold=b, b_c=bc, path=path), rad, u


def arm_g(eps, frac=0.3, R=200.0, dr=0.05, dt=0.2, T=400.0, tag=""):
    fd, rf, uf = fold(eps)
    q, a, _, bc = x3.params(eps)
    b = fd["b_fold"] + frac * (bc - fd["b_fold"])
    kap_s = x3.rates(q, b)[0]["decay"]
    rad = e.Radial(R, dr, "neumann", b)
    u = np.array([np.interp(rad.r, rf.r, uf[s] - uf[s][-1]) + rad.h0[s] for s in range(3)])
    I = sp.identity(3 * rad.n, format="csc")
    ts, rec, snaps = [], [], {}
    t = 0.0
    for k in range(int(round(T / dt)) + 1):
        if k % int(round(1.0 / dt)) == 0:
            phi, IG, IX, IY = budget(rad, u)
            fr = front(rad, phi)
            out_m = rad.r > fr + 2.0
            rec.append([t, fr, IG, IX, IY, rad.integral(phi[0] * out_m), float(phi[1, 0])])
            if k % int(round(20.0 / dt)) == 0:
                snaps[f"t{t:.0f}"] = phi.copy()
            if fr > R - 40.0 or phi[1, 0] > -1:
                break
        Fv = rad.F(u).ravel()
        u = u + spla.spsolve(I - dt * rad.J(u), dt * Fv).reshape(3, rad.n)
        t += dt
    A = np.array(rec)
    S = A[:, 3] + A[:, 4]
    dS = np.gradient(S, A[:, 0])
    pred_inst = -G * dS / (q - G * P)
    # window averages (20 time units) remove the shell-by-shell bursts
    win = []
    for t0 in np.arange(0, A[-1, 0] - 19.9, 20.0):
        m = (A[:, 0] >= t0) & (A[:, 0] <= t0 + 20)
        i0, i1 = np.nonzero(m)[0][[0, -1]]
        rate = (S[i1] - S[i0]) / (A[i1, 0] - A[i0, 0])
        win.append(dict(t0=float(t0), front0=float(A[i0, 1]), front1=float(A[i1, 1]), dSdt=float(rate),
                        IG_mean=float(np.mean(A[m, 2])), IG_pred=float(-G * rate / (q - G * P)),
                        IG_outside_mean=float(np.mean(A[m, 5]))))
    # far field outside the final front
    phi = u - rad.h0[:, None]
    fr = front(rad, phi)
    lo, hi = fr + 3.0, min(R - 5.0, fr + 3.0 + 6.0 / kap_s)
    tail = x3.slow_tail(rad.r, phi[0], kap_s, lo, hi) if hi > lo + 3 else dict(C=None, kappa=None)
    res = dict(eps=eps, q=q, a=a, b=b, b_c=bc, fold=fd, R=R, dr=dr, dt=dt, T_end=float(A[-1, 0]),
               kappa_slow_linear=kap_s, windows=win, final_front=fr, final_tail_G=tail,
               final_IG=float(A[-1, 2]), final_IG_outside=float(A[-1, 5]))
    np.savez_compressed(os.path.join(OUT, f"growth_eps{eps:.3g}{tag}.npz"), r=rad.r, series=A, dSdt=dS,
                        IG_pred_inst=pred_inst, **snaps, final=phi)
    print(f"arm G eps {eps}{tag}: b {b:.5f} T {A[-1, 0]:.0f} front {fr:.1f} final IG {A[-1, 2]:.3f} "
          f"outside {A[-1, 5]:.3f} tail {tail}", flush=True)
    for w_ in win:
        print("   ", {k: round(v, 4) for k, v in w_.items()}, flush=True)
    return res


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rec = dict(note="EXPLORATORY", armC=[], armG=[])
    for eps in (1.0, 0.03):
        rec["armC"].append(arm_c(eps))
        json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    for eps, kw in ((1.0, {}), (1.0, dict(dr=0.025, dt=0.1, tag="_fine")), (0.03, {}), (0.01, {})):
        rec["armG"].append(arm_g(eps, **kw))
        json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    rec["wall_seconds"] = time.time() - t0
    json.dump(rec, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)
    print("done", rec["wall_seconds"], flush=True)


if __name__ == "__main__":
    main()
