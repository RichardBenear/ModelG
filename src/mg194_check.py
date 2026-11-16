"""MG194 separate recomputation from saved trajectories (no use of the driver's fitted speeds or fit code).

- Re-fits every speed from trajectories.npz with its own least-squares code (Part A: t >= 100; Part B: t >= 20).
- Part A: recomputes the relative errors with the adjoint mu_a from results.json.
- Part B: recomputes the floor, resolved set, zeros and gaps; estimates k from the mean zero spacing and kappa from a
  log-linear regression of the per-half-wavelength peak |v_rel| (a different estimator from the driver's weighted fit).
Shares with the driver: the saved trajectories, results.json (for mu_a and run metadata) and NumPy.
"""
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "MG194_two_particle")


def slope(t, y):
    A = np.vstack([t, np.ones_like(t)]).T
    return float(np.linalg.lstsq(A, y, rcond=None)[0][0])


def main():
    rec = json.load(open(os.path.join(R, "results.json")))
    tr = np.load(os.path.join(R, "trajectories.npz"))
    mu = rec["soliton"]["mobilities"]["mu_a"]
    out = dict(partA=[], partB={})
    for r in rec["partA_runs"]:
        ts, xs, ds = tr[f"A_eps{r['eps']}_x0{r['x0']}"]
        sel = ts >= 100.0
        v = slope(ts[sel], xs[sel])
        dadx = 14.0 * r["eps"] * 2 * np.pi / 100.0 * np.cos(2 * np.pi * r["x0"] / 100.0)
        e = (v - mu * dadx) / (mu * dadx) if r["eps"] > 0 else None
        out["partA"].append(dict(eps=r["eps"], x0=r["x0"], v=v, rel_err=e, driver_v=r["v"],
                                 depth_dev=float(np.max(np.abs(ds - ds[0])) / abs(ds[0]))))
    out["A1_recomputed"] = all(abs(a["rel_err"]) <= 0.1 and a["depth_dev"] <= 0.05 for a in out["partA"] if a["eps"] > 0)
    runs = {(r["d"], r["L"]): r for r in rec["partB_runs"]}
    vb = {}
    for (d, L), r in runs.items():
        ts, sep, dep = tr[f"B_d{d}_L{int(L)}"]
        sel = ts >= 20.0
        vb[(d, L)] = slope(ts[sel], sep[sel])
    floor = max(10 * abs(vb[(0.0, 100.0)]), 2 * abs(vb[(30.0, 200.0)] - vb[(30.0, 100.0)]))
    ds_ = np.array(sorted(d for (d, L) in vb if L == 100.0 and d > 0))
    v = np.array([vb[(d, 100.0)] for d in ds_])
    merged = np.array([runs[(d, 100.0)]["merged"] for d in ds_])
    res = (~merged) & (np.abs(v) > floor)
    z = [float(ds_[i] - v[i] * (ds_[i + 1] - ds_[i]) / (v[i + 1] - v[i]))
         for i in range(len(ds_) - 1) if res[i] and res[i + 1] and v[i] * v[i + 1] < 0]
    gaps = np.diff(z)
    k_est = float(np.pi / np.mean(gaps)) if len(gaps) else None
    # peak |v| between consecutive zeros
    peaks = []
    edges = [ds_[0]] + z + [ds_[-1]]
    for lo, hi in zip(edges[:-1], edges[1:]):
        msk = res & (ds_ >= lo) & (ds_ <= hi)
        if msk.sum() >= 2:
            j = np.argmax(np.abs(v[msk]))
            peaks.append((ds_[msk][j], abs(v[msk][j])))
    kappa_est = None
    if len(peaks) >= 3:
        P = np.array(peaks)
        kappa_est = float(-slope(P[:, 0], np.log(P[:, 1])))
    out["partB"] = dict(floor=floor, n_resolved=int(res.sum()), zeros=z, gaps=gaps.tolist(), k_from_zeros=k_est,
                        kappa_from_peaks=kappa_est, n_peaks=len(peaks),
                        max_speed_diff_vs_driver=float(max(abs(vb[(d, 100.0)] - runs[(d, 100.0)]["v_rel"]) for d in ds_)),
                        all_attractive=bool(np.all(v[res] < 0)), all_repulsive=bool(np.all(v[res] > 0)))
    json.dump(out, open(os.path.join(R, "check.json"), "w"), indent=1, default=float)
    print(json.dumps(out["partB"], indent=1, default=float))
    print("A1 recomputed:", out["A1_recomputed"], [round(a["rel_err"], 4) for a in out["partA"] if a["eps"] > 0])


if __name__ == "__main__":
    main()
