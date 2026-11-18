"""EX005 check: sigma = 0.001 and 0.005 (kappa = 0) at finer grid (dr = 0.025) and larger domain; re-fit the far G tail
by a separate estimator (log-linear fit of r phi_G where the oscillating modes are below 1% of it)."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex005_model_g2 as m  # noqa: E402

out = []
for sigma in (0.001, 0.005):
    d = np.load(os.path.join(m.OUT, f"profile_k0_s{sigma:g}.npz"))
    beff = m.FRAC * m.b_c(sigma, kappa=0.0)
    rt = m.rates(beff, sigma, kappa=0.0)
    km = [x for x in rt if x["wavenumber"] < 1e-6][0]["decay"]
    for R, dr in ((d["r"][-1] + 0.025, 0.025), (min(200.0, 25 + 14 / km), 0.05)):
        rad = m.Radial2(R, dr, "neumann", beff, sigma, kappa=0.0)
        u = np.array([np.interp(rad.r, d["r"], d["u"][s] - d["h0"][s], right=0.0) + rad.h0[s] for s in range(3)])
        u, res, it, ok = rad.newton(u, tol=1e-9 * rad.h0[0] / 100, maxit=40)
        phi = u - rad.h0[:, None]
        r = rad.r
        lo = 22.0
        hi = min(R - 10, lo + 5 / km)
        mm = (r > lo) & (r < hi)
        y = r[mm] * phi[0][mm]
        c = np.polyfit(r[mm], np.log(np.abs(y)), 1) if (np.all(y < 0) or np.all(y > 0)) else [np.nan, np.nan]
        rec = dict(sigma=sigma, R=float(R), dr=dr, ok=bool(ok), res=float(res), int_G=rad.integral(phi[0]),
                   core_G=float(phi[0, 0]), sign_far=float(np.sign(y.mean())), kappa_fit=float(-c[0]),
                   kappa_linear=km, C_loglin=float(np.sign(y.mean()) * np.exp(c[1])),
                   halo_int=float(rad.integral(phi[0] * (r > 15))))
        print(rec, flush=True)
        out.append(rec)
json.dump(out, open(os.path.join(m.OUT, "check.json"), "w"), indent=1)
