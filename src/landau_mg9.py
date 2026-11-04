"""MG9 Part A1/A2 (spec/MG9_proposal.md draft 3, registered at 7810b8b): linear coefficient, cubic Landau coefficient,
Swift-Hohenberg code check, and the V2 fit.

Conventions (proposal §5): M(kappa) = J - kappa*D, kappa the positive Laplacian eigenvalue (k^2 in the continuum,
(4/h^2) sin^2(kh/2) on a lattice). Real cosine amplitude A of phi_X in original time t:
    dA/dt = sigma A + L3 A^3,   L3 = psi . [2 Q(phi, w0) + Q(phi, w2) + 3/4 C(phi)],
    w0 = -1/2 M(0)^-1 Q(phi, phi),  w2 = -1/2 M(kappa(2k))^-1 Q(phi, phi),
phi the right null vector of M(kappa_c) with phi_X = 1, psi the left null vector with psi.phi = 1.
L3 > 0: locally SUBCRITICAL.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

P17 = dict(a=14.0, b=29.0, g=0.1, p=1.0, q=1.0, dx=1.0, dy=12.0)
X0 = 140.0 / 9.0                     # G0 = X0 = a / (q - g p) with eq. 17; independent of b
LAM = 100.0 / 34.0                   # frozen physical pattern wavelength (mode 34 of L = 100)
K34 = 2 * np.pi * 34 / 100.0
DOMAINS = {"D02": dict(h=0.2, N=500, L=100.0, mode=34),
           "D02L": dict(h=0.2, N=1000, L=200.0, mode=68),
           "D01": dict(h=0.1, N=1000, L=100.0, mode=34)}
BRENT = dict(xtol=1e-12, rtol=4 * np.finfo(float).eps, maxiter=200)


def Y0(b):
    return b / X0


def jacobian(b, P=P17):
    """Reaction Jacobian at the uniform state (s = u = w = 0), order (G, X, Y)."""
    g, p, q = P["g"], P["p"], P["q"]
    y0 = Y0(b)
    return np.array([[-q, g, 0.0],
                     [p, -(1 + b) + 2 * X0 * y0, X0 ** 2],
                     [0.0, b - 2 * X0 * y0, -X0 ** 2]])


DMAT = np.diag([P17.get("dg", 1.0), P17["dx"], P17["dy"]])


def kappa(k, h=None):
    return k * k if h is None else (4.0 / h ** 2) * np.sin(k * h / 2.0) ** 2


def M(b, kap):
    return jacobian(b) - kap * DMAT


def sigma(b, kap):
    """Largest real part of the eigenvalues of M(b, kappa)."""
    return float(np.linalg.eigvals(M(b, kap)).real.max())


# ------------------------------------------------------------------ thresholds (A1, V1)
def sigma_max_continuous(b):
    """max over kappa >= 0 of sigma (continuous k; any lattice whose kappa range covers the optimum gives the same)."""
    r = minimize_scalar(lambda kap: -sigma(b, kap), bounds=(0.5, 12.0), method="bounded",
                        options=dict(xatol=1e-12))
    return -r.fun, r.x


def bc_continuous():
    bc = brentq(lambda b: sigma_max_continuous(b)[0], 29.0, 31.0, **BRENT)
    return bc, sigma_max_continuous(bc)[1]


def k_from_kappa(kap, h=None):
    return np.sqrt(kap) if h is None else (2.0 / h) * np.arcsin(np.sqrt(kap) * h / 2.0)


def bc_mode(k, h=None):
    """Threshold b at which the mode of wavenumber k (on mesh h, None = continuum) has sigma = 0."""
    kap = kappa(k, h)
    return brentq(lambda b: sigma(b, kap), 29.0, 40.0, **BRENT)


def box_thresholds(dom):
    d = DOMAINS[dom]
    rows = []
    for n in range(1, d["N"] // 2):
        k = 2 * np.pi * n / d["L"]
        kap = kappa(k, d["h"])
        if sigma(40.0, kap) <= 0:
            continue
        rows.append((bc_mode(k, d["h"]), n))
    onset, n_on = min(rows)
    sel = bc_mode(2 * np.pi * d["mode"] / d["L"], d["h"])
    return dict(onset=onset, onset_mode=n_on, selected_mode=d["mode"], selected_threshold=sel)


def alpha(k, h=None, db=1e-5):
    kap = kappa(k, h)
    b0 = bc_mode(k, h)
    return (sigma(b0 + db, kap) - sigma(b0 - db, kap)) / (2 * db)


# ------------------------------------------------------------------ cubic coefficient (A2)
def null_vectors(Mc):
    """Right and left null vectors of Mc (real), phi normalized phi[1] = 1, psi.phi = 1."""
    w, V = np.linalg.eig(Mc)
    i = int(np.argmin(np.abs(w)))
    phi = np.real(V[:, i])
    wl, U = np.linalg.eig(Mc.T)
    j = int(np.argmin(np.abs(wl)))
    psi = np.real(U[:, j])
    return phi, psi, float(np.abs(w[i]))


def cubic_coefficient(Mc, M0, M2, Q, C, norm_index=1):
    """Generic real-cosine cubic coefficient; returns L3 and its three terms."""
    phi, psi, resid = null_vectors(Mc)
    phi = phi / phi[norm_index]
    psi = psi / float(psi @ phi)
    qpp = Q(phi, phi)
    w0 = -0.5 * np.linalg.solve(M0, qpp)
    w2 = -0.5 * np.linalg.solve(M2, qpp)
    t1 = float(psi @ (2 * Q(phi, w0)))
    t2 = float(psi @ Q(phi, w2))
    t3 = float(psi @ (0.75 * C(phi)))
    return dict(L3=t1 + t2 + t3, terms=(t1, t2, t3), phi=phi, psi=psi, w0=w0, w2=w2, null_residual=resid)


E_XY = np.array([0.0, 1.0, -1.0])


def Q_mg(b):
    y0 = Y0(b)
    return lambda u, v: E_XY * (y0 * u[1] * v[1] + X0 * (u[1] * v[2] + u[2] * v[1]))


def C_mg(u):
    return E_XY * (u[1] ** 2 * u[2])


def l3_modelg(problem):
    """Problems of §5: 'i' D02 mode 34; 'ii' h = 0.2 continuous k_c; 'iii' continuum; 'iv' D01 mode 34."""
    if problem == "i":
        h, k = 0.2, K34
        b = bc_mode(k, h)
    elif problem == "iv":
        h, k = 0.1, K34
        b = bc_mode(k, h)
    elif problem in ("ii", "iii"):
        b, kap_c = bc_continuous()
        h = 0.2 if problem == "ii" else None
        k = k_from_kappa(kap_c, h)
    else:
        raise ValueError(problem)
    kap1, kap2 = kappa(k, h), kappa(2 * k, h)
    out = cubic_coefficient(M(b, kap1), M(b, 0.0), M(b, kap2), Q_mg(b), C_mg)
    out.update(problem=problem, b=b, k=k, h=h, kappa=kap1, kappa2=kap2)
    return out


def l3_swift_hohenberg(c2):
    """u_t = r u - (1 + d^2)^2 u + c2 u^2 - u^3 at r = 0, k = 1. Expected 19 c2^2/18 - 3/4."""
    Ms = lambda kap: np.array([[-(1.0 - kap) ** 2]])
    Q = lambda u, v: np.array([c2 * u[0] * v[0]])
    C = lambda u: np.array([-u[0] ** 3])
    return cubic_coefficient(Ms(1.0), Ms(0.0), Ms(4.0), Q, C, norm_index=0)["L3"]


def p1_outcome(code_check_ok, v2_ok, l3s, terms_i):
    """§5 ordered P1 outcome. l3s: dict problem -> L3; terms_i: the three terms of problem (i)."""
    if not (code_check_ok and v2_ok):
        return "NOT EVALUATED"
    if abs(l3s["i"]) < 1e-3 * sum(abs(t) for t in terms_i):
        return "DEGENERATE"
    if len({np.sign(v) for v in l3s.values()}) > 1:
        return "MIXED"
    return "SUBCRITICAL" if l3s["i"] > 0 else "SUPERCRITICAL"


# ------------------------------------------------------------------ V2 fit
def fit_l3(t, A, sig, t0=20.0):
    """Least squares of the integrated form on y = 1/A^2 over t >= t0, sigma fixed:
       y(t) = y0 e^{-2 sig tau} + L3 (e^{-2 sig tau} - 1)/sig,  tau = t - t0 (limit -2 tau as sig -> 0)."""
    t, A = np.asarray(t, float), np.asarray(A, float)
    m = t >= t0
    tau, y = t[m] - t0, 1.0 / A[m] ** 2
    e = np.exp(-2 * sig * tau)
    g = np.expm1(-2 * sig * tau) / sig if sig != 0 else -2 * tau
    X = np.column_stack([e, g])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    rel_change = float(abs(y[-1] - y[0]) / y[0])
    return dict(L3=float(coef[1]), y0=float(coef[0]), n=int(m.sum()), rel_change_invA2=rel_change)


def v2_decision(runs, l3_derived, sig, repeat=None):
    """runs: list of dict(A0, t, A, stop_t). Returns pass flag and details per §5."""
    usable = []
    for r in runs:
        if r["stop_t"] < 40.0:
            r["status"] = "unusable"
            continue
        f = fit_l3(r["t"], r["A"], sig)
        r.update(fit=f)
        ident = abs(f["L3"]) * r["A0"] ** 2 >= 100 * abs(sig) and f["rel_change_invA2"] >= 0.25
        r["status"] = "identifiable" if ident else "not identifiable"
        if ident:
            usable.append(f["L3"])
    ok = len(usable) >= 2
    if ok:
        ok = all(np.sign(v) == np.sign(l3_derived) and abs(v - l3_derived) <= 0.10 * abs(l3_derived) for v in usable)
        mean = np.mean(usable)
        ok = ok and (max(usable) - min(usable)) <= 0.10 * abs(mean)
    rep_ok = None
    if repeat is not None:
        base, rep = repeat
        rep_ok = abs(rep - base) <= 0.02 * abs(base)
        ok = ok and rep_ok
    return dict(passed=bool(ok), n_identifiable=len(usable), repeat_ok=rep_ok)


if __name__ == "__main__":
    import json
    out = {}
    bc, kc = bc_continuous()
    out["continuous"] = dict(bc=bc, kappa_c=kc)
    for d in DOMAINS:
        out[d] = box_thresholds(d)
    print(json.dumps(out, indent=1))
