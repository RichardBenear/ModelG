"""EXPLORATORY post-result diagnostic of MG9c (MG030): is the late drift a coherent translation, an amplitude
redistribution, or a movement of the envelope threshold crossings? Following `review/MG030_amplitude_center_prompt.md`.
No new integration; reads only results/MG9c/snapshots.npz and the run tracks. Changes no MG9c classification.
    python src/diag_mg9c_centre.py        -> results/MG9c/diagnostics/centre.json, centre.png

Conventions
- x_j = (j - N/2) h on D02 (N = 500, h = 0.2), x in [-50, 50). The patch sits near x = 0 (about +-29), so the periodic
  cut at x = +-50 lies in the uniform background; no unwrapping is needed (checked: weight at |x| > 45 is reported).
- C_edge = (x_L + x_R)/2 from the tracked envelope half-level crossings (the MG9c tracker, read from the run record at
  the snapshot times).
- C_amp = sum x E^2 / sum E^2, with E the frozen MG027 envelope of phi_X, over the full domain (no threshold mask).
  C_phi2 = the same with phi_X^2 (sensitivity check). Totals and weighted variances are reported. They are amplitude
  diagnostics, not conserved densities.
- Minima: local minima of phi_X with prominence >= E_pat/2 (as cell_count), sub-grid by a 3-point parabola. Interior
  minima = all but the outermost one on each side. Matched to the first snapshot by nearest position (|dx| < lam/4,
  else flagged ambiguous).
- Alignment: the shift s minimising ||phi_X(t2)(x) - phi_X(t1)(x - s)|| (Fourier translation, periodic), searched in
  [-lam/4, lam/4] around 0 (continuity), between successive snapshots; cumulative shift = sum. Remainder =
  ||phi(t2) - shifted phi(t1)|| / ||phi(t2) - phi(t1)|| (fraction of the change not explained by translation) and
  ||...|| / ||phi(t2)||.
"""
import os, sys, json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

sys.path.insert(0, os.path.dirname(__file__))
import run_mg9 as R
import run_mg9c as C

N, H = 500, 0.2
LAM = C.LAM
X = (np.arange(N) - N // 2) * H
OUTD = os.path.join(C.OUT, "diagnostics")
BS = ("29.25000", "29.30000", "29.35000", "29.22000")
TIMES = (1800, 3600, 6000, 12000, 15000, 30000)


def shift1(f, s):
    k = 2 * np.pi * np.fft.rfftfreq(N, d=H)
    return np.fft.irfft(np.fft.rfft(f) * np.exp(-1j * k * s), n=N)


def centres(phiX, E_pat):
    E = R.envelope(phiX, H)
    wE, w2 = E ** 2, phiX ** 2
    out = {}
    for name, w in (("amp", wE), ("phi2", w2)):
        c = float(np.sum(X * w) / np.sum(w))
        out[name] = dict(C=c, total=float(np.sum(w) * H), var=float(np.sum((X - c) ** 2 * w) / np.sum(w)),
                         edge_weight_frac=float(np.sum(w[np.abs(X) > 45]) / np.sum(w)))
    return out


def minima(phiX, E_pat):
    ext = np.concatenate([phiX, phiX, phiX])
    pk, _ = find_peaks(-ext, prominence=0.5 * E_pat, distance=max(1, int(np.ceil(LAM / 2 / H))))
    pk = pk[(pk >= N) & (pk < 2 * N)] - N
    xs = []
    for j in pk:
        a, b, c = phiX[(j - 1) % N], phiX[j], phiX[(j + 1) % N]
        d = a - 2 * b + c
        off = 0.5 * (a - c) / d if d != 0 else 0.0
        xs.append(X[j] + off * H)
    return np.sort(np.array(xs))


def align(f1, f2):
    obj = lambda s: float(np.sum((f2 - shift1(f1, s)) ** 2))
    r = minimize_scalar(obj, bounds=(-LAM / 4, LAM / 4), method="bounded", options=dict(xatol=1e-10))
    s = float(r.x)
    res = np.sqrt(obj(s))
    dn = np.linalg.norm(f2 - f1)
    edge = abs(abs(s) - LAM / 4) < 1e-3
    return dict(shift=s, rem_of_change=float(res / dn) if dn > 0 else 0.0,
                rem_of_field=float(res / np.linalg.norm(f2)), at_bound=bool(edge))


def analyse(series, edges, E_pat):
    """series: dict t -> phiX; edges: dict t -> (xL, xR)."""
    ts = sorted(series)
    rows, m0 = [], minima(series[ts[0]], E_pat)
    interior0 = m0[1:-1]
    cum, prev = 0.0, None
    for t in ts:
        f = series[t]
        c = centres(f, E_pat)
        m = minima(f, E_pat)
        d = []
        amb = 0
        for x0 in interior0:
            j = int(np.argmin(np.abs(m - x0)))
            if abs(m[j] - x0) < LAM / 4:
                d.append(m[j] - x0)
            else:
                amb += 1
        al = align(series[prev], f) if prev is not None else dict(shift=0.0, rem_of_change=None, rem_of_field=None,
                                                                    at_bound=False)
        cum += al["shift"]
        xl, xr = edges[t]
        rows.append(dict(t=t, C_edge=(xl + xr) / 2, width=xr - xl, C_amp=c["amp"]["C"], C_phi2=c["phi2"]["C"],
                         amp=c["amp"], phi2=c["phi2"], n_min=len(m),
                         interior_min_x=[float(v) for v in m],
                         interior_min_shift_mean=float(np.mean(d)) if d else None,
                         interior_min_shift_sd=float(np.std(d)) if d else None, ambiguous_minima=amb,
                         align_step=al, align_cum=cum))
        prev = t
    r0 = rows[0]
    for r in rows:
        for k in ("C_edge", "C_amp", "C_phi2", "width"):
            r["d" + k] = r[k] - r0[k]
    return rows


def late(rows, a=15000, b=30000):
    A = next(r for r in rows if r["t"] == a)
    B = next(r for r in rows if r["t"] == b)
    out = {k: B[k] - A[k] for k in ("C_edge", "C_amp", "C_phi2", "width")}
    out["interior_minima"] = (B["interior_min_shift_mean"] - A["interior_min_shift_mean"]
                              if A["interior_min_shift_mean"] is not None else None)
    ma, mb = np.array(A["interior_min_x"]), np.array(B["interior_min_x"])
    if len(ma) == len(mb):
        dd = (mb - ma)[1:-1]                                      # interior minima, matched by order
        out["interior_minima_late_sd"] = float(np.std(dd))
        out["interior_minima_late_range"] = [float(dd.min()), float(dd.max())]
    out["alignment"] = B["align_cum"] - A["align_cum"]
    out["alignment_rem_of_change"] = B["align_step"]["rem_of_change"]
    return out


# ------------------------------------------------------------------ synthetic validation
def synth(a=29.0, s=0.0, amp=lambda x: 1.0, phase=0.0):
    k = 2 * np.pi / LAM
    w = 0.5 * (np.tanh((X - s + a) / (LAM / 2)) - np.tanh((X - s - a) / (LAM / 2)))
    return -6.5 * amp(X) * w * (0.5 + 0.5 * np.cos(k * (X - s) + phase)) - 0.3 * w


def validate():
    E_pat = 6.5
    base = synth()
    edges0 = (-29.0, 29.0)
    cases = {}
    # (i) known sub-grid translation by 0.0083 (the late drift size) and by 0.05
    for s in (0.0083, -0.0083, 0.05):
        r = analyse({0: base, 1: synth(s=s)}, {0: edges0, 1: (-29 + s, 29 + s)}, E_pat)[-1]
        cases[f"translate_{s:+g}"] = dict(true=s, C_amp=r["dC_amp"], C_phi2=r["dC_phi2"],
                                          minima=r["interior_min_shift_mean"], align=r["align_cum"])
    # (ii) symmetric amplitude relaxation, no displacement
    sym = synth(amp=lambda x: 1 + 0.02 * np.exp(-(x / 15) ** 2))
    r = analyse({0: base, 1: sym}, {0: edges0, 1: edges0}, E_pat)[-1]
    cases["symmetric_amplitude"] = dict(true=0.0, C_amp=r["dC_amp"], C_phi2=r["dC_phi2"],
                                        minima=r["interior_min_shift_mean"], align=r["align_cum"],
                                        align_rem_of_change=r["align_step"]["rem_of_change"])
    # (iii) asymmetric amplitude change with fixed peak positions: a weighted centroid moves without translation
    asym = synth(amp=lambda x: 1 + 0.002 * x / 29.0)
    r = analyse({0: base, 1: asym}, {0: edges0, 1: edges0}, E_pat)[-1]
    cases["asymmetric_amplitude"] = dict(true=0.0, C_amp=r["dC_amp"], C_phi2=r["dC_phi2"],
                                         minima=r["interior_min_shift_mean"], align=r["align_cum"],
                                         align_rem_of_change=r["align_step"]["rem_of_change"])
    return cases


def main():
    os.makedirs(OUTD, exist_ok=True)
    z = np.load(os.path.join(C.OUT, "snapshots.npz"))
    out = dict(note="EXPLORATORY post-result diagnostic; no integration; changes no MG9c classification. Snapshot "
                    "times are sparse: " + ", ".join(map(str, TIMES)) + ".",
               validation=validate(), runs={})
    for b in BS:
        rec = C.load_rec(f"front_D02_{b}")
        tr = rec["track"]
        tt = np.array([s["t"] for s in tr])
        series, edges = {}, {}
        for t in TIMES:
            series[t] = z[f"front_D02_{b}_t{t}"][N:2 * N]
            i = int(np.argmin(np.abs(tt - t)))
            edges[t] = (tr[i]["xL"], tr[i]["xR"])
        rows = analyse(series, edges, rec["E_pat"])
        out["runs"][b] = dict(rows=rows, late=late(rows), descriptive=(b == "29.22000"))
    with open(os.path.join(OUTD, "centre.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 4, figsize=(15, 3.6), sharey=True)
        for ax, b in zip(axs, BS):
            rows = out["runs"][b]["rows"]
            t = [r["t"] for r in rows]
            for k, lab in (("dC_edge", "edge midpoint"), ("dC_amp", "E²-weighted"), ("dC_phi2", "φ_X²-weighted")):
                ax.plot(t, [r[k] for r in rows], marker="o", label=lab)
            ax.plot(t, [r["interior_min_shift_mean"] for r in rows], marker="s", label="interior minima")
            ax.plot(t, [r["align_cum"] for r in rows], marker="^", ls="--", label="profile alignment")
            ax.set_title(f"b = {float(b):.2f}" + (" (descriptive)" if b == "29.22000" else ""))
            ax.set_xlabel("t")
            ax.axhline(0, color="0.8", lw=0.8)
        axs[0].set_ylabel("change in centre since t = 1800")
        axs[0].legend(fontsize=7)
        fig.suptitle("MG9c exploratory: centre changes (sparse snapshots)", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(OUTD, "centre.png"), dpi=110)
    except ImportError:
        pass
    return out


if __name__ == "__main__":
    o = main()
    print(json.dumps(o["validation"], indent=1, default=float))
    for b, v in o["runs"].items():
        print(b, {k: (round(x, 5) if isinstance(x, float) else x) for k, x in v["late"].items()})
