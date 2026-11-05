"""MG9b post-hoc descriptive check (after ChatGPT MG035): event-shape mix by production cell.
Uses saved Part B run records only. Each scored first event is matched to its own run's confirmation (same sc).
Not part of the registered result. Usage: python src/shapes_mg9b.py"""
import os, json, glob, math
import numpy as np
from scipy import stats
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def wilson(k, n, z=1.959964):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


out = {}
for b in ("29.20", "29.35"):
    pb = json.load(open(f"results/MG9b/partB_{b}.json"))
    rows = []
    for c in pb["cells"]:
        Om = c["Om"]
        rs = sorted((json.load(open(f)) for f in glob.glob(f"results/MG9b/runs/b{b}_B_{Om:g}/*.json")),
                    key=lambda r: r["rep"])[: c["runs"]]
        nd = []
        for r in rs:
            if not (r["complete"] and r["event"]):
                continue
            m = [x for x in r["confirmations"] if x["sc"] == r["first"]["s_conf"]]
            assert len(m) == 1
            nd.append(m[0]["n_deep_end"])
        n, s = len(nd), sum(1 for k in nd if k == 1)
        lo, hi = wilson(s, n)
        rows.append(dict(Om=Om, events=n, single=s, multiple=n - s, pct=100 * s / n, wilson95=[100 * lo, 100 * hi],
                         n_deep_hist={str(k): nd.count(k) for k in sorted(set(nd))}))
    S, N = sum(r["single"] for r in rows), sum(r["events"] for r in rows)
    x = np.array([r["Om"] for r in rows]); k = np.array([r["single"] for r in rows]); n = np.array([r["events"] for r in rows])
    p = S / N
    z = np.sum(x * (k - n * p)) / math.sqrt(p * (1 - p) * (np.sum(n * x * x) - np.sum(n * x) ** 2 / N))
    out[b] = dict(cells=rows, single=S, events=N, trend_z=float(z), trend_p=float(2 * stats.norm.sf(abs(z))),
                  homog_p=float(stats.chi2_contingency(np.array([k, n - k]))[1]))
    for r in rows:
        print(b, f"{r['Om']:g}", r["events"], r["single"], r["multiple"], f"{r['pct']:.0f}%",
              "[%.0f, %.0f]" % tuple(r["wilson95"]), r["n_deep_hist"])
    print(b, "total", S, "/", N, "trend z %.2f p %.3f; homogeneity p %.3f" % (z, out[b]["trend_p"], out[b]["homog_p"]))
assert (out["29.20"]["single"], out["29.20"]["events"]) == (40, 75)
assert (out["29.35"]["single"], out["29.35"]["events"]) == (30, 76)
json.dump(out, open("results/MG9b/event_shapes_posthoc.json", "w"), indent=1)
