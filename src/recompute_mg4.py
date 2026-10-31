"""MG4 independent recompute (CLAUDE.md: recompute before stating a result). Re-derives the outcome classes (C1) and the
C2 shape errors from the saved fields with separately written code: brute-force periodic distances instead of k-d trees,
union-find instead of sparse components, np.histogram for the radial bins. Graphs and smoothing come from the
registered builders. Prints the counts and medians to compare with results/MG4/evaluation.json.
    python src/recompute_mg4.py
"""
import os, sys, json, glob
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import run_mg4_partA as pa


def pdist(pos, p, Lx, Ly):
    d = pos - p
    d[:, 0] -= Lx * np.round(d[:, 0] / Lx)
    d[:, 1] -= Ly * np.round(d[:, 1] / Ly)
    return np.hypot(d[:, 0], d[:, 1])


def find(par, i):
    while par[i] != i:
        par[i] = par[par[i]]
        i = par[i]
    return i


def count(g, s, thr):
    pos, Lx, Ly = g["pos"], g["Lx"], g["Ly"]
    below = np.flatnonzero(s < thr)
    cand = [i for i in below if s[i] <= s[pdist(pos, pos[i], Lx, Ly) <= 1.0].min() and
            i == min(j for j in np.flatnonzero(pdist(pos, pos[i], Lx, Ly) <= 1.0) if s[j] == s[i])]
    par = list(range(len(cand)))
    for a in range(len(cand)):
        d = pdist(pos[cand], pos[cand[a]], Lx, Ly)
        for b in np.flatnonzero(d <= 1.5):
            ra, rb = find(par, a), find(par, int(b))
            if ra != rb:
                par[ra] = rb
    groups = {}
    for a, i in enumerate(cand):
        groups.setdefault(find(par, a), []).append(i)
    return [min(m, key=lambda i: s[i]) for m in groups.values()]


def main():
    dec = json.load(open(pa.DEC))
    theta = json.load(open(os.path.join(pa.OUT, "threshold.json")))["theta"]
    cls, prof = {}, {}
    graphs = {}
    edges = np.concatenate([[0.0], np.arange(0.5, 8.0 + 1e-9, 0.25)])
    for p in sorted(glob.glob(os.path.join(pa.OUT, "runs", "*.json"))):
        r = json.load(open(p))
        key = (r["fam"], r["slot"], r["torus"])
        if key not in graphs:
            graphs[key] = pa.get_graph(dict(fam=r["fam"], slot=r["slot"], torus=r["torus"]), dec)[0]
        g = graphs[key]
        f = np.load(pa.field_file(r["key"]))
        ns = [count(g, f["checks_sX"][k].astype(float), theta) for k in range(11)]
        n = [len(x) for x in ns]
        if max(n) >= 2:
            c = "multiplied"
        elif n[-1] == 0:
            c = "decayed"
        elif min(n) == 1:
            c = "single"          # (track-link test not repeated; the evaluator reports any 'unresolved')
        else:
            c = "unresolved"
        cls[r["key"]] = c
        if c == "single":
            N = g["N"]
            X = f["final"][1].astype(float)
            i0 = ns[-1][0]
            rr = pdist(g["pos"], g["pos"][i0], g["Lx"], g["Ly"])
            sm, _ = np.histogram(rr, edges, weights=X)
            ct, _ = np.histogram(rr, edges)
            prof[r["key"]] = np.where(ct > 0, sm / np.maximum(ct, 1), np.nan)
        print(r["key"], c, n[-1], flush=True)
    out = dict(classes=cls)
    for fam in ("D", "R"):
        per = [sum(cls.get(f"{fam}{s}_p{i}_sol") == "single" for i in range(8)) for s in range(4)]
        out[fam] = dict(single_total=sum(per), per_graph=per)
    rc = 0.5 * (edges[1:] + edges[:-1]); rc[0] = 0.25
    T = np.nanmean([prof[f"T_p{i}_sol"] for i in range(4)], axis=0)
    m = (rc < 6)
    for fam in ("D", "R"):
        e = []
        for k, pr in prof.items():
            if k.startswith(fam):
                ok = m & np.isfinite(pr) & np.isfinite(T)
                e.append(np.sqrt((rc[ok] * (pr[ok] - T[ok]) ** 2).sum() / (rc[ok] * T[ok] ** 2).sum()))
        out[fam]["median_L2_centre_node"] = float(np.median(e)) if e else None
        out[fam]["depth_rel_median"] = float(np.median([abs(prof[k][0] / T[0] - 1) for k in prof if k.startswith(fam)])) if e else None
    print(json.dumps({k: v for k, v in out.items() if k != "classes"}, indent=1))
    json.dump(out, open(os.path.join(pa.ROOT, "recompute.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
