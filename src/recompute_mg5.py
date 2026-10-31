"""MG5 independent recompute (proposal §9): classes, tracks and outcomes from the committed sparse check records only,
with detection and classification written separately from measure_graph / evolve (brute-force periodic distances,
union-find merging). Compares with results/MG5/evaluation.json.
    python src/recompute_mg5.py
"""
import os, sys, json, glob
import numpy as np

HERE = os.path.dirname(__file__)
OUT = os.environ.get("MG5_OUT") or os.path.join(HERE, "..", "results", "MG5")
CAND, MERGE, CENT, STEP, SPLIT, MOBILE = 1.0, 1.5, 0.6, 2.0, 4.0, 1.5


def wrap(d, L):
    return d - L * np.round(d / L)


def dist(p, q, Lx, Ly):
    return np.hypot(wrap(q[..., 0] - p[0], Lx), wrap(q[..., 1] - p[1], Ly))


def detect(pos, s, rows, thr, Lx, Ly):
    order = np.argsort(rows, kind="stable")                  # row order (records are stored that way already)
    pos, s, rows = pos[order], s[order], rows[order]
    cands = []
    for i in np.flatnonzero(s < thr):
        nb = np.flatnonzero(dist(pos[i], pos, Lx, Ly) <= CAND)
        m = s[nb].min()
        if s[i] < m or (s[i] == m and i == nb[s[nb] == m].min()):
            cands.append(i)
    if not cands:
        return []
    parent = list(range(len(cands)))

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a in range(len(cands)):
        for b in range(a + 1, len(cands)):
            if dist(pos[cands[a]], pos[[cands[b]]], Lx, Ly)[0] <= MERGE:
                parent[root(a)] = root(b)
    groups = {}
    for a, c in enumerate(cands):
        groups.setdefault(root(a), []).append(c)
    out = []
    for mem in groups.values():
        mem = np.array(mem)
        i0 = mem[np.argmin(s[mem])]
        nb = np.flatnonzero(dist(pos[i0], pos, Lx, Ly) <= CENT)
        w = np.maximum(thr - s[nb], 0.0)
        d = np.column_stack([wrap(pos[nb, 0] - pos[i0, 0], Lx), wrap(pos[nb, 1] - pos[i0, 1], Ly)])
        c = pos[i0] + (w[:, None] * d).sum(0) / w.sum()
        out.append(dict(centre=np.mod(c, (Lx, Ly)), depth=float(s[i0])))
    return sorted(out, key=lambda p: p["depth"])


def classify(lists, c0, Lx, Ly):
    n = [len(x) for x in lists]
    if any(n[i] >= 2 and n[i + 1] >= 2 for i in range(len(n) - 1)):
        return dict(cls="multiplied")
    if any(n[i] == 0 and n[i + 1] == 0 for i in range(len(n) - 1)):
        return dict(cls="dissolved")
    if not all(m == 1 for m in n):
        return dict(cls="unresolved")
    cur, u = np.array(c0, float), np.array(c0, float)
    md = 0.0
    for x in lists:
        st = np.array([wrap(x[0]["centre"][0] - cur[0], Lx), wrap(x[0]["centre"][1] - cur[1], Ly)])
        if np.hypot(*st) > STEP:
            return dict(cls="unresolved")
        cur = np.array(x[0]["centre"], float)
        u = u + st
        md = max(md, float(np.hypot(*(u - np.array(c0)))))
    return dict(cls="single", max_d=md, mobility="mobile" if md >= MOBILE else "pinned")


def main():
    ev_ = json.load(open(os.path.join(OUT, "evaluation.json")))
    bad = []
    for path in sorted(glob.glob(os.path.join(OUT, "runs", "*.json"))):
        rec = json.load(open(path))
        if rec["status"] != "COMPLETE":
            continue
        d = np.load(os.path.join(OUT, "checks", rec["key"] + ".npz"))
        Lx, Ly = float(d["Lx"]), float(d["Ly"])
        lists = []
        for i in range(len(d["t"])):
            a, b = d["off"][i], d["off"][i + 1]
            lists.append(detect(d["pos"][a:b], d["s"][a:b], d["rows"][a:b], rec["theta"], Lx, Ly))
        r = classify(lists[1:], lists[0][0]["centre"], Lx, Ly) if len(lists[0]) == 1 else dict(cls="unresolved")
        e = ev_["runs"][rec["key"]]["cls"]
        same = r["cls"] == e["cls"] and (r["cls"] != "single" or (r["mobility"] == e["mobility"]
                                                                    and abs(r["max_d"] - e["max_d"]) < 1e-9))
        print(f"{rec['key']}: {r['cls']} {r.get('mobility', '')} {'agrees' if same else 'DIFFERS'}")
        if not same:
            bad.append(rec["key"])
    print("All classes agree with the evaluator." if not bad else f"DIFFER: {bad}")


if __name__ == "__main__":
    main()
