"""MG5 independent recompute (proposal §9): from the committed sparse check records only, with detection,
classification, tracks, D_p and outcome rules written separately from measure_graph / evolve / evaluate_mg5
(brute-force periodic distances, union-find merging). Compares with results/MG5/evaluation.json:
  per run: class, mobility, max displacement, unwrapped track (1e-9) and D_p (1e-12 relative plus 1e-15 absolute, for the frozen runs' rounding-level values);
  aggregate: Part A control validity, the gate, each rate's survival and mobility label, and the Part B outcome.
Exits with status 1 if anything differs (MG018 review item 1).
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
    path = [u.copy()]
    for x in lists:
        st = np.array([wrap(x[0]["centre"][0] - cur[0], Lx), wrap(x[0]["centre"][1] - cur[1], Ly)])
        if np.hypot(*st) > STEP:
            return dict(cls="unresolved")
        cur = np.array(x[0]["centre"], float)
        u = u + st
        path.append(u.copy())
        md = max(md, float(np.hypot(*(u - np.array(c0)))))
    return dict(cls="single", max_d=md, mobility="mobile" if md >= MOBILE else "pinned", track=np.array(path),
                D_p=dp(np.array(path)))


def dp(path, step_t=10.0, max_lag=10):
    """Descriptive D_p: mean squared displacement at lags 1..10 checks (10-100 time units), slope through the origin
    of MSD = 4 D_p tau by least squares."""
    num = den = 0.0
    for lag in range(1, max_lag + 1):
        if lag >= len(path):
            break
        tau = lag * step_t
        m = float(np.mean(np.sum((path[lag:] - path[:-lag]) ** 2, axis=1)))
        num += tau * m
        den += tau * tau
    return num / (4.0 * den) if den else None


def outcomes(res):
    """Proposal §8 rules, from {key: (part, nu, cls dict)} with every run COMPLETE."""
    def runs(part, nu):
        return [r for (p, n, r) in res.values() if p == part and n == nu]
    three = lambda xs, want: sum(x == want for x in xs) >= 3
    ctrl = runs("A", 0.0)
    valid = len(ctrl) == 4 and all(r["cls"] == "single" and r["max_d"] < 0.1 for r in ctrl)
    out = dict(validity="valid" if valid else "INVALID-CONTROL", rates={})
    for nu in (0.002, 0.02, 0.2):
        rs = runs("A", nu)
        surv = three([r["cls"] for r in rs], "single")
        mob = [r["mobility"] for r in rs if r["cls"] == "single"]
        out["rates"][f"{nu:g}"] = dict(survives=surv, mobility=None if not surv else "PINNED" if three(mob, "pinned")
                                       else "MOBILE" if three(mob, "mobile") else "MIXED")
    n = sum(v["survives"] for v in out["rates"].values())
    out["gate"] = out["validity"] if not valid else "PASS" if n == 3 else "FAIL" if n == 0 else "PARTIAL"
    fz, fa = [r["cls"] for r in runs("B", 0.0)], [r["cls"] for r in runs("B", 0.2)]
    out["partB"] = ("A-void" if not three(fz, "multiplied") else "A1" if three(fa, "single")
                    else "A2" if three(fa, "multiplied") else "A3")
    return out


def main():
    ev_ = json.load(open(os.path.join(OUT, "evaluation.json")))
    bad, res, n_complete = [], {}, 0
    for path in sorted(glob.glob(os.path.join(OUT, "runs", "*.json"))):
        rec = json.load(open(path))
        if rec["status"] != "COMPLETE":
            continue
        with np.load(os.path.join(OUT, "checks", rec["key"] + ".npz")) as z:
            d = {k: z[k] for k in z.files}
        Lx, Ly = float(d["Lx"]), float(d["Ly"])
        lists = []
        for i in range(len(d["t"])):
            a, b = d["off"][i], d["off"][i + 1]
            lists.append(detect(d["pos"][a:b], d["s"][a:b], d["rows"][a:b], rec["theta"], Lx, Ly))
        r = classify(lists[1:], lists[0][0]["centre"], Lx, Ly) if len(lists[0]) == 1 else dict(cls="unresolved")
        n_complete += 1
        res[rec["key"]] = (rec["spec"]["part"], rec["spec"]["nu"], r)
        e = ev_["runs"][rec["key"]]["cls"]
        same = r["cls"] == e["cls"]
        if same and r["cls"] == "single":
            et = np.array(e["track"], float)
            same = (r["mobility"] == e["mobility"] and abs(r["max_d"] - e["max_d"]) < 1e-9
                    and et.shape == r["track"].shape and float(np.abs(et - r["track"]).max()) < 1e-9
                    and (e.get("D_p") is None) == (r["D_p"] is None)
                    and (r["D_p"] is None or abs(r["D_p"] - e["D_p"]) <= 1e-15 + 1e-12 * abs(e["D_p"])))
        dp_txt = f" D_p {r['D_p']:.6e}" if r.get("D_p") is not None else ""
        print(f"{rec['key']}: {r['cls']} {r.get('mobility', '')}{dp_txt} {'agrees' if same else 'DIFFERS'}")
        if not same:
            bad.append(rec["key"])
    if n_complete == 24:
        o = outcomes(res)
        ea, eb = ev_["partA"], ev_["partB"]
        agg = [("validity", o["validity"], ea["validity"]), ("gate", o["gate"], ev_["gate"]),
               ("Part B", o["partB"], eb.get("outcome"))]
        for nu, v in o["rates"].items():
            er = ea["rates"][nu]
            agg += [(f"nu={nu} survives", v["survives"], er.get("survives")),
                    (f"nu={nu} mobility", v["mobility"], er.get("mobility"))]
        for name, mine, theirs in agg:
            ok = mine == theirs
            print(f"{name}: {mine} {'agrees' if ok else f'DIFFERS (evaluator: {theirs})'}")
            if not ok:
                bad.append(name)
    else:
        print(f"{n_complete} of 24 runs COMPLETE: aggregate outcomes not recomputed")
    print("Everything agrees with the evaluator." if not bad else f"DIFFER: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
