import os, sys, json, glob, csv
import numpy as np
sys.path.insert(0, "src")
import graphs as gr, graph_ode as go, measure_graph as mg, run_mg4_partA as pa, run_mg4_part0 as p0
from lattice import P17
R = "results/MG4c"
dec = json.load(open("results/MG4/part0/decisions.json"))
recs = {os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob(f"{R}/runs/*.json")}
tor = gr.torus(p0.H_T, p0.NX, p0.NY)
G = {}
def graph(r):
    k = (r["kind"], r["slot"])
    if k not in G:
        if r["kind"] == "T":
            g = gr.triangular(tor); h = float(dec["T"]["h"])
        else:
            pan = dec["panel"][f"D{r['slot']}"]; g = p0.build("D", pan["n_points"], tor["Lx"], tor["Ly"], pan["rep"]); h = float(pan["h"])
        G[k] = (g, h, mg._tree(g), mg.smoother(g))
    return G[k]
def cls(r, th):
    g, h, tree, S = graph(r)
    sX = np.load(f"{R}/fields/{r['key']}.npz")["checks_sX"].astype(float)
    return mg.classify(g, [mg.particles(g, sX[k], th, tree) for k in range(11)])[0]
out = {}
bad_bal = []
for k, r in sorted(recs.items()):
    g, h, tree, S = graph(r)
    f = np.load(f"{R}/fields/{k}.npz")
    y = f["final"].ravel()
    # independent: min of smoothed final from float64 final state
    m = float((S @ f["final"][1]).min())
    assert abs(m - r["min_smoothed_final"]) < 1e-9 * max(1, abs(m)), (k, m, r["min_smoothed_final"])
    op = w = None
    if r["op"] == "cotan":
        c = gr.cotan(g); op, w = c["L"], c["A"]
    rhs = go.make_system(g, h, {**P17, "b": r["b"]}, op=op)[0]
    d = rhs(0, y)
    bal = go.balances(y, d, g["N"], {**P17, "b": r["b"]}, w=w)
    if not bal["ok"]: bad_bal.append(k)
    assert abs(np.abs(d).max() - r["max_dphidt"]) <= 1e-6 * max(1e-12, r["max_dphidt"]) + 1e-12, (k, np.abs(d).max(), r["max_dphidt"])
# tables
TH29 = -5.117354421620525
bs = sorted({r["b"] for r in recs.values() if r["op"] is None}, reverse=True)
tab = {}
for b in bs:
    key = lambda s: f"b{b:.6f}_{s}"
    T = recs[key("T_p0")]
    mT = T["min_smoothed_final"]
    th = TH29; elig = False
    if cls(T, -1.0) == "single" and mT <= -2.0 and cls(T, 0.5 * mT) == "single":
        th, elig = 0.5 * mT, True
    D = [cls(recs[key(f"D{s}_p{p}")], th) for s, p in [(0,0),(0,1),(1,0),(1,1)]]
    D29 = [cls(recs[key(f"D{s}_p{p}")], TH29) for s, p in [(0,0),(0,1),(1,0),(1,1)]]
    tab[b] = (th, elig, D, D29)
    print(f"{b:.8f} th {th:.4f} elig {elig} {D} {'same@29' if D == D29 else D29}")
ev = json.load(open(f"{R}/evaluation.json"))
for b, (th, elig, D, _) in tab.items():
    e = ev["tables"][f"b{b:.6f}"]
    assert list(e["D"].values()) == D and e["T_eligible"] == elig and abs(e["theta"] - th) < 1e-12, b
p3 = {k: cls(recs[k], TH29) for k in ("w_b29.000000_D0_p0", "w_b29.000000_T_p0")}
yT = np.load(f"{R}/fields/w_b29.000000_T_p0.npz")["final"]; ref = np.load("results/MG4/part0/steps/Tstate_122x142.npy").reshape(yT.shape)
print("part3", p3, "maxdiff X", np.abs(yT[1] - ref[1]).max())
rows = list(csv.DictReader(open(f"{R}/budget_ledger.csv", encoding="utf-8")))
print("ledger rows", len(rows), "core-h", sum(float(r["wall_s"]) for r in rows) / 3600, "distinct", len({r["key"] for r in rows}))
print("float64 balance flags", bad_bal)
print("all-single b", [b for b, t in tab.items() if all(c == "single" for c in t[2])])
