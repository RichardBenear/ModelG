"""MG9 evaluator (spec/MG9_proposal.md draft 3, §§5, 6.4-6.6). Reads results/MG9/ and writes results/MG9/evaluation.json.
    python src/evaluate_mg9.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import run_mg9 as R

TOL = 0.10
D01_BENCH = dict(b_s=(28.97815, 28.97827), b_r=(29.4117, 29.4121))       # MG003, h = 0.1 (not rechecked)


def side_bracket(seq, side):
    """For a D01 side: the adjacent pair (non-A outside, A inside). seq sorted (b, class)."""
    for (b1, c1), (b2, c2) in zip(seq, seq[1:]):
        if side == "lower" and c1 in ("M", "N") and c2 == "A":
            return (b1, b2)
        if side == "upper" and c1 == "A" and c2 in ("M", "N"):
            return (b1, b2)
    return None


def p2_outcome(sweep_start_ok, bench, ed, stages_complete=True):
    """§6.6 ordered outcome. bench: dict name -> dict(status, bracket); ed: run_mg9.final_edges()."""
    flags = []
    pts = ed.get("points", {})
    seq = R.grid_classes(pts)
    if not sweep_start_ok or not stages_complete:
        return "UNRESOLVED", flags, "sweep start failed or a required stage is incomplete"
    if any(bench.get(n, {}).get("status") != "OK" for n in ("b_s", "b_r")):
        return "UNRESOLVED", flags, "BENCHMARK-UNRESOLVED"
    if any(bench[n].get("sensitive") for n in ("b_s", "b_r")):
        flags.append("BENCHMARK-SENSITIVE")
    bs_lo, bs_hi = bench["b_s"]["bracket"]
    br_lo, br_hi = bench["b_r"]["bracket"]
    n_int = ed.get("intervals", 0)
    if n_int > 1:
        return "UNRESOLVED", flags, "disconnected arrest intervals"
    if n_int == 0:
        if any(c == "U" for _, c in seq):
            return "UNRESOLVED", flags, "no A point, some U points"
        return "NO ARREST OBSERVED", flags, "N points: " + ", ".join(f"{b}" for b, c in seq if c == "N")
    lo, up = ed["lower"], ed["upper"]
    if "unresolved" in (lo["status"], up["status"]):
        return "UNRESOLVED", flags, "an edge is unresolved"
    if "censored" in (lo["status"], up["status"]):
        flags.append("CENSORED")
    l2, u1 = lo["inner"], up["inner"]
    # ordering (ANOMALOUS): below the interval only RETREATS or N; above only INVADES
    for b, lab in pts.items():
        c = R.front_class(lab)
        if b < l2 and not (lab == "RETREATS" or c in ("N", "U")):
            flags.append("ANOMALOUS"); break
        if b > u1 and not (lab == "INVADES" or c == "U"):
            flags.append("ANOMALOUS"); break
    if u1 < bs_hi or l2 > br_lo:
        return "OUTSIDE", flags, ""
    if lo["status"] == "resolved" and up["status"] == "resolved":
        l1, u2 = lo["outer"], up["outer"]
        L_ok = bs_hi - TOL <= min(l1, l2) and max(l1, l2) <= bs_lo + TOL
        U_ok = br_hi - TOL <= min(u1, u2) and max(u1, u2) <= br_lo + TOL
        if L_ok and U_ok:
            return "COINCIDES", flags, ""
    return "OVERLAPS", flags, ""


def mesh_flag(ed, bench, d01_points):
    """MESH-SENSITIVE per §6.4: edge midpoints, widths and offsets from each mesh's own benchmarks."""
    if ed.get("intervals") != 1 or not d01_points:
        return "NOT EVALUATED"
    seq1 = R.grid_classes(d01_points)
    lb, ub = side_bracket(seq1, "lower"), side_bracket(seq1, "upper")
    if lb is None or ub is None or ed["lower"]["status"] != "resolved" or ed["upper"]["status"] != "resolved":
        return "NOT EVALUATED"
    m2 = lambda e: 0.5 * (e["inner"] + e["outer"])
    l02, u02 = m2(ed["lower"]), m2(ed["upper"])
    l01, u01 = 0.5 * sum(lb), 0.5 * sum(ub)
    bs02, br02 = 0.5 * sum(bench["b_s"]["bracket"]), 0.5 * sum(bench["b_r"]["bracket"])
    bs01, br01 = 0.5 * sum(D01_BENCH["b_s"]), 0.5 * sum(D01_BENCH["b_r"])
    w02, w01 = u02 - l02, u01 - l01
    sens = (abs(w01 - w02) > 0.5 * abs(w02) or abs((l01 - bs01) - (l02 - bs02)) > TOL
            or abs((u01 - br01) - (u02 - br02)) > TOL)
    return "MESH-SENSITIVE" if sens else "not set"


def main():
    out = {}
    p1p = os.path.join(R.OUT, "p1.json")
    out["P1"] = json.load(open(p1p))["P1"] if os.path.exists(p1p) else "NOT RUN"
    sweep = {b: R.load_rec(R._key("sweep", f"{b:.2f}")) for b in R.SWEEP_B}
    out["sweep_status"] = {f"{b:.2f}": (r["status"] if r else "NOT RUN") for b, r in sweep.items()}
    seqs = [out["sweep_status"][f"{b:.2f}"] for b in R.SWEEP_B]
    run = 0
    while run < len(seqs) and seqs[run] == "PERSISTS":
        run += 1
    out["persistence_bracket"] = (R.SWEEP_B[run] if run < len(seqs) else None, R.SWEEP_B[run - 1] if run else None)
    out["regained"] = [f"{b:.2f}" for i, b in enumerate(R.SWEEP_B) if i > run and seqs[i] == "PERSISTS"]
    sweep_ok = seqs[0] == "PERSISTS"
    bp = os.path.join(R.OUT, "benchmarks.json")
    bench = json.load(open(bp)) if os.path.exists(bp) else {}
    ed = R.final_edges()
    grid_done = all(R.load_rec(R._key("front", "grid", f"{b:.4f}")) for b in R.GRID_B)
    label, flags, note = p2_outcome(sweep_ok, bench, ed, stages_complete=grid_done and bool(bench))
    d01 = {}
    if os.path.isdir(R.RUNS):
        for fn in os.listdir(R.RUNS):
            if fn.startswith("front_D01_"):
                r = R.load_rec(fn[:-5]); d01[round(r["b"], 4)] = r["label"]
    mf = mesh_flag(ed, bench, d01) if bench else "NOT EVALUATED"
    if mf == "MESH-SENSITIVE":
        flags.append(mf)
    out.update(P2=label, P2_flags=flags, P2_note=note, mesh=mf, edges={k: v for k, v in ed.items() if k != "points"},
               points={str(k): v for k, v in ed.get("points", {}).items()}, d01_points={str(k): v for k, v in d01.items()},
               benchmarks=bench, ledger=dict(n=len(R.ledger()["integrations"]),
                                             core_h=sum(r["wall_s"] for r in R.ledger()["integrations"]) / 3600))
    import platform, scipy
    out["platform"] = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                           system=platform.platform())
    with open(os.path.join(R.OUT, "evaluation.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps({k: out[k] for k in ("P1", "P2", "P2_flags", "mesh", "persistence_bracket")}, default=float))


if __name__ == "__main__":
    main()
