"""MG7 evaluation (proposal §4, §5; GATES.md MG7 pre-run specification). Reads results/MG7/runs/*.json and
checks/*.npz, plus MG6's recorded b = 28 CP2 and SH2 runs (midpoint, reported only), re-detects every check from the
sparse records (stops if counts differ from the driver), and applies the registered outcome rules.
Writes results/MG7/evaluation.json.
    python src/evaluate_mg7.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import evolve as ev
import mg6
import evaluate_mg6 as E6
import run_mg7 as R

T_CRIT, T_975 = E6.T_CRIT, E6.T_975
MG6_OUT = os.path.join(R.ROOT, "MG6")
WIN = (100.0, 200.0)
EARLY, LATE = (100.0, 140.0), (160.0, 200.0)
STARTS = [f"D{s}_p{p}" for s, p in R.DRUNS]


def load(path):
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def _mask(t, lo, hi):
    return (t >= lo - 1e-9) & (t <= hi + 1e-9)


def run_row(rec, d, theta_alt=None):
    """Per-run table row from the run record and its check record (all quantities at r = 1.5 unless named)."""
    lists = E6.redetect(d, rec["theta"])
    mismatch = [len(x) for x in lists] != d["n_parts"].tolist()
    cls = ev.classify_run(lists[1:], lists[0][0]["centre"], rec["Lx"], rec["Ly"]) if len(lists[0]) == 1 \
        else dict(cls="unresolved", sub="start not single")
    if cls.get("cls") == "single":
        cls["D_p"] = ev.msd_dp(cls["track"])
    if "summary" in rec and cls.get("cls") != rec["summary"]["cls"].get("cls"):
        mismatch = True
    t = d["t"]
    win = _mask(t, *WIN)
    single = bool(np.all(d["n_parts"][win] == 1))
    row = dict(status=rec["status"], cls=dict(cls=cls.get("cls"), sub=cls.get("sub"), mobility=cls.get("mobility"),
                                              D_p=cls.get("D_p")),
               single_in_window=single, n_window=int(win.sum()), P_start=float(d["m_P0"][0, 1]),
               M_minus_start=float(d["m_M_minus"][0]))
    if single:
        row.update(Y_D0=float(d["m_D0"][win, 1].mean()), Y_Dc=float(d["m_Dc"][win, 1].mean()),
                   Y_P0=float(d["m_P0"][win, 1].mean()), M_minus=float(d["m_M_minus"][win].mean()),
                   early=float(d["m_D0"][_mask(t, *EARLY), 1].mean()), late=float(d["m_D0"][_mask(t, *LATE), 1].mean()),
                   ann_dn0=d["m_ann_dn0"][win].mean(0).tolist(), ann_phiG=d["m_ann_phiG"][win].mean(0).tolist(),
                   fine_dn0=d["m_fine_dn0"][win].mean(0).tolist(), fine_phiG=d["m_fine_phiG"][win].mean(0).tolist())
        row["early_late_absdiff"] = abs(row["late"] - row["early"])
    if theta_alt is not None:
        alt = E6.redetect(d, theta_alt)
        row["alt_threshold"] = dict(theta=theta_alt, counts=[len(x) for x in alt],
                                    single_in_window=bool(all(len(x) == 1 for x, w in zip(alt, win) if w)))
    return row, mismatch


def fieller(a, b, tq=T_975):
    """Paired Fieller 95% interval for mean(a) / mean(b) over the n starts (conditional working model, df = n - 1).
    Unbounded when the denominator is indistinguishable from zero at that level."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = len(a)
    ma, mb = a.mean(), b.mean()
    C = np.cov(a, b, ddof=1) / n
    A = mb ** 2 - tq ** 2 * C[1, 1]
    B = ma * mb - tq ** 2 * C[0, 1]
    Cq = ma ** 2 - tq ** 2 * C[0, 0]
    out = dict(ratio=float(ma / mb) if mb != 0 else None, sum_ratio=float(a.sum() / b.sum()) if b.sum() != 0 else None)
    disc = B ** 2 - A * Cq
    if A <= 0 or disc < 0:
        out.update(ci95=None, bounded=False)
    else:
        r = np.sqrt(disc)
        out.update(ci95=[float((B - r) / A), float((B + r) / A)], bounded=True)
    return out


def main():
    dec = R.part0_decision()
    table, mismatches = {}, []
    for s in R.run_specs():
        path = os.path.join(R.OUT, "runs", s["key"] + ".json")
        if not os.path.exists(path):
            table[s["key"]] = dict(spec=s, status=R.INCOMPLETE)
            continue
        rec = json.load(open(path))
        if rec["b"] != s["b"] or rec.get("P_b") != s["b"] or rec["start_field"]["shift_applied"]:
            raise SystemExit(f"{s['key']}: run record b or start field does not match the registered spec")
        row = dict(spec=s, status=rec["status"], note=rec.get("note"), flags=len(rec.get("flags", [])),
                   max_phiG=rec.get("max_phiG"))
        if rec["status"] == R.COMPLETE:
            r, bad = run_row(rec, load(os.path.join(R.OUT, "checks", s["key"] + ".npz")), R.THETA_B[s["b"]])
            row.update(r)
            if bad:
                mismatches.append(s["key"])
        table[s["key"]] = row
    mid = {}
    for arm in ("CP2", "SH2"):
        for st in STARTS:
            rec = json.load(open(os.path.join(MG6_OUT, "runs", f"{arm}_{st}.json")))
            r, bad = run_row(rec, load(os.path.join(MG6_OUT, "checks", f"{arm}_{st}.npz")))
            mid[f"{arm}_b28.00_{st}"] = r
            if bad:
                mismatches.append(f"MG6 {arm}_{st}")
    if mismatches:
        raise SystemExit(f"re-detection from sparse records differs from the driver for: {sorted(set(mismatches))}")
    out = score(table, mid)
    out["runs"] = table
    out["midpoint_runs_mg6"] = mid
    out["window"] = list(WIN)
    out["part0"] = dict(proceed=dec.get("proceed"), stop=dec.get("stop"))
    ch, done = R.spent()
    out["core_h_ledger"] = ch
    json.dump(out, open(os.path.join(R.OUT, "evaluation.json"), "w"), indent=1, default=float)
    print("Gate MG7:", out["gate"])
    if "primary" in out:
        p = out["primary"]
        print(f"  Delta_i = {np.round(p['values'], 2).tolist()}, mean {p['mean']:.2f}, t {p['t']:.3f}")
    print("Feedback:", {k: v["outcome"] for k, v in out["feedback"].items()})
    print(f"core-h (ledger) {ch:.2f}")


def score(table, mid):
    """Pure function of the per-run table and the MG6 midpoint rows (proposal §5)."""
    row = lambda arm, b, st: (mid if b == 28.0 else table)[f"{arm}_b{b:.2f}_{st}"]
    arm = lambda a, b: [row(a, b, st) for st in STARTS]
    out = dict(statuses={k: v["status"] for k, v in table.items()})
    if any(v["status"] == R.INVALID_BOUND for v in table.values()):
        out["campaign"] = "STOPPED (INVALID-BOUND): PI decides"
    feedback = {}
    for b in R.B_END:
        rs = arm("CP2", b)
        if any(r["status"] != R.COMPLETE for r in rs):
            feedback[f"b{b:.2f}"] = dict(outcome="INCOMPLETE")
            continue
        cl = [r["cls"]["cls"] for r in rs]
        feedback[f"b{b:.2f}"] = dict(outcome="SURVIVES" if sum(c == "single" for c in cl) >= 3 else "CHANGED",
                                    classes=cl, mobility=[r["cls"].get("mobility") for r in rs])
    out["feedback"] = feedback
    allr = [r for b in R.B_END for a in ("CP2", "SH2") for r in arm(a, b)]
    if any(r["status"] != R.COMPLETE for r in allr):                  # INCOMPLETE precedes every other rule
        out["gate"] = "INCOMPLETE"
        return out
    if not all(r["single_in_window"] for r in allr):
        out["gate"] = "NOT SCORED (feedback)"
        return out
    c = {b: np.array([r1["Y_D0"] - r2["Y_D0"] for r1, r2 in zip(arm("CP2", b), arm("SH2", b))]) for b in R.B_END}
    D = c[28.5] - c[27.5]
    prim = E6.one_sample(D)
    out["primary"] = prim
    out["gate"] = "LARGER AT HIGHER b" if prim["passes"] else "NOT DEMONSTRATED"
    rep = {}
    for b in (27.5, 28.0, 28.5):
        if all(r.get("single_in_window") for a in ("CP2", "SH2") for r in arm(a, b)):
            rep[f"b{b:.2f}"] = E6.paired([r["Y_D0"] for r in arm("CP2", b)], [r["Y_D0"] for r in arm("SH2", b)])
    out["contrasts"] = rep
    out["trajectories"] = {st: {f"b{b:.2f}": rep[f"b{b:.2f}"]["c"][i] for b in (27.5, 28.0, 28.5)
                                if f"b{b:.2f}" in rep} for i, st in enumerate(STARTS)}
    resid = {}
    for name, key in (("start_field", "P_start"), ("in_run", "Y_P0")):
        pc = {b: np.array([r1[key] - r2[key] for r1, r2 in zip(arm("CP2", b), arm("SH2", b))]) for b in R.B_END}
        DP = pc[28.5] - pc[27.5]
        e = D - DP
        m, s = float(e.mean()), float(e.std(ddof=1))
        resid[name] = dict(DeltaP=DP.tolist(), e=e.tolist(), mean=m, sd=s,
                           ci95=[m - T_975 * s / 2, m + T_975 * s / 2], fieller_ratio=fieller(D, DP))
    out["residuals"] = resid
    xs, ys = [], []
    for b in (27.5, 28.0, 28.5):
        if f"b{b:.2f}" in rep:
            xs += [r["M_minus"] for r in arm("CP2", b)]
            ys += rep[f"b{b:.2f}"]["c"]
    if len(xs) >= 3:
        A = np.column_stack([np.ones(len(xs)), xs])
        coef, *_ = np.linalg.lstsq(A, np.asarray(ys), rcond=None)
        out["slope_vs_M_minus"] = dict(intercept=float(coef[0]), slope=float(coef[1]), n=len(xs), M_minus=xs, c=ys,
                                       note="descriptive; read beside the four within-start trajectories")
    out["falloff"] = {}
    for b in R.B_END:
        prof = np.mean([r["ann_phiG"] for r in arm("CP2", b)], axis=0)
        out["falloff"][f"b{b:.2f}"] = E6.falloff([r["ann_dn0"] for r in arm("CP2", b)],
                                                 [r["ann_dn0"] for r in arm("SH2", b)], prof, mg6.ANN)
    out["threshold_sensitivity"] = {k: v["alt_threshold"] for k, v in table.items() if "alt_threshold" in v}
    out["early_late"] = {k: dict(early=v["early"], late=v["late"], absdiff=v["early_late_absdiff"],
                                 P_start_ref=abs(v["P_start"])) for k, v in table.items() if "early" in v}
    out["start_field_P0"] = {k: v["P_start"] for k, v in table.items() if "P_start" in v}
    return out


if __name__ == "__main__":
    main()
