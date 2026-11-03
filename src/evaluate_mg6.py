"""MG6 evaluation (proposal §4.4, §5; GATES.md MG6 pre-run specification). Reads results/MG6/runs/*.json and
checks/*.npz, re-detects every check from the sparse records (stops if counts differ from the driver), and applies the
registered outcome rules. Writes results/MG6/evaluation.json.
    python src/evaluate_mg6.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import evolve as ev
import mg6
import run_mg6 as R

T_CRIT = 2.353                                  # one-sided t, alpha = 0.05, df = 3
T_975 = 3.182                                   # two-sided 95% t interval, df = 3


def load(key):
    with np.load(os.path.join(R.OUT, "checks", key + ".npz")) as z:
        return {k: z[k] for k in z.files}


def redetect(d, theta):
    out = []
    for i in range(len(d["t"])):
        a, b = d["off"][i], d["off"][i + 1]
        rec = dict(rows=d["rows"][a:b], pos=d["pos"][a:b], s=d["s"][a:b], Lx=float(d["Lx"]), Ly=float(d["Ly"]))
        out.append(ev.detect_sparse(rec, theta))
    return out


def paired(x, y):
    """Run-level paired contrasts c = x - y (by start): values, mean, sd, t, conditional 95% interval."""
    c = np.asarray(x, float) - np.asarray(y, float)
    m, s = float(c.mean()), float(c.std(ddof=1))
    t = m / (s / 2.0) if s > 0 else (-np.inf if m < 0 else np.inf if m > 0 else 0.0)
    return dict(c=c.tolist(), mean=m, sd=s, t=float(t), ci95=[m - T_975 * s / 2, m + T_975 * s / 2],
                all_negative=bool(np.all(c < 0)), passes=bool(np.all(c < 0) and t <= -T_CRIT))


def one_sample(x):
    x = np.asarray(x, float)
    m, s = float(x.mean()), float(x.std(ddof=1))
    t = m / (s / 2.0) if s > 0 else (-np.inf if m < 0 else np.inf)
    return dict(values=x.tolist(), mean=m, sd=s, t=float(t), ci95=[m - T_975 * s / 2, m + T_975 * s / 2],
                all_negative=bool(np.all(x < 0)), passes=bool(np.all(x < 0) and t <= -T_CRIT))


def falloff(per_run_dn, per_run_null, phi_profile, edges):
    """§4.4 rule on the coarse annuli: range from r = 1.0 to the last annulus inside the first zero crossing of the
    arm's mean radial phi_G profile; >= 3 annuli with paired |t| >= 2.353, else UNRESOLVED; if resolved, exponential
    and power-law least-squares fits on run-level means weighted by run-level SEs."""
    lo, hi = np.array(edges[:-1]), np.array(edges[1:])
    mid = 0.5 * (lo + hi)
    prof = np.asarray(phi_profile, float)
    r0 = None
    for i in range(1, len(prof)):
        if prof[i - 1] < 0 <= prof[i]:
            r0 = float(mid[i - 1] + (mid[i] - mid[i - 1]) * (-prof[i - 1]) / (prof[i] - prof[i - 1]))
            break
    if prof[0] >= 0:
        r0 = 0.0
    elig = [j for j in range(len(lo)) if lo[j] >= 1.0 - 1e-12 and (r0 is not None and hi[j] <= r0 + 1e-12)]
    C = np.asarray(per_run_dn, float) - np.asarray(per_run_null, float)       # (runs, annuli)
    stats = []
    for j in elig:
        m, s = C[:, j].mean(), C[:, j].std(ddof=1)
        stats.append(dict(annulus=[float(lo[j]), float(hi[j])], mean=float(m), se=float(s / 2),
                          t=float(m / (s / 2)) if s > 0 else None))
    sig = [st for st in stats if st["t"] is not None and abs(st["t"]) >= T_CRIT]
    out = dict(first_zero_crossing=r0, eligible=stats, n_significant=len(sig))
    if len(sig) < 3:
        out["result"] = "UNRESOLVED"
        return out
    r = np.array([0.5 * sum(st["annulus"]) for st in stats])
    y = np.array([st["mean"] for st in stats])
    w = 1.0 / np.maximum(np.array([st["se"] for st in stats]), 1e-12)
    from scipy.optimize import curve_fit
    fits = {}
    for name, fn, p0 in (("exponential", lambda x, A, lam: A * np.exp(-x / lam), (y[0], 1.0)),
                         ("power", lambda x, A, p: A * x ** (-p), (y[0], 1.0))):
        try:
            p, cov = curve_fit(fn, r, y, p0=p0, sigma=1 / w, absolute_sigma=True, maxfev=20000)
            fits[name] = dict(params=p.tolist(), se=np.sqrt(np.diag(cov)).tolist(),
                              resid=(y - fn(r, *p)).tolist())
        except Exception as e:
            fits[name] = dict(error=str(e))
    out.update(result="RESOLVED", fits=fits)
    return out


def main():
    dec = R.timing_params()
    w0, t_end = dec["window"]
    table, mismatches = {}, []
    for s in R.run_specs():
        path = os.path.join(R.OUT, "runs", s["key"] + ".json")
        if not os.path.exists(path):
            table[s["key"]] = dict(spec=s, status=R.INCOMPLETE)
            continue
        rec = json.load(open(path))
        row = dict(spec=s, status=rec["status"], note=rec.get("note"), flags=len(rec.get("flags", [])))
        if rec["status"] == R.COMPLETE:
            d = load(s["key"])
            lists = redetect(d, rec["theta"])
            if [len(x) for x in lists] != d["n_parts"].tolist():
                mismatches.append(s["key"])
            cls = ev.classify_run(lists[1:], lists[0][0]["centre"], rec["Lx"], rec["Ly"]) if len(lists[0]) == 1 \
                else dict(cls="unresolved", sub="start not single")
            if cls.get("cls") == "single":
                cls["D_p"] = ev.msd_dp(cls["track"])
            if cls.get("cls") != rec["summary"]["cls"].get("cls"):
                mismatches.append(s["key"])
            t = d["t"]
            win = (t >= w0 - 1e-9) & (t <= t_end + 1e-9)
            single_win = bool(np.all(d["n_parts"][win] == 1))
            row.update(cls=cls, single_in_window=single_win, n_window=int(win.sum()))
            if single_win:
                row["Y_D0"] = float(d["m_D0"][win, 1].mean())
                row["Y_Dc"] = float(d["m_Dc"][win, 1].mean())
                row["Y_P0"] = float(d["m_P0"][win, 1].mean())
                row["Y_Pc"] = float(d["m_Pc"][win, 1].mean())
                row["ann_dn0"] = d["m_ann_dn0"][win].mean(0).tolist()
                row["ann_phiG"] = d["m_ann_phiG"][win].mean(0).tolist()
                row["fine_dn0"] = d["m_fine_dn0"][win].mean(0).tolist()
                row["fine_phiG"] = d["m_fine_phiG"][win].mean(0).tolist()
                row["fine_P0"] = d["m_fine_P0"][win].mean(0).tolist()
                row["obs_over_pred"] = (row["Y_D0"] / row["Y_P0"]) if row["Y_P0"] != 0 else None
            row["M_minus"] = d["m_M_minus"].tolist()
            row["depth"] = [min(p["depth"] for p in x) if x else None for x in lists]
        table[s["key"]] = row
    if mismatches:
        raise SystemExit(f"re-detection from sparse records differs from the driver for: {sorted(set(mismatches))}")
    out = score(table)
    out["runs"] = table
    out["window"] = [w0, t_end]
    ch, done = R.spent()
    out["core_h_ledger"] = ch
    json.dump(out, open(os.path.join(R.OUT, "evaluation.json"), "w"), indent=1, default=float)
    print("Gate MG6:", out["gate"])
    print("Feedback:", {k: v["outcome"] for k, v in out["feedback"].items()})
    print("Absolute deficit (CP2):", out.get("absolute_cp2", {}).get("passes"))
    print("Falloff:", out.get("falloff", {}).get("result"), "| CP4:", out.get("falloff_cp4", {}).get("result"))
    print(f"core-h (ledger) {ch:.2f}")


def score(table):
    """Pure function of the per-run table (proposal §5)."""
    arm = lambda a: [table[f"{a}_D{sl}_p{pl}"] for sl, pl in R.DRUNS]
    out = dict(statuses={k: v["status"] for k, v in table.items()})
    if any(v["status"] == R.INVALID_BOUND for v in table.values()):
        out["campaign"] = "STOPPED (INVALID-BOUND): PI decides"
    feedback = {}
    for a in ("CP2", "CP4"):
        rs = arm(a)
        if any(r["status"] != R.COMPLETE for r in rs):
            feedback[a] = dict(outcome="INCOMPLETE")
            continue
        cl = [r["cls"]["cls"] for r in rs]
        feedback[a] = dict(outcome="SURVIVES" if sum(c == "single" for c in cl) >= 3 else "CHANGED", classes=cl)
    out["feedback"] = feedback
    prim = arm("CP2") + arm("SH2") + arm("DEC")
    if any(r["status"] != R.COMPLETE for r in prim):
        out["gate"] = "INCOMPLETE"
        return out
    if not all(r["single_in_window"] for r in prim):
        out["gate"] = "NOT SCORED (feedback)"
        return out
    Y = lambda a, k="Y_D0": [r[k] for r in arm(a)]
    vs_sh, vs_dec = paired(Y("CP2"), Y("SH2")), paired(Y("CP2"), Y("DEC"))
    out["primary"] = dict(CP2_minus_SH2=vs_sh, CP2_minus_DEC=vs_dec)
    out["gate"] = "DEMONSTRATED" if vs_sh["passes"] and vs_dec["passes"] else "NOT DEMONSTRATED"
    out["absolute_cp2"] = one_sample(Y("CP2"))
    desc = {}
    if all(r["status"] == R.COMPLETE and r.get("single_in_window") for r in arm("CP4") + arm("SH4")):
        desc["CP4_minus_SH4"] = paired(Y("CP4"), Y("SH4"))
        desc["CP4_minus_DEC"] = paired(Y("CP4"), Y("DEC"))
        desc["Dc_CP4_minus_SH4"] = paired(Y("CP4", "Y_Dc"), Y("SH4", "Y_Dc"))
        desc["Dc_CP4_minus_DEC"] = paired(Y("CP4", "Y_Dc"), Y("DEC", "Y_Dc"))            # MG022 review item 1
    desc["Dc_CP2_minus_SH2"] = paired(Y("CP2", "Y_Dc"), Y("SH2", "Y_Dc"))
    desc["Dc_CP2_minus_DEC"] = paired(Y("CP2", "Y_Dc"), Y("DEC", "Y_Dc"))
    desc["SH2_minus_DEC"] = paired(Y("SH2"), Y("DEC"))
    if all(r["status"] == R.COMPLETE and r.get("single_in_window") for r in arm("SH4")):
        desc["SH4_minus_DEC"] = paired(Y("SH4"), Y("DEC"))
    desc["obs_over_pred_CP2"] = [r["obs_over_pred"] for r in arm("CP2")]
    out["descriptive"] = desc
    prof = np.mean([r["ann_phiG"] for r in arm("CP2")], axis=0)
    out["falloff"] = falloff([r["ann_dn0"] for r in arm("CP2")], [r["ann_dn0"] for r in arm("SH2")], prof, mg6.ANN)
    if all(r["status"] == R.COMPLETE and r.get("single_in_window") for r in arm("CP4") + arm("SH4")):
        prof4 = np.mean([r["ann_phiG"] for r in arm("CP4")], axis=0)                 # proposal 4.4: kappa = 4 too
        out["falloff_cp4"] = falloff([r["ann_dn0"] for r in arm("CP4")], [r["ann_dn0"] for r in arm("SH4")], prof4,
                                     mg6.ANN)
    out["fine_profile"] = dict(edges=list(mg6.FINE),
                               contrast=np.mean(np.subtract([r["fine_dn0"] for r in arm("CP2")],
                                                            [r["fine_dn0"] for r in arm("SH2")]), 0).tolist(),
                               phiG=np.mean([r["fine_phiG"] for r in arm("CP2")], 0).tolist(),
                               P0=np.mean([r["fine_P0"] for r in arm("CP2")], 0).tolist())
    return out


if __name__ == "__main__":
    main()
