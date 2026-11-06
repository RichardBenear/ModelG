"""Independent recompute of MG8 (proposal draft 4 §7). Reads only the committed records (results/MG8/runs, fields, snap,
status.json, part0/decisions.json) and re-derives every label, interval and fit with its own code: its own Fourier shift
and Nyquist filter (rfft), its own alignment search (scipy bounded minimisation plus a dense scan), its own fate, window,
settling, label and refinement rules written from the proposal text, and scipy's own exact binomial and Welch intervals.
It shares no measurement code with mg8.py / evaluate_mg8.py. Compares with evaluation.json and writes
results/MG8/recompute.json. Usage: python src/recompute_mg8.py [outdir]"""
import os, sys, json
import numpy as np
from scipy import stats, optimize

ROOT = os.path.join(os.path.dirname(__file__), "..")
H, N, LD, TAU0 = 0.2, 500, 100.0, 6.667e-4


# ------------------------------------------------------------------ alignment (own implementation)
def rshift(f, s, filt):
    """f(x + s) by band-limited interpolation; filt drops the Nyquist coefficient."""
    F = np.fft.rfft(f, axis=-1)
    k = 2 * np.pi * np.fft.rfftfreq(N, d=H)
    ph = np.exp(1j * k * s)
    if filt:
        F[..., -1] = 0.0
    else:
        ph[-1] = np.cos(k[-1] * s)                 # the real part is all a real-valued shift keeps of that mode
    return np.fft.irfft(F * ph, n=N, axis=-1)


def nyq(f):
    a = (-1.0) ** np.arange(N)
    return np.outer((f * a).sum(-1) / N, a)


def E_shape(A, B, filt=True):
    A, B = np.asarray(A, float), np.asarray(B, float)
    if filt:
        A, B = A - nyq(A), B - nyq(B)
    nB = np.sqrt((B ** 2).sum(-1))
    obj = lambda s: float(np.max(np.sqrt(((rshift(A, s, False) - B) ** 2).sum(-1)) / nB))
    ints = [(obj(m * H), abs(m * H), m * H) for m in range(-N // 2, N // 2)]
    e0, _, s0 = min(ints)
    r = optimize.minimize_scalar(obj, bounds=(s0 - H / 2, s0 + H / 2), method="bounded", options=dict(xatol=1e-12))
    scan = min((obj(s), s) for s in np.linspace(s0 - H / 2, s0 + H / 2, 101))
    best = min([(r.fun, r.x), scan, (e0, s0)])
    if best[1] != r.x:                             # refine around a scan point that beat the minimiser
        r2 = optimize.minimize_scalar(obj, bounds=(best[1] - H / 100, best[1] + H / 100), method="bounded",
                                      options=dict(xatol=1e-12))
        best = min(best, (r2.fun, r2.x))
    return best[0], best[1]


# ------------------------------------------------------------------ B1
def load(out):
    runs = {}
    for f in sorted(os.listdir(os.path.join(out, "runs"))):
        if f.endswith(".json"):
            runs[f[:-5]] = json.load(open(os.path.join(out, "runs", f)))
    return runs


def flds(out, key):
    z = np.load(os.path.join(out, "fields", key + ".npz"))
    return {int(k[1:]): z[k] for k in z.files}


def unwrap(c):
    c = np.asarray(c, float)
    return c[0] + np.concatenate([[0.0], np.cumsum((np.diff(c) + LD / 2) % LD - LD / 2)])


def vel(rec, u0, u1):
    u = np.asarray(rec["unit"]); t = np.asarray(rec["t"]); c = unwrap(rec["centre"])
    w = (u >= u0) & (u <= u1)
    assert w.sum() == u1 - u0 + 1
    return float(np.polyfit(t[w], c[w], 1)[0])


def single(rec, U):
    u = np.asarray(rec["unit"]); n = np.asarray(rec["n_dets"]); f = np.asarray(rec["finite_units"])
    if not set(range(1, U + 1)) <= set(u.tolist()):
        return False
    return bool(np.all(n[(u >= 50) & (u <= U)] == 1) and np.all(f[u <= U]))


def is_settled(E, E1, v, v1):
    return abs(E - E1) <= 0.01 * E + 1e-7 and abs(v - v1) <= 0.1 * abs(v) + 1e-7


def labels(sg, E, Ef, v, vf, st):
    if not sg:
        return "CHANGED", None, None, "SENSITIVE"
    shape = ("INSUFFICIENT RESOLUTION" if 10 * Ef >= 1e-3 else "UNSETTLED" if not st else
             "MATERIAL" if E >= 1e-2 else "SMALL" if E >= 1e-3 else "NEGLIGIBLE")
    drift = ("INSUFFICIENT RESOLUTION" if 10 * vf >= 1e-6 else "UNSETTLED" if not st else
             "DRIFTING" if abs(v) >= 1e-6 else "NOT DRIFTING AT THIS HORIZON")
    if shape == "MATERIAL" or drift == "DRIFTING":
        s = "SENSITIVE"
    elif shape == "NEGLIGIBLE" and drift == "NOT DRIFTING AT THIS HORIZON":
        s = "PRACTICALLY INSENSITIVE"
    elif shape == "SMALL" and drift == "NOT DRIFTING AT THIS HORIZON":
        s = "MINOR"
    elif "INSUFFICIENT RESOLUTION" in (shape, drift):
        s = "INSUFFICIENT RESOLUTION"
    else:
        s = "UNRESOLVED"
    return "SAME", shape, drift, s


def b1(out, runs):
    lev = {}
    for k in range(7):
        Jr = runs[f"det_J_k{k}_r0"]; Jf = flds(out, f"det_J_k{k}_r0")
        Ef = E_shape(Jf[250], Jf[300])[0] + 1e-9
        vf = abs(vel(Jr, 150, 300)) + 1e-8
        reps = []
        for j in range(3):
            rec = runs[f"det_R_k{k}_r{j}"]; Rf = flds(out, f"det_R_k{k}_r{j}")
            sg = single(rec, 300)
            ok = all(u in Rf and np.all(np.isfinite(Rf[u])) for u in (250, 300))
            E, sig = E_shape(Rf[300], Jf[300]) if ok else (np.nan, np.nan)
            E1 = E_shape(Rf[250], Jf[250])[0] if ok else np.nan
            Eraw = E_shape(Rf[300], Jf[300], filt=False)[0] if ok else np.nan
            v, v1 = vel(rec, 150, 300), vel(rec, 100, 250)
            st = bool(ok and is_settled(E, E1, v, v1))
            reps.append(dict(single=sg, E=E, E1=E1, E_raw=Eraw, sigma=sig, v=v, v1=v1, settled=st,
                             labels=labels(sg, E, Ef, v, vf, st)))
        Es = [r["E"] for r in reps]
        Ek = float(np.median(Es))
        sm = [r["labels"][3] for r in reps]
        lev[k] = dict(Efloor=Ef, vfloor=vf, reps=reps, Ek=Ek, S=(max(Es) - min(Es)) / max(Ek, 10 * Ef),
                      all_single=all(r["single"] for r in reps), all_settled=all(r["settled"] for r in reps),
                      summary=sm[0] if len(set(sm)) == 1 else "SEAM-DEPENDENT")
    rec = runs["det_R_k6_r0"]; Rf = flds(out, "det_R_k6_r0"); Jf = flds(out, "det_J_k6_r0")
    E600, E550 = E_shape(Rf[600], Jf[600])[0], E_shape(Rf[550], Jf[550])[0]
    v, v1 = vel(rec, 450, 600), vel(rec, 400, 550)
    hz = dict(single=single(rec, 600), settled=bool(is_settled(E600, E550, v, v1)), E300=lev[6]["reps"][0]["E"],
              E600=E600, Efloor=lev[6]["Efloor"])
    hz_rel = hz["single"] and hz["settled"] and abs(E600 - hz["E300"]) <= 0.05 * hz["E300"] + 1e-7
    elig = [k for k in range(2, 7) if lev[k]["all_single"] and lev[k]["all_settled"] and lev[k]["S"] <= 0.25
            and lev[k]["Ek"] >= 10 * lev[k]["Efloor"]]
    ref = "NOT RESOLVED BY THIS LADDER"
    info = dict(eligible=elig)
    if len(elig) >= 4 and TAU0 / 2 ** min(elig) / (TAU0 / 2 ** max(elig)) >= 8 - 1e-9:
        x = np.log([TAU0 / 2 ** k for k in elig]); y = np.log([lev[k]["Ek"] for k in elig])
        X = np.vstack([x, np.ones_like(x)]).T
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ coef
        se = np.sqrt(res @ res / (len(x) - 2) / ((x - x.mean()) ** 2).sum())
        lo = coef[0] - stats.t.ppf(0.975, len(x) - 2) * se
        info.update(p=float(coef[0]), p_lo=float(lo), resid=res.tolist())
        if coef[0] >= 0.7 and lo >= 0.5 and np.all(np.abs(res) <= 0.15) and \
                lev[max(elig)]["Ek"] <= 0.1 * lev[min(elig)]["Ek"] and hz_rel:
            ref = "CONVERGES TO J"
    if ref != "CONVERGES TO J" and len(elig) < 4 and any(lev[k]["Ek"] < 10 * lev[k]["Efloor"] for k in range(2, 7)):
        hzb = hz["single"] and hz["settled"] and (hz_rel if hz["E300"] >= 10 * hz["Efloor"] else E600 < 1e-4)
        if all(10 * lev[k]["Efloor"] < 1e-4 for k in range(2, 7)) and \
                all(lev[k]["all_single"] and lev[k]["all_settled"] for k in range(7)) and \
                all(lev[k]["S"] <= 0.25 for k in range(7) if lev[k]["Ek"] >= 10 * lev[k]["Efloor"]) and \
                all(lev[k + 1]["Ek"] <= lev[k]["Ek"] for k in range(6)) and lev[6]["Ek"] < 1e-4 and hzb:
            ref = "BELOW MEASUREMENT RESOLUTION"
    mirror = {}
    for k, TL in ((0, 300), (6, 50)):
        Lf, Rf = flds(out, f"det_L_k{k}_r0"), flds(out, f"det_R_k{k}_r{0}")
        mirror[k] = max(float(np.abs(np.flip(Lf[u], -1) - Rf[u]).max() / np.abs(Rf[u]).max())
                        for u in range(50, TL + 1, 50))
    return dict(levels=lev, horizon=hz, production=lev[0]["summary"], refinement=ref, refinement_info=info,
                mirror=mirror)


# ------------------------------------------------------------------ B2
def b2(out, runs):
    res = {}
    prof = {}
    for a in ("J", "R", "P"):
        rs = [runs[f"sto_{a}_k0_r{i}"] for i in range(8)]
        ev = {t: sum(any(e[0] == t for e in r["events"]) for r in rs) for t in ("loss", "replication", "merge")}
        el = []
        for i, r in enumerate(rs):
            phys = any(e[0] in ("loss", "replication", "merge") for e in r["events"])
            assert r["eligible"] == ((not phys) and not r["track_fail"])
            if r["eligible"]:
                z = np.load(os.path.join(out, "snap", f"sto_{a}_k0_r{i}.npz"))
                nb = z["n_per_block"].astype(float); bl = np.nan_to_num(z["blocks"])
                p = np.tensordot(nb, bl, axes=(0, 0)) / nb.sum()
                assert np.allclose(p, np.asarray(r["profile"]), rtol=0, atol=1e-12)
                el.append((r, p))
        ci = {t: list(stats.binomtest(ev[t], 8).proportion_ci(0.95, method="exact")) for t in ev}
        res[a] = dict(runs_with_event=ev, event_ci=ci, n_eligible=len(el), v=[r["v"] for r, _ in el])
        prof[a] = [p for _, p in el]
    for a in ("R", "P"):
        c = {}
        if res[a]["n_eligible"] >= 5 and res["J"]["n_eligible"] >= 5:
            tt = stats.ttest_ind(res[a]["v"], res["J"]["v"], equal_var=False).confidence_interval(0.95)
            c["drift_diff"] = [float(np.mean(res[a]["v"]) - np.mean(res["J"]["v"])), float(tt.low), float(tt.high)]
            c["E_stoch"] = E_shape(np.mean(prof[a], 0), np.mean(prof["J"], 0))[0]
        res[f"{a}_vs_J"] = c
    pj = prof["J"]
    h1 = (len(pj) + 1) // 2
    res["J_split_half"] = E_shape(np.mean(pj[:h1], 0), np.mean(pj[h1:], 0))[0] if len(pj) >= 2 else None
    res["mirror_stoch"] = runs["sto_L_k0_r0"]["counts_sha256"] == runs["sto_R_k0_r0"]["counts_sha256"]
    return res


def close(a, b, rel=1e-6, ab=1e-12):
    return a is None and b is None or (a is not None and b is not None and abs(a - b) <= rel * abs(b) + ab)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "MG8")
    runs = load(out)
    rc = dict(B1=b1(out, runs), B2=b2(out, runs))
    ev = json.load(open(os.path.join(out, "evaluation.json")))
    cmp = []
    e1 = ev["B1"]
    for k in range(7):
        lv, ln = rc["B1"]["levels"][k], e1["levels"][str(k)]
        cmp.append((f"k{k} summary", lv["summary"] == ln["summary"]))
        cmp.append((f"k{k} Ek", close(lv["Ek"], ln["Ek"])))
        cmp.append((f"k{k} Efloor", close(lv["Efloor"], ln["Efloor"])))
        for j in range(3):
            a, b = lv["reps"][j], ln["reps"][j]
            cmp.append((f"k{k} r{j} single/settled", a["single"] == b["single"] and a["settled"] == b["settled"]))
            cmp.append((f"k{k} r{j} v", close(a["v"], b["v"], 1e-6, 1e-12)))
            cmp.append((f"k{k} r{j} E_raw", close(a["E_raw"], b["E_raw"], 1e-6, 1e-12)))
    cmp.append(("production", rc["B1"]["production"] == e1["production"]))
    cmp.append(("refinement", rc["B1"]["refinement"] == e1["refinement"]))
    cmp.append(("horizon E600", close(rc["B1"]["horizon"]["E600"], e1["horizon"]["E600"])))
    for k in (0, 6):
        cmp.append((f"mirror k{k}", (rc["B1"]["mirror"][k] <= 1e-12) == (e1["mirror"][str(k)]["status"] == "PASS")))
    e2 = ev["B2"]
    for a in ("J", "R", "P"):
        cmp.append((f"{a} events", rc["B2"][a]["runs_with_event"] == e2["arms"][a]["runs_with_event"]))
        cmp.append((f"{a} n_eligible", rc["B2"][a]["n_eligible"] == e2["arms"][a]["n_eligible"]))
        for t in ("loss", "replication", "merge"):
            cmp.append((f"{a} {t} CI", np.allclose(rc["B2"][a]["event_ci"][t], e2["arms"][a]["event_ci"][t], atol=1e-9)))
    for a in ("R", "P"):
        c, d = rc["B2"][f"{a}_vs_J"], e2.get(f"{a}_vs_J", {})
        if "drift_diff" in c:
            dd = d["drift_diff"]
            cmp.append((f"{a} drift diff", np.allclose(c["drift_diff"], [dd["d"], dd["lo"], dd["hi"]], rtol=1e-6, atol=1e-15)))
            cmp.append((f"{a} E_stoch", close(c["E_stoch"], d["E_stoch"]["E"])))
    cmp.append(("J split half", close(rc["B2"]["J_split_half"], e2["J_split_half"])))
    cmp.append(("mirror stoch", rc["B2"]["mirror_stoch"] == (e2["mirror_stoch"]["status"] == "PASS")))
    rc["comparison"] = dict(n=len(cmp), mismatches=[n for n, ok in cmp if not ok])
    json.dump(rc, open(os.path.join(out, "recompute.json"), "w"), indent=1, default=float)
    print(json.dumps(dict(production=rc["B1"]["production"], refinement=rc["B1"]["refinement"],
                          comparison=rc["comparison"]), indent=1))


if __name__ == "__main__":
    main()
