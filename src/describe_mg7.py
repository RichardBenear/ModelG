"""MG7 descriptive summary (post hoc, NOT registered; written after the data were seen). The registered outcome is
NOT SCORED (feedback) (results/MG7/evaluation.json), so the evaluator stops before its descriptive section. This script
reports, from the committed check records only:
  - per run: particle counts, first check with no particle, window mean D0(1.5) where the run stays single;
  - the CP2 - SH2 contrast at b = 28.5 (MG7) and b = 28.0 (MG6 records), per start. The two use DIFFERENT candidate
    streams (MG7 shares one stream across b; MG6 used its own), so their difference is not the registered paired test
    and is not a substitute for it.
Writes results/MG7/descriptive.json.
    python src/describe_mg7.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))

ROOT = os.path.join(os.path.dirname(__file__), "..", "results")
STARTS = ["D0_p0", "D0_p1", "D1_p0", "D1_p1"]


def run(path):
    with np.load(path) as z:
        t, n, D0, P0, Mm = z["t"], z["n_parts"], z["m_D0"][:, 1], z["m_P0"][:, 1], z["m_M_minus"]
    w = (t >= 100 - 1e-9) & (t <= 200 + 1e-9)
    zero = t[n == 0]
    out = dict(n_parts=n.tolist(), first_t_no_particle=float(zero[0]) if len(zero) else None,
               single_in_window=bool(np.all(n[w] == 1)), P0_start=float(P0[0]))
    if out["single_in_window"]:
        out.update(Y_D0=float(D0[w].mean()), M_minus=float(Mm[w].mean()))
    return out


def contrast(rows):
    c = np.array([rows["CP2"][s]["Y_D0"] - rows["SH2"][s]["Y_D0"] for s in STARTS])
    m, sd = float(c.mean()), float(c.std(ddof=1))
    return dict(c=c.tolist(), mean=m, sd=sd, t=m / (sd / 2), all_negative=bool(np.all(c < 0)))


def main():
    out = dict(note="post hoc descriptive summary; not registered; the registered MG7 outcome is NOT SCORED (feedback)")
    runs = {}
    for b in ("27.50", "28.50"):
        for a in ("SH2", "CP2"):
            for s in STARTS:
                runs[f"{a}_b{b}_{s}"] = run(os.path.join(ROOT, "MG7", "checks", f"{a}_b{b}_{s}.npz"))
    out["runs"] = runs
    r285 = {a: {s: runs[f"{a}_b28.50_{s}"] for s in STARTS} for a in ("CP2", "SH2")}
    r280 = {a: {s: run(os.path.join(ROOT, "MG6", "checks", f"{a}_{s}.npz")) for s in STARTS} for a in ("CP2", "SH2")}
    out["contrast_b28.50"] = contrast(r285)
    out["contrast_b28.00_mg6"] = contrast(r280)
    out["difference_28.50_minus_28.00_unpaired_streams"] = (np.array(out["contrast_b28.50"]["c"])
                                                            - np.array(out["contrast_b28.00_mg6"]["c"])).tolist()
    out["M_minus_CP2"] = {"b28.50": [r285["CP2"][s]["M_minus"] for s in STARTS],
                          "b28.00": [r280["CP2"][s]["M_minus"] for s in STARTS]}
    ev = json.load(open(os.path.join(ROOT, "MG7", "evaluation.json")))   # MG025 review: carry what the evaluator stored
    out["threshold_sensitivity_single_in_window"] = {k: v["alt_threshold"]["single_in_window"]
                                                     for k, v in ev["runs"].items() if "alt_threshold" in v}
    el = {}
    for s in STARTS:
        c, h = ev["runs"][f"CP2_b28.50_{s}"], ev["runs"][f"SH2_b28.50_{s}"]
        el[s] = dict(CP2_late_minus_early=c["late"] - c["early"], SH2_late_minus_early=h["late"] - h["early"],
                     contrast_late_minus_early=(c["late"] - h["late"]) - (c["early"] - h["early"]))
    out["early_late_b28.50"] = el
    out["unavailable"] = ["contrast at b = 27.5", "endpoint change Delta and its t", "residuals (start-field, in-run)",
                          "Fieller ratio", "twelve-run slope against M-", "falloff at b = 27.5"]
    json.dump(out, open(os.path.join(ROOT, "MG7", "descriptive.json"), "w"), indent=1)
    for k, v in runs.items():
        print(f"{k}: first t with no particle {v['first_t_no_particle']}, single in window {v['single_in_window']}")
    for k in ("contrast_b28.50", "contrast_b28.00_mg6"):
        v = out[k]
        print(f"{k}: c = {np.round(v['c'], 2).tolist()}, mean {v['mean']:.2f}, t {v['t']:.2f}")


if __name__ == "__main__":
    main()
