"""MG003c: calibration of the COMPLETE pre-registered decision procedures (protocol v1.3 §5.3), including
survival filtering with a minimum of 12 eligible runs, A1 diagnostics on block profiles, and the single
12-run extension; plus the Welch two-sample bound with unequal noise / unequal M and deterministic
normalisation; plus HELD-OUT families that violate A1 (rare run-level and block-level mixtures,
including the reviewer's 0.9/0.1 example), to document behaviour when A1 fails.
Synthetic run i: K=6 block profiles B_ib = phi + beta*u + r_i + e_ib (reference |phi| = 1, d = 100).
Writes results/MG003c/procedure_calibration.txt and .json."""
import json, sys
import numpy as np
sys.path.insert(0, "src")
from stats_mg2 import decide_cell, welch_diff_bounds, a1_diagnostics

rng = np.random.default_rng(20261030)
d, K, EPS = 100, 6, 0.01
phi = rng.standard_normal(d); phi /= np.linalg.norm(phi)
u = rng.standard_normal(d); u -= (u @ phi) * phi; u /= np.linalg.norm(u)
Q, _ = np.linalg.qr(np.column_stack([u, rng.standard_normal((d, 2))]))

def make_blocks(fam, n, beta, s_run):
    """n runs x K blocks. s_run = RMS norm of the run-mean fluctuation."""
    sb = s_run * np.sqrt(K)                                  # block noise RMS so run mean has RMS s_run
    if fam == "gauss_iso":
        e = rng.standard_normal((n, K, d)) * sb / np.sqrt(d); off = np.zeros((n, d))
    elif fam == "gauss_rank3":
        e = (rng.standard_normal((n, K, 3)) * sb / np.sqrt(3)) @ Q.T; off = np.zeros((n, d))
    elif fam == "t3_blocks":
        g = rng.standard_normal((n, K, 3)); w = rng.chisquare(3, (n, K, 1)) / 3
        e = (g / np.sqrt(w) / np.sqrt(3.0) * sb / np.sqrt(3)) @ Q.T; off = np.zeros((n, d))
    elif fam in ("mix_run_0.9", "mix_run_0.97"):            # A1-violating: rare run-level regime along u
        p = 0.9 if fam.endswith("0.9") else 0.97
        big = 0.06 if p == 0.9 else 0.2                       # reviewer: 0.5% vs 6.5% (offset 6%)
        off = np.outer(np.where(rng.random(n) < p, 0.0, big) - (1 - p) * big, u)
        e = rng.standard_normal((n, K, d)) * 0.0005 * np.sqrt(K) / np.sqrt(d)
    elif fam == "mix_block_0.1":                            # A1-violating: rare block-level jumps
        jump = np.where(rng.random((n, K)) < 0.1, 0.06 * np.sqrt(K), 0.0) - 0.1 * 0.06 * np.sqrt(K)
        e = jump[..., None] * u + rng.standard_normal((n, K, d)) * 0.0005 * np.sqrt(K) / np.sqrt(d)
        off = np.zeros((n, d))
    return phi + beta * u + off[:, None, :] + e

def runs_source(fam, beta, s_run, surv):
    def get(nreq, cap=24, state={"att": 0}):
        out = []
        while len(out) < nreq and state["att"] < cap:
            state["att"] += 1
            if rng.random() < surv:
                out.append(make_blocks(fam, 1, beta, s_run)[0])
        return (np.array(out) if out else np.zeros((0, K, d))), state["att"]
    return get

lines, js = [], {"D": [], "E": [], "F": []}
R = 600
GFAM, HFAM = ["gauss_iso", "gauss_rank3", "t3_blocks"], ["mix_run_0.9", "mix_run_0.97", "mix_block_0.1"]
# ---------------- Part D: full single-cell procedure
for fam in GFAM + HFAM:
    for surv in (1.0, 0.7):
        for beta in (0.0, 0.005, 0.009, 0.010, 0.011, 0.02):
            for s_run in ((0.001, 0.003, 0.01) if fam in GFAM else (None,)):
                cnt = {"PASS": 0, "FAIL": 0, "INCONCLUSIVE": 0}; reasons = {}
                for _ in range(R):
                    get = runs_source(fam, beta, s_run if s_run else 0.0, surv)
                    get.__defaults__[1]["att"] = 0
                    dcs, info = decide_cell(get, phi, EPS)
                    cnt[dcs] += 1
                    if "reason" in info: reasons[info["reason"]] = reasons.get(info["reason"], 0) + 1
                js["D"].append(dict(fam=fam, surv=surv, beta=beta, s_run=s_run, **{k: v / R for k, v in cnt.items()}, reasons=reasons))
summ = {}
for fam in GFAM + HFAM:
    rows = [r for r in js["D"] if r["fam"] == fam]
    fp = max(r["PASS"] for r in rows if r["beta"] >= EPS)
    ff = max(r["FAIL"] for r in rows if r["beta"] <= EPS)
    summ[fam] = (fp, ff)
lines.append(f"Part D: full single-cell procedure (min 12 eligible of <= 24 attempted; A1 diagnostics; one 12-run extension), {R} reps/cell")
lines.append("  A1-satisfying families (Gaussian; t3 blocks is a stress case):")
for fam in GFAM:
    lines.append(f"    {fam:12s}: worst false-PASS (true bias >= 1%) {summ[fam][0]:.3f}; worst false-FAIL (true bias <= 1%) {summ[fam][1]:.3f}")
lines.append("  HELD-OUT A1-violating families (documented, not covered by the confidence claim):")
for fam in HFAM:
    lines.append(f"    {fam:12s}: worst false-PASS {summ[fam][0]:.3f}; worst false-FAIL {summ[fam][1]:.3f}")
g0 = [r for r in js["D"] if r["fam"] in ("gauss_iso", "gauss_rank3") and r["beta"] in (0.0, 0.005) and r["s_run"] == 0.001 and r["surv"] == 1.0]
lines.append(f"  power check (Gaussian, bias 0-0.5%, run noise 0.1%): P(PASS) min {min(r['PASS'] for r in g0):.3f}")
print("\n".join(lines), flush=True)

# ---------------- Part F: A1 diagnostic false-alarm and detection
Rf = 1000
for fam in GFAM + HFAM:
    fl = 0
    for _ in range(Rf):
        fl += not a1_diagnostics(make_blocks(fam, 12, 0.005, 0.003), phi)["ok"]
    js["F"].append(dict(fam=fam, flag_rate=fl / Rf))
    lines.append(f"  A1 diagnostics flag rate, {fam:12s} (12 runs, bias 0.5%, run noise 0.3%): {fl / Rf:.3f}")
print("\n".join(lines[-6:]), flush=True)

if "--skip-e" in sys.argv:
    prev = json.load(open("results/MG003c/procedure_calibration_partE.json"))
    js["E"] = prev["E"]
    lines += prev["lines"]
    import os; os.makedirs("results/MG003c", exist_ok=True)
    open("results/MG003c/procedure_calibration.txt", "w").write("\n".join(lines) + "\n")
    json.dump(js, open("results/MG003c/procedure_calibration.json", "w"), indent=1)
    sys.exit(0)
# ---------------- Part E: Welch two-sample, unequal noise and M, deterministic normalisation
lines.append("Part E: Welch two-sample bound, deterministic normalisation (|phi| = 1), unequal noise (x3) and M (8/12)")
Re = 2000
for (MA, MB) in ((8, 12), (12, 8), (12, 12)):
    for ratio in (1 / 3, 1.0, 3.0):
        for margin, deltas, sig in ((0.10, (0.0, 0.05, 0.10, 0.15), (0.01, 0.03)), (0.0025, (0.0, 0.0015, 0.0025, 0.004), (0.0005, 0.001))):
            worst_cov, worst_fe = 1.0, 0.0
            for fam in ("gauss_iso", "gauss_rank3", "t3_blocks"):
                for delta in deltas:
                    for s in sig:
                        cov = fe = 0
                        for _ in range(Re):
                            A = make_blocks(fam, MA, delta, s).mean(1)
                            B = make_blocks(fam, MB, 0.0, s * ratio).mean(1)
                            b = welch_diff_bounds(A, B, 1.0)
                            cov += b["lower"] - 1e-15 <= delta <= b["upper"] + 1e-15
                            fe += (b["upper"] < margin) and (delta >= margin)
                        worst_cov = min(worst_cov, cov / Re); worst_fe = max(worst_fe, fe / Re)
            js["E"].append(dict(MA=MA, MB=MB, sd_ratio=ratio, margin=margin, worst_cov=worst_cov, worst_false_equiv=worst_fe))
            lines.append(f"  MA={MA} MB={MB} sdB/sdA={ratio:.2f} margin {100*margin:.2f}%: worst coverage {worst_cov:.3f}, worst false-equivalence {worst_fe:.3f}")
            print(lines[-1], flush=True)
import os; os.makedirs("results/MG003c", exist_ok=True)
open("results/MG003c/procedure_calibration.txt", "w").write("\n".join(lines) + "\n")
json.dump(js, open("results/MG003c/procedure_calibration.json", "w"), indent=1)
json.dump({"E": js["E"], "lines": [l for l in lines if l.startswith("Part E") or l.startswith("  MA=")]},
          open("results/MG003c/procedure_calibration_partE.json", "w"), indent=1)
