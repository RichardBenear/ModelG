"""Independent records audit. Usage: python audit.py PATH_TO_REVIEWED_MODELG_CHECKOUT."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(sys.argv[1]) / "src"))
import recompute_mg6 as replay

repo = Path(sys.argv[1])
root = repo / "results/MG7"
ev = json.loads((root / "evaluation.json").read_text())
rows, streams = {}, {}
for path in sorted((root / "runs").glob("*.json")):
    rec = json.loads(path.read_text())
    key = path.stem
    assert rec["status"] == "COMPLETE" and rec["t_reached"] == 200 and not rec["flags"]
    assert rec["b"] == rec["P_b"] == rec["spec"]["b"]
    assert (rec["delta"], rec["dt"], rec["grid"]) == (.1, .05, [398, 402])
    assert not rec["start_field"]["shift_applied"]
    rng = np.random.default_rng(rec["replay"]["seed"])
    life = rng.exponential(5, 18901)
    size = rng.poisson(.2 * 18901 * np.exp(.4) * 200)
    born = np.sort(rng.random(size) * 200)
    x = rng.random(size) * rec["Lx"]
    y = rng.random(size) * rec["Ly"]
    mark = rng.random(size)
    lifetime = rng.exponential(5, size)
    digest = hashlib.sha256(np.array([18901], dtype="<i8").tobytes())
    for array in (life, born, x, y, mark, lifetime):
        digest.update(np.asarray(array, dtype="<f8").tobytes())
    assert digest.hexdigest() == rec["replay"]["sha256"]
    start = (rec["spec"]["slot"], rec["spec"]["pl"])
    streams.setdefault(start, set()).add(digest.hexdigest())
    with np.load(root / "checks" / (key + ".npz")) as d:
        accepted = np.unpackbits(d["accepted_bits"])[:size].astype(bool)
        assert accepted.sum() == rec["n_accepted"] == d["tick_n_acc"].sum()
        deaths = born + lifetime
        assert not np.any(accepted & (deaths <= np.ceil(born/.1-1e-12)*.1))
        times = d["tick_t"]
        np.testing.assert_allclose(times, np.arange(1, 2000)*.1, atol=1e-12)
        population = (18901 - np.searchsorted(np.sort(life), times, side="right")
                      + np.searchsorted(born[accepted], times, side="right")
                      - np.searchsorted(np.sort(deaths[accepted]), times, side="right"))
        np.testing.assert_array_equal(population, d["tick_N"])
        np.testing.assert_array_equal(np.diff(population, prepend=18901),
                                      d["tick_n_acc"]-d["tick_n_dead"])
        assert population[-1] == rec["N_final"]
        assert rec["max_phiG"] <= .2 and rec["max_p"] <= 1
        np.testing.assert_array_equal(d["t"], np.arange(0,201,10))
        n, t = d["n_parts"], d["t"]
        assert np.all((n == 0) | (n == 1)) and n[0] == 1
        zero = t[n == 0]
        assert not len(zero) or np.all(n[t >= zero[0]] == 0)
        assert np.all(d["smin"][n == 0] >= rec["theta"])
        win = t >= 100
        row = dict(first_zero=float(zero[0]) if len(zero) else None,
                   final_smoothed_X_min=float(d["smin"][-1]), single_window=bool(np.all(n[win] == 1)),
                   alt_single_window=ev["runs"][key]["alt_threshold"]["single_in_window"])
        # Independently verify non-detection for each alternative-threshold zero check.
        alt = ev["runs"][key]["alt_threshold"]
        alt_counts = [len(replay.detect(d["pos"][a:b],d["s"][a:b],alt["theta"],rec["Lx"],rec["Ly"]))
                      for a,b in zip(d["off"][:-1],d["off"][1:])]
        assert alt_counts == alt["counts"]
        assert np.all(d["smin"][np.array(alt["counts"]) == 0] >= alt["theta"])
        for name in ("P0", "Pc", "M_minus", "M_plus"):
            a, b = d["m_"+name], d["x_"+name]
            if name.startswith("P"):
                a = a[:,1]
            mask = n == 1
            assert np.all(abs(a[mask]-b[mask]) <= .01*np.maximum(abs(a[mask]),abs(b[mask])) +
                          (.5 if name.startswith("P") else 0))
        if row["single_window"]:
            baseline = 18901/(rec["Lx"]*rec["Ly"])*np.pi*1.5**2
            row["Y"] = float((d["m_N_R"][win,1]-baseline).mean())
            np.testing.assert_allclose(row["Y"],ev["runs"][key]["Y_D0"],atol=1e-10)
            row["late_minus_early"] = float(d["m_N_R"][(t>=160),1].mean() -
                                            d["m_N_R"][(t>=100)&(t<=140),1].mean())
        rows[key] = row
assert all(len(v) == 1 for v in streams.values()) and len(streams) == 4
assert len(rows) == 16 and not all(r["single_window"] for r in rows.values())
assert ev["gate"] == "NOT SCORED (feedback)"
starts = [f"D{s}_p{p}" for s in range(2) for p in range(2)]
contrast = np.array([rows[f"CP2_b28.50_{s}"]["Y"]-rows[f"SH2_b28.50_{s}"]["Y"] for s in starts])
desc = json.loads((root/"descriptive.json").read_text())
np.testing.assert_allclose(contrast,desc["contrast_b28.50"]["c"])
part0 = {}
for b in ("27.50","28.50"):
    with np.load(root/f"part0/check_b{b}_d0.1_dt0.05.npz") as a:
        comparisons = {}
        for label, name, tol, dtol in [
            ("step",f"check_b{b}_d0.1_dt0.025_forced.npz",2,.02),
            ("tick",f"check_b{b}_d0.05_dt0.05.npz",6,.03)]:
            with np.load(root/"part0"/name) as z:
                np.testing.assert_array_equal(a["n_parts"],z["n_parts"])
                depth = float(np.max(abs(a["smin"]-z["smin"])/abs(a["smin"])))
                diff = float(abs(a["m_D0"][a["t"]>=50,1].mean()-z["m_D0"][z["t"]>=50,1].mean()))
                assert diff <= tol and depth <= dtol
                comparisons[label] = dict(count_mean_difference=diff, relative_depth_difference=depth)
        part0[b] = comparisons
with (root/"budget_ledger.csv").open() as f:
    ledger = list(csv.DictReader(f))
assert len(ledger) == 16 and {r["key"] for r in ledger} == set(rows)
print(json.dumps(dict(rows=rows, gate=ev["gate"], contrast_b285=contrast.tolist(),
                      mean_b285=float(contrast.mean()), t_b285=float(contrast.mean()/(contrast.std(ddof=1)/2)),
                      part0=part0, core_hours=sum(float(r["wall_s"]) for r in ledger)/3600),indent=2))
