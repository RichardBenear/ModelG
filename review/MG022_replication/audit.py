"""Run with Python, numpy and scipy. Argument: path to a ModelG checkout at 9032882."""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import t as student


def stats(values):
    mean = values.mean()
    se = values.std(ddof=1) / np.sqrt(len(values))
    return dict(values=values.tolist(), mean=mean, t=mean/se,
                ci95=(mean + np.array([-1, 1]) * student.ppf(.975, 3) * se).tolist())


root = Path(sys.argv[1]) / "results/MG6"
evaluation = json.loads((root / "evaluation.json").read_text())
rows = {}
for path in sorted((root / "runs").glob("*.json")):
    key = path.stem
    rec = json.loads(path.read_text())
    assert rec["status"] == "COMPLETE" and not rec["flags"] and rec["t_reached"] == 200
    assert (rec["delta"], rec["dt"], rec["grid"]) == (.1, .05, [398, 402])
    with np.load(root / "checks" / (key + ".npz")) as d:
        np.testing.assert_array_equal(d["t"], np.arange(0, 201, 10))
        assert np.all(d["n_parts"] == 1) and len(d["tick_t"]) == 1999
        np.testing.assert_array_equal(np.diff(d["tick_N"], prepend=18901),
                                      d["tick_n_acc"] - d["tick_n_dead"])
        assert np.unpackbits(d["accepted_bits"]).sum() == rec["n_accepted"] == d["tick_n_acc"].sum()
        assert rec["kappa"] * rec["max_phiG"] <= .4 and rec["max_p"] <= 1
        area = rec["Lx"] * rec["Ly"]
        disc = np.pi * 1.5**2
        counts = d["m_N_R"][:, 1]
        window = d["t"] >= 100
        y = (counts - 18901/area * disc)[window].mean()
        yc = (counts - d["m_N"]/area * disc)[window].mean()
        p = d["m_P0"][window, 1].mean()
        L = np.array([rec["Lx"], rec["Ly"]])
        steps = np.diff(d["centre"], axis=0)
        track = np.vstack([np.zeros(2), np.cumsum(steps - L*np.round(steps/L), axis=0)])
        mobile = bool(np.linalg.norm(track, axis=1).max() >= 1.5)
        lags = np.arange(1, 11)
        msd = [np.square(track[j:] - track[:-j]).sum(axis=1).mean() for j in lags]
        dp = np.dot(msd, lags*10) / (4*np.sum((lags*10)**2))
        pub = evaluation["runs"][key]
        np.testing.assert_allclose([y, yc, dp], [pub["Y_D0"], pub["Y_Dc"], pub["cls"]["D_p"]], atol=1e-10)
        assert mobile == (pub["cls"]["mobility"] == "mobile")
        for name in ["M_minus", "M_plus", "P0", "Pc"]:
            a, b = d["m_"+name], d["x_"+name]
            if name.startswith("P"):
                a = a[:, 1]
            assert np.all(abs(a-b) <= .01*np.maximum(abs(a), abs(b)) + (.5 if name.startswith("P") else 0))
        edges = np.array([0,.5,1,1.5,2,2.5,3,4,5,6,8])
        ann = d["m_ann_count"][window].mean(axis=0)/(np.pi*np.diff(edges**2))-18901/area
        rows[key] = dict(y=y, yc=yc, ratio=y/p if p else None, mobile=mobile, dp=dp,
                         ann=ann.tolist(), phi=d["m_ann_phiG"][window].mean(axis=0).tolist())


def vector(arm, field):
    return np.array([rows[f"{arm}_D{s}_p{p}"][field] for s in range(2) for p in range(2)])


contrasts, falloff = {}, {}
for strength in (2, 4):
    cp, sh = f"CP{strength}", f"SH{strength}"
    for null in (sh, "DEC"):
        for field in ("y", "yc"):
            contrasts[f"{cp}-{null}_{field}"] = stats(vector(cp, field) - vector(null, field))
    profile = vector(cp, "phi").mean(axis=0)
    mids = (edges[1:] + edges[:-1])/2
    j = np.where((profile[:-1] < 0) & (profile[1:] >= 0))[0][0]
    zero = np.interp(0., profile[j:j+2], mids[j:j+2])
    eligible = np.flatnonzero((edges[:-1] >= 1) & (edges[1:] <= zero))
    falloff[cp] = dict(zero=zero, eligible=[edges[i:i+2].tolist() for i in eligible],
                       annular_t=[stats(v)["t"] for v in (vector(cp,"ann")-vector(sh,"ann"))[:,eligible].T])
    if strength == 2:
        for null in (sh, "DEC"):
            a = contrasts[f"{cp}-{null}_y"]
            b = evaluation["primary"][f"{cp}_minus_{null}"]
            np.testing.assert_allclose(a["values"], b["c"])
            np.testing.assert_allclose(a["t"], b["t"])
            assert max(a["values"]) < 0 and a["t"] <= -2.353

with (root / "budget_ledger.csv").open() as f:
    ledger = list(csv.DictReader(f))
assert len(rows) == len(ledger) == 20 and {r["key"] for r in ledger} == set(rows)
print(json.dumps(dict(contrasts=contrasts, falloff=falloff, absolute_cp2=stats(vector("CP2","y")),
                      ratios=vector("CP2","ratio").tolist(), mobile=[k for k,v in rows.items() if v["mobile"]],
                      core_hours=sum(float(r["wall_s"]) for r in ledger)/3600), indent=2))
