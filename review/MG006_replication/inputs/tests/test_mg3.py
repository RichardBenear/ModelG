"""Tests for src/mg3.py (MG3 library). Run: python tests/test_mg3.py  (the resume test takes ~1 min)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import mg3
from runner import reference
from measure import smooth_periodic, detect


def test_first_passage_cases():
    s = mg3.first_passage(20.5, 70.5)
    assert s["startup"] and not s["event"] and s["exposure"] == 0
    s = mg3.first_passage(120.0, 170.0)
    assert s["event"] and s["exposure"] == 70.0 and s["in_w1"] and s["expo_w1"] == 70.0 and s["expo_w2"] == 0
    s = mg3.first_passage(300.0, 350.0)
    assert s["event"] and s["exposure"] == 250.0 and not s["in_w1"] and s["expo_w1"] == 150.0 and s["expo_w2"] == 100.0
    s = mg3.first_passage(480.0, 500.0)               # censoring boundary: detected at T - 20 counts
    assert s["event"] and s["exposure"] == 430.0
    s = mg3.first_passage(None, 500.0)                # right-censored at 480
    assert not s["event"] and s["exposure"] == 430.0 and s["complete"]
    s = mg3.first_passage(None, 300.0)                # unfinished run is not complete
    assert not s["complete"]


def test_garwood():
    lo, hi = mg3.garwood(0)
    assert lo == 0 and abs(hi - 3.6889) < 1e-3
    lo, hi = mg3.garwood(10)
    assert abs(lo - 4.7954) < 1e-3 and abs(hi - 18.3904) < 1e-3


def _synthetic_runs(k, n, rng, hazard_t=None):
    out = []
    for _ in range(n):
        t = mg3.T_B + rng.exponential(1 / (k * mg3.L)) if hazard_t is None else hazard_t(rng)
        tf = t if t <= mg3.T_C else None
        out.append(mg3.first_passage(tf, 500.0))
    return out


def test_cell_rate_recovers_k():
    rng = np.random.default_rng(1)
    runs = _synthetic_runs(1e-4, 4000, rng)
    r = mg3.cell_rate(runs)
    assert abs(r["k"] / 1e-4 - 1) < 0.05 and r["k_lo"] < 1e-4 < r["k_hi"]
    assert r["stationary_ok"]


def test_stationary_check_rejects_nonconstant_hazard():
    rng = np.random.default_rng(2)
    runs = _synthetic_runs(None, 400, rng, hazard_t=lambda g: 50 + g.uniform(0, 30))   # all early
    runs += _synthetic_runs(None, 400, rng, hazard_t=lambda g: 1000.0)                  # never
    assert not mg3.cell_rate(runs)["stationary_ok"]


def test_rate_ratio():
    r = mg3.rate_ratio(30, 1000.0, 30, 1000.0)
    assert r["lo"] < 1 < r["hi"] and abs(r["ratio"] - 1) < 1e-12
    r = mg3.rate_ratio(90, 1000.0, 30, 1000.0)
    assert r["lo"] > 1


def test_activated_fit_exact():
    cells = [(om, 50, np.exp(-2.0 - 0.03 * om)) for om in (100, 150, 200)]
    f = mg3.activated_fit(cells)
    assert abs(f["S"] - 0.03) < 1e-12 and f["r2"] > 0.999999 and f["ci"][1] < 0


def _c(om, cnt, exc=False, ev=None):
    return dict(Om=om, count=cnt, excluded=exc, any_event=(cnt > 0) if ev is None else ev)


def test_select_partial_plus_neighbours():
    sel, _ = mg3.select_stage_b([_c(300, 0), _c(250, 0), _c(200, 2), _c(150, 4)])
    assert sel == [250, 200, 150], sel


def test_select_transition_without_partial():
    sel, _ = mg3.select_stage_b([_c(300, 0), _c(250, 0), _c(200, 4)])
    assert sel == [300, 250, 200], sel


def test_select_top_and_none_and_excluded():
    assert mg3.select_stage_b([_c(300, 4)])[0] == [300]
    assert mg3.select_stage_b([_c(300, 0), _c(250, 0)])[0] == []
    assert mg3.select_stage_b([_c(300, 0, exc=True)])[0] == []
    sel, _ = mg3.select_stage_b([_c(300, 0), _c(250, 1), _c(200, 2), _c(150, 0, exc=True)])
    assert 150 not in sel and sel == [300, 250, 200], sel


def test_select_more_than_four_partials():
    sel, _ = mg3.select_stage_b([_c(o, 2) for o in (300, 250, 200, 150, 120, 100)])
    assert sel == [200, 150, 120, 100], sel


def test_select_events_only_after_150():
    sel, _ = mg3.select_stage_b([_c(300, 0), _c(250, 0, ev=True), _c(200, 0, ev=True)])
    assert sel == [300, 250, 200], sel


def test_wavelength_core_on_reference():
    ref = reference(29.0, 0.2)
    lam = mg3.wavelength_core(ref["x"], ref["phi"][1])
    assert abs(lam / ref["lam"] - 1) < 0.02, (lam, ref["lam"])


def test_nucleated_profile_synthetic():
    """A soliton appearing at x = 7 at t = 100 in an otherwise empty field: eligible profile equal to phi*."""
    ref = reference(29.0, 0.2)
    x = ref["x"]
    st = dict(times=[], dets=[], phi=[])
    for s in range(1, 401):
        t = 0.5 * s
        phi = np.zeros((3, 500)) if t < 100 else np.roll(ref["phi"], int(round(7 / 0.2)), axis=1)
        st["times"].append(t); st["phi"].append(phi.astype(np.float32))
        st["dets"].append(detect(smooth_periodic(phi[1], 0.2), x, 0.2, ref["thr"], ref["lam"]))
    from measure import track
    tr = track(st["dets"], np.array(st["times"]), 100.0, x0=None)
    q = tr["tracks"][0]
    st["first"] = dict(t_start=q["t_start"], track=q["id"])
    p = mg3.nucleated_profile(st)
    assert p["eligible"] and p["n_per_block"] == [10, 10, 10, 11], p["n_per_block"]
    win = np.abs(x) < 10
    err = np.linalg.norm(p["blocks"].mean(0)[1] - ref["phi"][1][win]) / np.linalg.norm(ref["phi"][1][win])
    assert err < 1e-3, err


def test_resume_is_deterministic():
    key = (99, 0, 0, 0, 7)
    a = mg3.sim_run(mg3.new_state(29.0, 1e6, 6.666666666666667e-4, key), until=20.0)
    b = mg3.new_state(29.0, 1e6, 6.666666666666667e-4, key)
    b = mg3.sim_run(b, until=10.0); b = mg3.sim_run(b, until=20.0)
    assert np.array_equal(a["nX"], b["nX"]) and a["times"] == b["times"] and a["births"] == [] == b["births"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for f in fns:
        f(); print("ok", f.__name__, flush=True)
    print(f"{len(fns)} tests passed")
