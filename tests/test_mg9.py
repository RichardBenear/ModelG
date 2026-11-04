"""MG9 tests (spec/MG9_proposal.md draft 3 §9; amendment 1 for the envelope). No Model G nonlinear MG9 run is made:
the only integration is the MG003 b = 29 reference (wrapper check). Synthetic fields and tracks otherwise.
    python tests/test_mg9.py
"""
import os, sys, time, tempfile
import numpy as np
HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import landau_mg9 as L
import run_mg9 as R
import evaluate_mg9 as E

LAM = L.LAM
K = 2 * np.pi / LAM


def test_seed_scale_must_be_explicit():
    try:
        R.mg9_ode(29.0, 500, 0.2, 1.0, np.zeros(1500))
    except TypeError:
        return
    raise AssertionError("mg9_ode accepted a call without seed_scale")


def test_wrapper_matches_lattice_ode_1d():
    from lattice import lattice_ode_1d, P17
    x, pG, pX, pY, res = lattice_ode_1d(P17, 500, 0.2, 400.0)
    r = R.mg9_ode(29.0, 500, 0.2, 400.0, np.zeros(1500), seed_scale=1.0)
    assert r["complete"] and r["success"] and r["t_reached"] == 400.0
    assert np.abs(r["y_final"] - np.concatenate([pG, pX, pY])).max() <= 1e-12
    assert abs(r["residual"] - res) <= 1e-12
    r0 = R.mg9_ode(29.0, 500, 0.2, 30.0, np.zeros(1500), seed_scale=0.0)
    assert np.abs(r0["y_final"]).max() <= 1e-12                   # no seed: the uniform state stays uniform


def test_incomplete_run_is_not_complete():
    """A solver that stops early (success False, t short of t_end) must give complete = False after all retries."""
    real = R.solve_ivp
    class Fake:
        success, t, y = False, np.array([0.0, 2.0]), np.zeros((150, 2))
    calls = []
    R.solve_ivp = lambda *a, **k: (calls.append(k.get("first_step")), Fake)[1]
    try:
        r = R.mg9_ode(29.0, 50, 0.2, 5.0, np.zeros(150), seed_scale=0.0)
    finally:
        R.solve_ivp = real
    assert not r["complete"] and not r["success"] and r["t_reached"] == 2.0 and calls == [None, 1e-3, 1e-2]


def test_carry_rules():
    N = 7
    y = np.arange(3 * N, dtype=float)
    z = R.carry_concentration(y, N, 29.0, 28.5)
    assert np.array_equal(z[:2 * N], y[:2 * N])
    assert np.allclose(L.Y0(28.5) + z[2 * N:], L.Y0(29.0) + y[2 * N:], atol=1e-14)


def test_swift_hohenberg_coefficient():
    for c in (0.0, 0.5, 1.0, np.sqrt(27 / 38)):
        assert abs(L.l3_swift_hohenberg(c) - (19 * c * c / 18 - 0.75)) <= 1e-8
    assert abs(L.l3_swift_hohenberg(np.sqrt(27 / 38))) <= 1e-12


def test_v1_thresholds():
    bc, kap = L.bc_continuous()
    assert abs(bc - 30.04407748) <= 1e-4
    t = {d: L.box_thresholds(d) for d in L.DOMAINS}
    assert abs(t["D02"]["onset"] - 30.04577) <= 1e-4 and t["D02"]["onset_mode"] == 34
    assert abs(t["D02L"]["onset"] - 30.044468) <= 1e-4 and t["D02L"]["onset_mode"] == 67
    assert abs(t["D02L"]["selected_threshold"] - 30.045770) <= 1e-4
    assert abs(t["D01"]["onset"] - 30.04757) <= 1e-4 and t["D01"]["onset_mode"] == 33
    assert abs(L.alpha(L.K34, 0.2) - 0.1981) <= 1e-3
    assert abs(L.alpha(2 * np.pi * 33 / 100, 0.1) - 0.1906) <= 1e-3


def _patch(x, a, ell, ph):
    w = 0.5 * (np.tanh((x + a) / ell) - np.tanh((x - a) / ell))
    return np.cos(K * x + ph) * w


def _envelope_validation(env):
    """§6.2 synthetic validation; returns (cosine range, worst edge error / lam) over both meshes."""
    lo, hi, worst = 9.0, 0.0, 0.0
    for N, h in ((500, 0.2), (1000, 0.1)):
        x = (np.arange(N) - N // 2) * h
        for ph in np.linspace(0, 2 * np.pi, 10, endpoint=False):
            e = env(np.cos(K * x + ph), h)
            lo, hi = min(lo, e.min()), max(hi, e.max())
            for ell in (LAM / 2, LAM / 8, 2 * LAM):
                for a in (25.0, 17.3):
                    E_ = env(_patch(x, a, ell, ph), h)
                    comps = R.components(E_ >= 0.5, h)
                    c = [c for c in comps if N // 2 in set(c.tolist())][0]
                    tr = R.Tracker(N, h, 1.0)
                    xl, xr = tr._ends(E_, c)
                    worst = max(worst, abs(xl + a), abs(xr - a))
    return lo, hi, worst / LAM


def test_registered_rms_envelope_fails_its_validation():
    """Documents why amendment 1 exists: draft 3's single-boxcar RMS envelope misses the lam/10 edge criterion."""
    lo, hi, worst = _envelope_validation(R.envelope_rms)
    assert 0.98 <= lo and hi <= 1.02
    assert worst > 0.10


def test_envelope_validation_amendment_1():
    lo, hi, worst = _envelope_validation(R.envelope)
    assert 0.98 <= lo and hi <= 1.02, (lo, hi)
    assert worst <= 0.10, worst


def test_tracker_on_moving_patch():
    N, h = 500, 0.2
    x = (np.arange(N) - N // 2) * h
    tr = R.Tracker(N, h, 1.0)
    r0 = tr.sample(_patch(x, 25.0, LAM / 2, 0.0))
    r1 = tr.sample(_patch(x, 27.0, LAM / 2, 0.0))
    assert r0["status"] == r1["status"] == "ok"
    assert abs((r1["xR"] - r0["xR"]) - 2.0) < 0.2 and abs((r0["xL"] - r1["xL"]) - 2.0) < 0.2
    r2 = tr.sample(_patch(x, 27.0, LAM / 2, 0.0) + np.where(np.abs(x - 45) < 3, np.cos(K * x), 0))
    assert r2["remote"]
    assert R.Tracker(N, h, 1.0).sample(np.cos(K * x))["status"] == "fill"
    assert R.Tracker(N, h, 1.0).sample(np.zeros(N))["status"] == "loss"


def _track(xLout, xRout, t_end=600.0, dt=5.0, status=None, n=None, breakup_at=None, remote_at=None):
    ts = np.arange(0, t_end + 1e-9, dt)
    out = []
    for i, t in enumerate(ts):
        st = status(t) if status else "ok"
        out.append(dict(t=float(t), status=st, xL=-25.0 - xLout(t), xR=25.0 + xRout(t),
                        n=(n(t) if n else 17), breakup=bool(breakup_at is not None and t >= breakup_at),
                        remote=bool(remote_at is not None and t >= remote_at)))
    return out


def test_classification_procedure():
    c = lambda tr, T=600.0, res=1e-7, ok=True: R.classify(tr, T, res, ok)[0]
    zero = lambda t: 0.0
    assert c(_track(zero, zero)) == "ARRESTED"
    assert R.classify(_track(zero, zero), 600, 1e-7, True)[1]["margin"]
    assert c(_track(zero, zero), ok=False) == "UNRESOLVED-NUMERICAL"
    assert c(_track(zero, zero), res=1e-3) == "SLOW"
    inv = lambda t: 4 * LAM * t / 600
    assert c(_track(inv, inv)) == "INVADES"
    assert c(_track(lambda t: -inv(t), lambda t: -inv(t))) == "RETREATS"
    # ChatGPT r2 example: 2 lam out by t = 300, back to lam by t = 600 -> reversal, not invasion
    rev = lambda t: 2 * LAM * min(t, 300) / 300 - LAM * max(0, t - 300) / 300
    assert c(_track(rev, rev)) == "UNRESOLVED-REVERSAL"
    assert c(_track(inv, zero)) == "UNRESOLVED-ASYMMETRIC"
    assert c(_track(zero, lambda t: -inv(t))) == "UNRESOLVED-ASYMMETRIC"
    assert c(_track(inv, lambda t: -inv(t))) == "UNRESOLVED-ASYMMETRIC"          # drift of the whole patch
    assert c(_track(zero, zero, status=lambda t: "fill" if t >= 50 else "ok")) == "INVADES"
    assert c(_track(zero, zero, status=lambda t: "loss" if t >= 50 else "ok")) == "RETREATS"
    assert c(_track(zero, zero, status=lambda t: "loss" if 200 <= t < 300 else "ok")) == "UNRESOLVED-TEMPORARY"
    assert c(_track(zero, zero, status=lambda t: "fill" if 200 <= t < 300 else "ok")) == "UNRESOLVED-TEMPORARY"
    assert c(_track(zero, zero, breakup_at=400)) == "UNRESOLVED-INTERNAL"
    assert c(_track(zero, zero, remote_at=400)) == "UNRESOLVED-INTERNAL"
    assert c(_track(zero, zero, n=lambda t: 17 + int(t > 300) + int(t > 400))) == "SLOW"
    slow = lambda t: LAM / 1500 * t * 0.5                     # 0.2 lam over [100, 600], v = lam/3000: arrest, no margin
    lab, info = R.classify(_track(slow, slow), 600, 1e-7, True)
    assert lab == "ARRESTED" and not info["margin"] and R.needs_extension(lab, info)
    fast = lambda t: LAM / 1000 * t                            # v = lam/1000: too fast for arrest, too slow for motion
    assert c(_track(fast, fast)) == "SLOW"


def test_edge_and_outcome_logic():
    lab = {"A": "ARRESTED", "R": "RETREATS", "I": "INVADES", "N": "UNAVAILABLE-LOST", "U": "UNRESOLVED-SLOW"}
    def pts(s, b0=28.80, db=0.05):
        return {round(b0 + db * i, 4): lab[ch] for i, ch in enumerate(s)}
    bench = dict(b_s=dict(status="OK", bracket=(28.9290, 28.9298), sensitive=False),
                 b_r=dict(status="OK", bracket=(29.3961, 29.3969), sensitive=False))
    def out(s, **kw):
        p = pts(s, **kw)
        ed = R.arrest_edges(R.grid_classes(p)); ed["points"] = p
        return E.p2_outcome(True, bench, ed)
    # 28.80 R, 28.85 R, 28.90 R, 28.95 A ... 29.40 A, 29.45 I
    assert out("RRRAAAAAAAAAAI")[0] == "COINCIDES"
    assert out("RRRAAAAAAAAAAI")[1] == []
    assert out("NNNAAAAAAAAAAI")[0] == "OVERLAPS" and "CENSORED" in out("NNNAAAAAAAAAAI")[1]
    assert out("RRUAAAAAAAAAAI")[0] == "UNRESOLVED"
    assert out("AAAAAAAAAAAAAI")[0] == "UNRESOLVED"                               # edge at the scan end
    assert out("RRRAAARAAAAAAI")[0] == "UNRESOLVED"                               # disconnected
    assert out("RRRRRRRRRRRRRI")[0] == "NO ARREST OBSERVED"
    assert out("RRRRRRRRRRRRRU")[0] == "UNRESOLVED"
    assert out("RRRRRRRRRRRRRRRRAAI")[0] == "OUTSIDE"                               # arrest at 29.60-29.65
    assert out("RRRAAAAAAAAI")[0] == "COINCIDES"                                   # upper bracket 29.30-29.35
    assert out("RRRAAAAAAAI")[0] == "OVERLAPS"                                     # upper bracket 29.25-29.30: too far
    o = out("IRRAAAAAAAAAAI")
    assert "ANOMALOUS" in o[1]
    bad = dict(bench, b_s=dict(status="BENCHMARK-UNRESOLVED"))
    p = pts("RRRAAAAAAAAAAI"); ed = R.arrest_edges(R.grid_classes(p)); ed["points"] = p
    assert E.p2_outcome(True, bad, ed)[0] == "UNRESOLVED"
    assert E.p2_outcome(False, bench, ed)[0] == "UNRESOLVED"


def test_v2_fit_and_decision():
    sig, L3 = 1e-8, 3e-2
    def run(A0):
        t = np.arange(0, 3000.0, 1.0)
        y = (1 / A0 ** 2 + L3 / sig) * np.exp(-2 * sig * t) - L3 / sig
        keep = y >= 0.25 / A0 ** 2                                   # amplitude stop at 2 A0
        return dict(A0=A0, t=t[keep], A=1 / np.sqrt(y[keep]), stop_t=float(t[keep][-1]))
    f = L.fit_l3(run(0.2)["t"], run(0.2)["A"], sig)
    assert abs(f["L3"] - L3) <= 1e-9 and f["rel_change_invA2"] > 0.25
    runs = [run(a) for a in (0.05, 0.1, 0.2)]
    d = L.v2_decision([dict(r) for r in runs], L3 * 1.05, sig, repeat=(L3, L3 * 1.01))
    assert d["passed"] and d["n_identifiable"] >= 2
    assert not L.v2_decision([dict(r) for r in runs], -L3, sig)["passed"]
    assert not L.v2_decision([dict(r) for r in runs], L3 * 1.2, sig)["passed"]
    assert not L.v2_decision([dict(r) for r in runs], L3, sig, repeat=(L3, L3 * 1.05))["passed"]
    short = [dict(r, stop_t=30.0) for r in runs]
    assert not L.v2_decision(short, L3, sig)["passed"]


def test_p1_outcome_order():
    l3 = dict(i=1.0, ii=1.2, iii=0.9, iv=1.1)
    assert L.p1_outcome(True, True, l3, (1.0, 0.5, -0.5)) == "SUBCRITICAL"
    assert L.p1_outcome(True, True, {k: -v for k, v in l3.items()}, (1.0, 0.5, -0.5)) == "SUPERCRITICAL"
    assert L.p1_outcome(False, True, l3, (1.0, 0.5, -0.5)) == "NOT EVALUATED"
    assert L.p1_outcome(True, True, dict(l3, ii=-1), (1.0, 0.5, -0.5)) == "MIXED"
    assert L.p1_outcome(True, True, dict(l3, i=1e-6), (1.0, 0.5, -0.5)) == "DEGENERATE"


def test_fourier_resample_and_shift():
    N, h = 500, 0.2
    x = (np.arange(N) - N // 2) * h
    y = np.concatenate([np.cos(K * x + 0.3), 2 * np.cos(K * x), np.sin(2 * K * x)])
    y1 = R.fourier_resample(y, N, 1000)
    x1 = (np.arange(1000) - 500) * 0.1
    assert np.abs(y1[1000:2000] - 2 * np.cos(K * x1)).max() <= 1e-10
    ys = R.fourier_shift(y, N, h, LAM / 2)
    assert np.abs(ys[N:2 * N] + 2 * np.cos(K * x)).max() <= 1e-10


def test_kernel_weights_sum_to_lam_over_h():
    for h in (0.2, 0.1):
        assert abs(sum(R.kernel_weights(h).values()) - LAM / h) <= 1e-12


def test_ledger_cap_enforced():
    old = R.OUT
    with tempfile.TemporaryDirectory() as d:
        R.OUT = d
        try:
            for i in range(R.RUN_CAP):
                R._ledger_add(f"k{i}", 0.0)
            try:
                R._ledger_check()
            except R.Stop:
                pass
            else:
                raise AssertionError("cap not enforced")
        finally:
            R.OUT = old


if __name__ == "__main__":
    t0 = time.time()
    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]
    for name, fn in tests:
        t = time.time()
        fn()
        print(f"  ok  {name} ({time.time() - t:.1f} s)", flush=True)
    print(f"{len(tests)} tests passed in {time.time() - t0:.1f} s")
