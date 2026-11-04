"""MG9c tests (spec/MG9c_proposal.md draft 2 §§4-8). No Model G integration is made: synthetic tracks, and a fake
solver for the segmentation, snapshot, ledger and replay logic.
    python tests/test_mg9c.py
"""
import os, sys, time, tempfile, json
import numpy as np
HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import run_mg9c as C

LAM = C.LAM


def track(fL, fR, T, n=lambda t: 18, status=lambda t: "ok", breakup=lambda t: False, x0=29.0):
    """Synthetic track every 5 units: outward displacements fL(t), fR(t) from edges -x0, +x0."""
    out = []
    for t in np.arange(0.0, T + 1e-9, 5.0):
        st = status(t)
        rec = dict(t=float(t), n=n(t), status=st, breakup=breakup(t), remote=False, xL=None, xR=None)
        if st == "ok":
            rec["xL"], rec["xR"] = -x0 - fL(t), x0 + fR(t)
        out.append(rec)
    return out


def classify_main(tr, res=1e-7):
    return C.late_classify(tr, tr[-1]["t"], res, True, C.MAIN_W)[0]


def expo(t, tau=715.0, a=4.0):
    return a * (1 - np.exp(-t / tau))


def test_exponential_settling_is_settled():
    tr = track(expo, expo, 12000)
    assert classify_main(tr) == "SETTLED"


def test_constant_slow_speed_advances():
    v = 3.6e-4                                     # MG9's late speed: 2.2 over W3 >= lam/4
    tr = track(lambda t: v * t, lambda t: v * t, 12000)
    assert classify_main(tr) == "ADVANCES"


def test_constant_recession_recedes():
    v = -3.6e-4
    assert classify_main(track(lambda t: v * t, lambda t: v * t, 12000)) == "RECEDES"


def test_creep_between_thresholds_is_slow():
    v = (LAM / 20) / 6000.0                        # D3 = lam/20: above lam/40, below lam/4
    tr = track(lambda t: v * t, lambda t: v * t, 12000)
    assert classify_main(tr) == "SLOW"


def test_power_law_creep_at_threshold():
    # x = c t^0.5 with D3 just above / below lam/40
    d3 = lambda c: c * (np.sqrt(12000) - np.sqrt(6000))
    c_hi = 1.05 * (LAM / 40) / (np.sqrt(12000) - np.sqrt(6000))
    c_lo = 0.5 * (LAM / 40) / (np.sqrt(12000) - np.sqrt(6000))
    assert d3(c_hi) > LAM / 40
    assert classify_main(track(lambda t: c_hi * np.sqrt(t), lambda t: c_hi * np.sqrt(t), 12000)) == "SLOW"
    # below: speed over the last 1000 units c/(2 sqrt t) ~ 6e-5 > lam/40000 = 7.4e-5? check it settles only if v small
    lab = classify_main(track(lambda t: c_lo * np.sqrt(t), lambda t: c_lo * np.sqrt(t), 12000))
    v_end = c_lo / (2 * np.sqrt(11500))
    assert lab == ("SETTLED" if v_end <= LAM / 40000 else "SLOW")


def test_early_stop_settles():
    f = lambda t: min(t, 1000.0) * 1e-3
    assert classify_main(track(f, f, 12000)) == "SETTLED"


def test_one_interface_moving_is_asymmetric():
    v = 3.6e-4
    assert classify_main(track(lambda t: v * t, lambda t: 0.0, 12000)) == "UNRESOLVED-ASYMMETRIC"


def test_opposite_sign_translation_is_asymmetric():
    # (+0.3 lam, -0.3 lam) over W3: the patch translates
    v = 0.3 * LAM / 6000
    f = lambda t: v * max(t - 6000, 0)
    assert classify_main(track(f, lambda t: -f(t), 12000)) == "UNRESOLVED-ASYMMETRIC"


def test_early_fill_invades_and_temporary_loss():
    tr = track(expo, expo, 12000, status=lambda t: "fill" if t >= 50 else "ok")
    assert classify_main(tr) == "INVADES"
    tr = track(expo, expo, 12000, status=lambda t: "loss" if 50 <= t < 80 else "ok")
    assert classify_main(tr) == "UNRESOLVED-TEMPORARY"
    tr = track(expo, expo, 12000, status=lambda t: "loss" if t >= 20 else "ok")
    assert classify_main(tr) == "RETREATS"


def test_breakup_and_reversal():
    tr = track(expo, expo, 12000, breakup=lambda t: 4000 <= t < 4010)
    assert classify_main(tr) == "UNRESOLVED-INTERNAL"
    tr = track(expo, expo, 12000, breakup=lambda t: t < 50)        # before t = 100: not a disqualifier
    assert classify_main(tr) == "SETTLED"
    rev = lambda t: (LAM * 0.6 * t / 1000 if t < 1000 else LAM * 0.6 - LAM * 0.6 * min(t - 1000, 1000) / 1000)
    assert classify_main(track(rev, expo, 12000)) == "UNRESOLVED-REVERSAL"


def test_residual_and_cell_count_block_settled():
    tr = track(expo, expo, 12000)
    assert classify_main(tr, res=2e-5) == "SLOW"
    tr = track(expo, expo, 12000, n=lambda t: 18 if t < 5000 else 20)
    assert classify_main(tr) == "SLOW"
    tr = track(expo, expo, 12000, n=lambda t: 18 if t < 1000 else 19)  # change before W1: allowed
    assert classify_main(tr) == "SETTLED"


def test_speed_must_decrease():
    # a front that is slower in W1 than in W3 (but below lam/40 in W3 and above lam/400): not SETTLED
    f = lambda t: 0.0 if t < 6000 else (LAM / 50) * (t - 6000) / 6000
    tr = track(f, f, 12000)
    lab, info = C.late_classify(tr, 12000, 1e-7, True, C.MAIN_W)
    assert lab == "SLOW" and not info["decreasing"]["L"]


def test_extension_windows_and_perturbation_variant():
    tr = track(expo, expo, 30000)
    assert C.late_classify(tr, 30000, 1e-7, True, C.EXT_W)[0] == "SETTLED"
    tr = track(lambda t: expo(t, a=0.01), lambda t: expo(t, a=0.01), 6000)
    assert C.late_classify(tr, 6000, 1e-7, True, C.PERT_W)[0] == "SETTLED"
    assert C.late_classify(tr, 6000, 1e-7, False, C.PERT_W)[0] == "UNRESOLVED-NUMERICAL"


def test_trial_and_point_status():
    rec_s = dict(label="SETTLED", track=[dict(xL=-29.1, xR=29.0)])
    assert C.trial_status("a", rec_s) == "PASS"
    assert C.trial_status("a", dict(label="ADVANCES")) == "FAIL"
    assert C.trial_status("a", dict(label="UNRESOLVED-SLOW")) == "UNRESOLVED"
    assert C.trial_status("a", dict(label="UNRESOLVED-PREP")) == "UNRESOLVED"
    assert C.trial_status("a", None) == "UNRESOLVED"
    assert C.trial_status("pert", rec_s, (-29.1, 29.0 + 0.5 * LAM / 10)) == "PASS"
    assert C.trial_status("pert", rec_s, (-29.1, 29.0 + 2 * LAM / 10)) == "FAIL"
    assert C.point_status(dict(a="PASS", b="PASS", **{"c+": "PASS", "c-": "PASS"})) == "SETTLED-ROBUST"
    assert C.point_status(dict(a="PASS", b="FAIL", **{"c+": "UNRESOLVED", "c-": "PASS"})) == "NOT-ROBUST"
    assert C.point_status(dict(a="PASS", b="UNRESOLVED", **{"c+": "PASS", "c-": "PASS"})) == "ROBUSTNESS-UNRESOLVED"


def test_q2_interval_rule():
    lo, hi = C.BR
    assert C.q2_compare(29.39375, 29.40) == "NEAR b_r"              # contains the hypothetical 29.398 edge
    assert C.q2_compare(29.38750, 29.39375) == "NEAR b_r"
    assert C.q2_compare(29.35, 29.35625) == "BELOW b_r"
    assert C.q2_compare(29.375, 29.38125) == "BELOW b_r"
    assert C.q2_compare(29.38125, 29.3875) == "UNRESOLVED RELATIVE TO THE TOLERANCE"
    assert C.q2_compare(29.41, 29.42) == "ABOVE b_r"
    # exact interval arithmetic at the boundary
    assert C.q2_compare(hi - 0.01, lo + 0.01) == "NEAR b_r"


def test_q1_outcome_order():
    S = lambda st="SETTLED-ROBUST": dict(label="SETTLED", replay_ok=True, status=st)
    A = dict(label="ADVANCES", replay_ok=True, status=None)
    I = dict(label="INVADES", replay_ok=True, status=None)
    assert C.q1_outcome({1: S(), 2: S(), 3: S()}) == "BOUNDED"
    assert C.q1_outcome({1: A, 2: I, 3: A}) == "SLOW INVASION"
    assert C.q1_outcome({1: S(), 2: S(), 3: A}) == "MIXED"
    assert C.q1_outcome({1: S(), 2: S("NOT-ROBUST"), 3: S()}) == "MIXED"
    assert C.q1_outcome({1: S(), 2: S("ROBUSTNESS-UNRESOLVED"), 3: S()}) == "UNRESOLVED"
    assert C.q1_outcome({1: S(), 2: dict(label="UNRESOLVED-SLOW", replay_ok=True), 3: S()}) == "UNRESOLVED"
    assert C.q1_outcome({1: dict(S(), replay_ok=False), 2: S(), 3: S()}) == "UNRESOLVED"
    assert C.q1_outcome({1: dict(label="RECEDES", replay_ok=True, status=None), 2: S(), 3: S()}) == "MIXED"


def test_replay_check():
    old = track(expo, expo, 1800)
    new = track(expo, expo, 12000)
    assert C.replay_check(new, old)["ok"]
    bad = [dict(r) for r in new]
    bad[100]["xR"] += 2e-6
    r = C.replay_check(bad, old)
    assert not r["ok"] and r["max_abs_diff"] > 1e-6


def test_mg9_tracks_replay_against_themselves_and_states_exist():
    for b in C.Q1_B:
        rec = C.mg9_front_track(b)
        assert C.replay_check(rec["track"], rec["track"])["ok"]
        y, prov = C.d02_prepared(b)
        assert y.shape == (1500,) and prov == rec["prep"]
        assert abs(C.e_pat(y, 500, 0.2) - rec["E_pat"]) <= 1e-12
    y, prov = C.d02_prepared(C.DESC_B)
    assert prov == "sweep_29.22" and y.shape == (1500,)


class FakeSolver:
    """Stands in for mg9_ode: field stays constant; records the segments requested."""
    def __init__(self, residual=1e-7):
        self.calls, self.residual = [], residual

    def __call__(self, b, N, h, t_end, y0, *, seed_scale, t_eval=None, t_start=0.0, **kw):
        assert seed_scale == 0.0
        self.calls.append((t_start, t_end))
        te = np.asarray(t_eval if t_eval is not None else [t_end])
        Y = np.repeat(np.asarray(y0, float)[:, None], te.size, axis=1)
        return dict(t=te, y=Y, y_final=np.asarray(y0, float), residual=self.residual, success=True,
                    t_reached=t_end, complete=True, retry=None)


def _with_tmp(fn):
    old = (C.OUT, C.RUNS, C.SNAP, C.R.mg9_ode)
    with tempfile.TemporaryDirectory() as d:
        C.OUT, C.RUNS, C.SNAP = d, os.path.join(d, "runs"), os.path.join(d, "snap")
        try:
            fn(d)
        finally:
            C.OUT, C.RUNS, C.SNAP, C.R.mg9_ode = old


def test_segmentation_snapshots_and_ledger():
    def body(d):
        y, _ = C.d02_prepared(29.30)
        E = C.e_pat(y, 500, 0.2)
        fake = FakeSolver()
        C.R.mg9_ode = fake
        rec = C.long_front("t_main", 29.30, 500, 0.2, C.patch(y, 500, 0.2, 25.0), E)
        assert fake.calls == list(zip(C.MAIN_SEGS, C.MAIN_SEGS[1:]))
        assert [s["t"] for s in rec["snapshots"]] == [1800.0, 3600.0, 6000.0, 12000.0]
        assert rec["label"] == "SETTLED" and len(rec["track"]) == 2401
        assert len(C.ledger()["integrations"]) == 1
        fake2 = FakeSolver(residual=1e-3)                          # rule 6 -> one extension, then UNRESOLVED-SLOW
        C.R.mg9_ode = fake2
        rec = C.long_front("t_ext", 29.30, 500, 0.2, C.patch(y, 500, 0.2, 25.0), E)
        assert fake2.calls[-2:] == [(12000.0, 15000.0), (15000.0, 30000.0)]
        assert rec["extended"] and rec["label"] == "UNRESOLVED-SLOW"
        assert [s["t"] for s in rec["snapshots"]][-2:] == [15000.0, 30000.0]
        fake3 = FakeSolver(residual=1e-3)
        C.R.mg9_ode = fake3
        rec = C.long_front("t_pert", 29.30, 500, 0.2, C.patch(y, 500, 0.2, 25.0), E, kind="pert")
        assert fake3.calls == list(zip(C.PERT_SEGS, C.PERT_SEGS[1:]))
        assert rec["label"] == "UNRESOLVED-SLOW" and not rec["extended"]
        assert [s["t"] for s in rec["snapshots"]] == [1500.0, 3000.0, 6000.0]
        assert len(C.ledger()["integrations"]) == 3
        assert C.long_front("t_main", 29.30, 500, 0.2, None, E)["key"] == "t_main"   # resumable, no new count
        assert len(C.ledger()["integrations"]) == 3
    _with_tmp(body)


def test_cap_enforced():
    def body(d):
        with open(os.path.join(d, "ledger.json"), "w") as f:
            json.dump(dict(integrations=[dict(key=str(i), wall_s=1.0) for i in range(C.RUN_CAP)]), f)
        try:
            C._ledger_check()
        except C.Stop:
            return
        raise AssertionError("cap not enforced")
    _with_tmp(body)


if __name__ == "__main__":
    t0 = time.time()
    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]
    for name, fn in tests:
        t = time.time()
        fn()
        print(f"  ok  {name} ({time.time() - t:.1f} s)", flush=True)
    print(f"{len(tests)} tests passed in {time.time() - t0:.1f} s")
