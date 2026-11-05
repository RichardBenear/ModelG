"""MG9b tests (spec/MG9b_proposal.md draft 3; registered, frozen at 185aceb). Run: python tests/test_mg9b.py
No MG9b trajectory is simulated: the driver tests use a fake run generator and a fake confirmation. The confirmation
tests integrate MG003's b = 29 reference over the standard 200-300 window (about 10 s each; not MG9b data)."""
import os, sys, json, shutil, tempfile, math, zlib
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import mg9b
import confirm
import calib_mg9b as C
import run_mg9b as D
import evaluate_mg9b as EV

MG3B = [(130.0, 25, 15461.0), (120.0, 44, 12891.5), (110.0, 28, 4464.0), (100.0, 21, 1865.5)]


# ------------------------------------------------------------------ confirmation path (V-b2, V-b3)
def _ref_field():
    ref = np.load(os.path.join(D.ROOT, "results", "MG003", "ref_profiles.npz"))["h0.2_b29.0"]
    return ref[1:4]


def test_confirm_defaults_unchanged_and_tolerance_args():
    phi = _ref_field()
    a = confirm.confirm_field(phi, 29.0, 0.2, float(phi[1, 250]), 3.0)
    b = confirm.confirm_field(phi, 29.0, 0.2, float(phi[1, 250]), 3.0, checks=confirm.CHECKS, rtol=1e-8, atol=1e-10)
    assert a == b and a["complete"] and a["cls"] == "persistent"


def test_confirm_late_window_and_tolerances_are_passed():
    seen = {}
    old = confirm.lattice_ode_1d

    def fake(P, N, h, t_end, seed_scale=1.0, y0=None, rtol=1e-8, atol=1e-10, t_eval=None):
        seen.update(t_end=t_end, t_eval=np.array(t_eval), rtol=rtol, atol=atol)
        return None, t_eval, np.repeat(y0[:, None], len(t_eval), 1), 0.0
    confirm.lattice_ode_1d = fake
    try:
        phi = _ref_field()
        c = confirm.confirm_field(phi, 29.0, 0.2, float(phi[1, 250]), 3.0, checks=mg9b.VB2_LATE, **mg9b.VB3_TOL)
    finally:
        confirm.lattice_ode_1d = old
    assert seen["t_end"] == 2000.0 and seen["t_eval"][0] == 1900.0 and len(seen["t_eval"]) == 11
    assert seen["rtol"] == 1e-10 and seen["atol"] == 1e-12 and c["cls"] == "persistent"


def test_confirm_incomplete_is_never_a_class():
    old = confirm.lattice_ode_1d
    confirm.lattice_ode_1d = lambda *a, t_eval=None, **k: (None, t_eval[:4], np.zeros((1500, 4)), 0.0)
    try:
        c = confirm.confirm_field(_ref_field(), 29.0, 0.2, -1.0, 3.0)
    finally:
        confirm.lattice_ode_1d = old
    assert c["cls"] == "incomplete" and not c["complete"] and c["n_checks_returned"] == 4


# ------------------------------------------------------------------ fits
def test_fit_b_reproduces_mg3b():
    f = mg9b.fit_b(MG3B)
    assert abs(f["S"] - 0.06420908708993153) < 1e-12 and abs(f["deviance"] - 0.22729586691512227) < 1e-9
    assert abs(f["ci"][0] - 0.04651507938637795) < 1e-8 and abs(f["ci"][1] - 0.08167736659606697) < 1e-8
    assert f["df"] == 2 and abs(f["p_lof"] - 0.8925721453604557) < 1e-9


def test_lack_of_fit_df():
    assert mg9b.fit_b(MG3B[:3])["df"] == 1 and mg9b.fit_b(MG3B[:3])["p_lof"] is not None
    f2 = mg9b.fit_b(MG3B[:2])
    assert f2["df"] == 0 and f2["p_lof"] is None


def test_identifiability_rule():
    assert mg9b.identifiable([(150, 1, 100.0), (135, 1, 100.0), (120, 0, 50.0)])
    assert not mg9b.identifiable([(150, 0, 100.0), (135, 5, 100.0), (120, 0, 50.0)])
    assert not mg9b.identifiable([(150, 3, 100.0)])
    assert mg9b.fit_b([(150, 0, 100.0), (135, 5, 100.0)]) == dict(identifiable=False)
    f = mg9b.fit_b([(150, 0, 4000.0), (135, 2, 3000.0), (120, 9, 2000.0)])
    assert f["ok"] and math.isfinite(f["S"])


def test_profile_unbounded_side_is_infinite_not_failure():
    f = mg9b.fit_b([(150.0, 1, 10.0), (149.0, 1, 10.0)])
    assert f["ok"] and (math.isinf(f["ci"][0]) or math.isinf(f["ci"][1]))


def _exact(S, Om0, Oms, E=5000.0):
    return [(o, 3.41e-5 * math.exp(-S * (o - Om0)) * mg9b.L * E, E) for o in Oms]


def test_trend_recovers_linear_S1_and_lr_only_with_three_b():
    S = lambda b: 0.064 - 0.1 * (b - 29.0)
    by = {29.0: _exact(S(29.0), 120, (130, 120, 110, 100)), 29.2: _exact(S(29.2), 140, (150, 135, 120)),
          29.35: _exact(S(29.35), 160, (170, 155, 140))}
    t = mg9b.fit_trend(by)
    assert t["status"] == "ok" and abs(t["S1"] + 0.1) < 1e-6 and t["lr"]["df"] == 1 and t["lr"]["stat"] < 1e-6
    t2 = mg9b.fit_trend({29.0: by[29.0], 29.2: by[29.2]})
    assert "lr" not in t2 and mg9b.branch_of(t2["b_used"]) == "B2-20"
    assert mg9b.fit_trend({29.0: by[29.0]})["status"] == "insufficient"
    assert mg9b.branch_of([29.0, 29.35]) == "B2-35" and mg9b.branch_of([29.0, 29.2, 29.35]) == "B3"


def test_nonmonotone_case_flags_misfit():
    rng = np.random.default_rng(1)
    def cells(S, Om0, Oms):
        return [(o, int(rng.poisson(3.41e-5 * math.exp(-S * (o - Om0)) * mg9b.L * 2e5)), 2e5) for o in Oms]
    by = {29.0: cells(0.064, 120, (130, 120, 110, 100)), 29.2: cells(0.080, 140, (150, 135, 120)),
          29.35: cells(0.020, 160, (170, 155, 140))}
    t = mg9b.fit_trend(by)
    assert t["S1"] < 0 and t["lr"]["p"] < 0.05


def test_primary_rules():
    ok_neg = dict(status="ok", ci=(-0.2, -0.01))
    assert mg9b.primary(dict(status="insufficient"), True, False)["label"] == "NOT RESOLVED (insufficient)"
    assert mg9b.primary(dict(status="numerical"), True, False)["label"] == "NOT RESOLVED (NUMERICAL)"
    assert mg9b.primary(ok_neg, True, False)["confirmatory"] == "NEGATIVE LINEAR TREND"
    p = mg9b.primary(ok_neg, False, False)
    assert p["label"] == "NEGATIVE LINEAR TREND (EXPLORATORY)" and p["confirmatory"] == "NOT RESOLVED"
    p = mg9b.primary(ok_neg, True, True)
    assert p["exploratory"] and p["why"] == ["BUDGET-TRUNCATED"]
    assert mg9b.primary(dict(status="ok", ci=(0.01, math.inf)), True, False)["label"] == "POSITIVE LINEAR TREND"
    assert mg9b.primary(dict(status="ok", ci=(-math.inf, 0.1)), True, False)["label"] == "NOT RESOLVED"


# ------------------------------------------------------------------ design rules
def test_parta_rules():
    c = lambda Om, n, ok=True: dict(Om=Om, screen_ok=ok, count=n)
    assert mg9b.parta_next([])["next_Om"] == 250.0
    r = mg9b.parta_next([c(250, 0), c(200, 2), c(170, 1), c(150, 3)])
    assert r["stop"] and r["Om_T"] == 200
    r = mg9b.parta_next([c(250, 0), c(200, 1), c(170, None, False)])
    assert r["stop"] and r["Om_T"] is None and "screen failed" in r["reason"]
    r = mg9b.parta_next([c(Om, 0) for Om in mg9b.GRID])
    assert r["stop"] and r["Om_T"] is None and r["reason"] == "end of grid"


def test_partb_candidates_unique_and_ordered():
    assert mg9b.candidates(150.0) == [135.0, 120.0, 165.0, 180.0]
    acc, nxt = mg9b.partb_cells(150.0, [False, False])
    assert nxt == 165.0
    acc, nxt = mg9b.partb_cells(150.0, [False, False, True, True])
    assert acc == [150.0, 165.0, 180.0] and nxt is None and len(set(acc)) == 3
    acc, nxt = mg9b.partb_cells(150.0, [True, True])
    assert acc == [150.0, 135.0, 120.0] and nxt is None
    acc, nxt = mg9b.partb_cells(150.0, [False, False, False, True])
    assert acc == [150.0, 180.0] and nxt is None


def test_stop20():
    assert mg9b.stop20_need([1] * 10) == (False, 16)
    assert mg9b.stop20_need([1] * 16 + [1] * 4 + [0] * 12) == (True, 32)
    assert mg9b.stop20_need([1] * 20 + [0] * 12) == (True, 32)
    assert mg9b.stop20_need([1] * 16) == (False, 32)
    assert mg9b.stop20_need([0] * 64) == (True, 64)


# ------------------------------------------------------------------ calibration generator
def test_calibration_generator():
    rng = np.random.default_rng(0)
    ev, ex = C.runs(rng, 1.0, 1.0, 50)
    assert ev.sum() == 0 and ex.sum() == 0.0                              # startup: no event, zero exposure
    kept = sum(C.historical(np.random.default_rng(i))[0] for i in range(200))
    assert 140 <= kept <= 195                                              # P(EQUIVALENT | ratio 1) ~ 0.85
    cs = C.cases()
    assert len(cs) == 9 and sum(c["S1"] is not None for c in cs) == 8
    for i, c in enumerate(cs):
        c["idx"] = i
    r = C.replicate(cs[0], 3)
    assert r == C.replicate(cs[0], 3)                                      # reproducible


def test_calibration_labels_are_per_branch():
    def case(S1, n, cov, neg=0.0, pos=0.0):
        b = dict(n=n, coverage=dict(p=cov), neg=dict(p=neg), pos=dict(p=pos))
        return dict(case="x", S1_true=S1, branches={"B3": b, "B2-20": dict(b, n=50), "B2-35": dict(b, n=50)})
    lab = C.labels([case(0.0, 1000, 0.95), case(-0.1, 1000, 0.94)])
    assert lab["B3"]["label"] == "calibrated" and lab["B2-20"]["label"] == "not calibrated"
    lab = C.labels([case(0.0, 1000, 0.95, neg=0.04), case(-0.1, 1000, 0.94)])
    assert lab["B3"]["label"] == "not calibrated"


# ------------------------------------------------------------------ driver with a fake simulator
class SyncPool:
    """In-process pool: fake runs never reach real worker processes (spawn on Windows would re-import the driver
    with the real simulate() and the real output folder)."""
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def imap(self, f, it, chunksize=1):
        return map(f, it)

    def imap_unordered(self, f, it, chunksize=1):
        return map(f, it)


REAL_OUT = os.path.join(D.ROOT, "results", "MG9b")


def _real_tree():
    return sorted((os.path.relpath(os.path.join(d, f), REAL_OUT), os.path.getmtime(os.path.join(d, f)))
                  for d, _, fs in os.walk(REAL_OUT) for f in fs)


class Tmp:
    def __enter__(self):
        self.tmp = tempfile.mkdtemp()
        self.real = _real_tree()
        self.old = (D.OUT, EV.OUT, D.simulate, D._confirm_job, D.vb1, D.part_a, D.production, D.vb2, D.POOL_FACTORY)
        D.OUT = EV.OUT = self.tmp
        D.POOL_FACTORY = SyncPool
        return self

    def __exit__(self, *a):
        D.OUT, EV.OUT, D.simulate, D._confirm_job, D.vb1, D.part_a, D.production, D.vb2, D.POOL_FACTORY = self.old
        shutil.rmtree(self.tmp)
        assert _real_tree() == self.real, "a test touched the real results/MG9b folder"


def fake_simulate(a):
    """Deterministic fake run: event probability falls with Omega (Omega0 140 at 29.20, 160 at 29.35)."""
    Om0 = 140.0 if a["b"] < 29.3 else 160.0
    rng = np.random.default_rng(zlib.crc32(f"{a['cell']}/{a['rep']}".encode()))
    lam = 3.41e-5 * math.exp(-0.05 * (a["Om"] - Om0)) * 100
    st = rng.random() < 0.2
    u = rng.exponential() / lam
    t_first = 30.0 if st else (50.0 + u if u <= 430 else None)
    if a["until"] < 500:
        t_first = None if (t_first is None or t_first > a["until"]) else t_first
    event = (not st) and t_first is not None
    conf = [dict(sc=7, t_conf=(t_first or 0) + 20, cls="persistent", n_persistent=1, n_deep_end=1,
                 min_phiX_end=-1.0, maxdt_end=0.0, complete=True)] if t_first is not None else []
    rec = dict(key=list(a["key"]), b=a["b"], b_det=a["b"], Om=a["Om"], tau=mg9b.TAU0, t_end=a["until"], wall=1.0,
               clip_frac=0.0, mean_dets_t_lt_50=0.1, startup=st, event=event,
               exposure=0.0 if st else (min(u, 430.0) if a["until"] >= 500 else 0.0), complete=True,
               t_first=t_first, first=dict(s_conf=7, t_start=t_first) if event else None, births=[],
               confirmations=conf)
    return rec, ({"7": np.zeros((3, 500), np.float32)} if conf else {})


def fake_confirm(a):
    return dict(cell=a["cell"], rep=a["rep"], sc=a["sc"], cls="persistent", complete=True, n_persistent=1, wall=0.5)


def test_vb1_fresh_folder_with_static_reference():
    """V-b1 on a fresh output folder with a fake solver returning the reference unchanged: PASS, files written."""
    old = D.lattice_ode_1d
    def fake(P_, N, h, t_end, seed_scale=1.0, y0=None, t_eval=None, **k):
        if t_eval is None:
            return None, None, None, None, 0.0
        return None, t_eval, np.repeat(y0[:, None], len(t_eval), 1), 0.0
    with Tmp() as t:
        D.OUT = os.path.join(t.tmp, "fresh")
        D.lattice_ode_1d = fake
        try:
            r = D.vb1(29.35)
        finally:
            D.lattice_ode_1d = old
        assert r["result"] == "PASS" and r["drift"] < 1e-12
        assert os.path.exists(os.path.join(D.OUT, "vb1_29.35_final.npz"))


def test_driver_guard_requires_go():
    with Tmp():
        try:
            D.stage_runs()
            assert False
        except SystemExit as e:
            assert e.code == 3


def test_driver_end_to_end_fake_and_evaluate():
    with Tmp() as t:
        D.simulate = fake_simulate
        D._confirm_job = fake_confirm
        D.vb1 = lambda b: dict(b=b, pass_=True, result="PASS")
        open(os.path.join(t.tmp, "GO"), "w").write("test")
        D.stage_runs()
        pa = json.load(open(os.path.join(t.tmp, "partA_29.35.json")))
        pb = json.load(open(os.path.join(t.tmp, "partB_29.35.json")))
        assert pa["Om_T"] is not None and len(pb["accepted"]) == 3
        assert all(c["runs"] % 16 == 0 and (c["events"] >= 20 or c["runs"] == 64) for c in pb["cells"])
        v2 = json.load(open(os.path.join(t.tmp, "vb2_29.35.json")))
        assert v2["result"] == "PASS" and v2["n_agree"] == 10
        lg = json.load(open(os.path.join(t.tmp, "ledger.json")))
        assert lg["budget_stop"] is None and len(lg["entries"]) > 50
        ev = EV.evaluate()
        assert ev["branch"] == "B3" and ev["primary"]["label"] is not None
        assert ev["per_b"]["29.35"]["identifiable"] and "landmarks" in next(iter(ev["per_b"]["29.35"]["cells"].values()))
        # no calibration file in the temp folder: any directional result must be exploratory
        assert ev["primary"]["confirmatory"] == "NOT RESOLVED"


def test_vb2_nine_agree_plus_one_incomplete_fails():
    with Tmp() as t:
        D.simulate = fake_simulate
        calls = {"n": 0}

        def cj(a):
            calls["n"] += 1
            inc = calls["n"] == 10
            return dict(cell=a["cell"], rep=a["rep"], sc=a["sc"], cls="incomplete" if inc else "persistent",
                        complete=not inc, n_persistent=0 if inc else 1, wall=0.5)
        D._confirm_job = cj
        for r in range(12):
            D.run_job(dict(b=29.35, Om=150.0, cell="b29.35_B_150", rep=r, item="B", until=500.0,
                           key=(9, 2, 150, 3, r)))
        pa = dict(cells=[])
        pb = dict(accepted=[150.0])
        assert len(D.confirmation_order(29.35, pa, pb)) >= 10
        v = D.vb2(None, 29.35, pa, pb)
        assert v["n_agree"] == 9 and v["result"] == "CONTROL-NUMERICAL-FAILURE"


def test_budget_stop_marks_cell_and_halts():
    with Tmp() as t:
        D.simulate = fake_simulate
        D._confirm_job = fake_confirm
        D.vb1 = lambda b: dict(b=b, pass_=True, result="PASS")
        open(os.path.join(t.tmp, "GO"), "w").write("test")
        old = D.mean_run_s
        D.mean_run_s = lambda: 15 * 3600.0           # 15 core-h per run: the cutoff is hit at the first Part A set
        try:
            D.stage_runs()
        finally:
            D.mean_run_s = old
        lg = json.load(open(os.path.join(t.tmp, "ledger.json")))
        assert lg["budget_stop"] is not None
        try:
            D.stage_runs()
            assert False
        except SystemExit as e:
            assert e.code == 2


# ------------------------------------------------------------------ amendment 1: study-level budget truncation
def _fake_setup(t):
    D.simulate = fake_simulate
    D._confirm_job = fake_confirm
    D.vb1 = lambda b: dict(b=b, pass_=True, result="PASS")
    open(os.path.join(t.tmp, "GO"), "w").write("test")
    os.makedirs(os.path.join(t.tmp, "calibration"))
    json.dump(dict(labels={br: dict(label="calibrated") for br in ("B3", "B2-20", "B2-35")}),
              open(os.path.join(t.tmp, "calibration", "calibration.json"), "w"))


def test_budget_before_29_20_starts_is_study_truncation():
    with Tmp() as t:
        _fake_setup(t)
        real = D.part_a
        def pa(pool, b):
            if b < 29.3:
                raise D.BudgetStop("Part A 29.20")
            return real(pool, b)
        D.part_a = pa
        D.stage_runs()
        lg = json.load(open(os.path.join(t.tmp, "ledger.json")))
        assert lg["budget_stop"]["rate_complete"] is False and lg["budget_stop"]["budget_truncated"]
        ev = EV.evaluate()
        assert ev["branch"] == "B2-35" and ev["budget_truncated"]
        assert not ev["per_b"]["29.35"]["budget_short"]                     # no short cell anywhere
        assert ev["primary"]["confirmatory"] == "NOT RESOLVED" and "BUDGET-TRUNCATED (study)" in ev["qualifiers"]
        if ev["primary"]["label"].startswith(("NEGATIVE", "POSITIVE")):
            assert ev["primary"]["exploratory"] and "BUDGET-TRUNCATED" in ev["primary"]["why"]


def test_budget_before_third_cell_is_study_truncation():
    with Tmp() as t:
        _fake_setup(t)
        real = D.production
        seen = {}
        def prod(pool, b, Om):
            seen[b] = seen.get(b, 0) + 1
            if b < 29.3 and seen[b] == 3:
                raise D.BudgetStop("third cell")
            return real(pool, b, Om)
        D.production = prod
        D.stage_runs()
        ev = EV.evaluate()
        d = ev["per_b"]["29.20"]
        done = [c for c in d["cells"].values() if c["runs"] > 0]
        assert len(done) == 2 and not any(c["budget_short"] for c in done) and d["identifiable"]
        assert len(d["unstarted_cells"]) == 1 and not d["budget_short"]      # unstarted, not short
        assert ev["budget_truncated"] and ev["primary"]["confirmatory"] == "NOT RESOLVED"


def test_failed_reference_then_control_stop_is_not_rate_truncation():
    """29.20 fails V-b1; 29.35 finishes Part B; the budget stops its V-b3: controls only."""
    with Tmp() as t:
        _fake_setup(t)
        D.vb1 = lambda b: dict(b=b, pass_=b > 29.3, result="PASS" if b > 29.3 else "REFERENCE-INVALID")
        old3 = D.vb3
        def v3(pool, b, Om_T):
            raise D.BudgetStop("V-b3 29.35")
        D.vb3 = v3
        try:
            D.stage_runs()
        finally:
            D.vb3 = old3
        lg = json.load(open(os.path.join(t.tmp, "ledger.json")))
        assert lg["budget_stop"]["rate_complete"] is True and lg["budget_stop"]["budget_truncated"] is False
        ev = EV.evaluate()
        assert not ev["budget_truncated"] and "BUDGET STOP AFTER RATE COLLECTION (controls only)" in ev["qualifiers"]


def test_event_shapes_invariant_to_cell_order():
    def run(cell, nd):
        return dict(cell=cell, rep=0, event=True, first=dict(s_conf=80),
                    confirmations=[dict(sc=80, n_deep_end=nd, n_persistent=nd)])
    a, b = run("b29.35_B_150", 1), run("b29.35_B_135", 2)
    for order in ([a, b], [b, a]):
        e = EV.event_shapes(order)
        assert e["events"] == 2 and e["single_minimum"] == 1


def test_vb3_with_no_confirmations_is_insufficient():
    with Tmp() as t:
        os.makedirs(os.path.join(t.tmp, "runs", "b29.35_B_150"))
        json.dump(dict(rep=0, confirmations=[], event=False, complete=True),
                  open(os.path.join(t.tmp, "runs", "b29.35_B_150", "0.json"), "w"))
        v = D.vb3(None, 29.35, 150.0)
        assert v["n"] == 0 and v["result"] == "CONTROL-INSUFFICIENT"


def test_budget_stop_during_last_control_is_not_rate_truncation():
    with Tmp() as t:
        _fake_setup(t)
        real = D.vb2
        def v2(pool, b, pa, pb):
            if b < 29.3:
                raise D.BudgetStop("V-b2 29.20")
            return real(pool, b, pa, pb)
        D.vb2 = v2
        D.stage_runs()
        ev = EV.evaluate()
        assert ev["budget_stop"]["rate_complete"] is True and not ev["budget_truncated"]
        assert ev["per_b"]["29.20"]["vb2"] is None
        assert any("CONTROL-INSUFFICIENT" in f for f in ev["per_b"]["29.20"]["lattice_flags"])
        assert "BUDGET STOP AFTER RATE COLLECTION (controls only)" in ev["qualifiers"]
        if ev["primary"]["label"].startswith(("NEGATIVE", "POSITIVE")):
            assert not ev["primary"]["exploratory"]


def test_branch_catch_all_without_historical():
    assert mg9b.branch_of([29.2, 29.35]) == "B2-new"
    f = mg9b.fit_b([(150.0, 1, 10.0), (149.0, 1, 10.0)])
    assert f["ci_note"] and "not located within the registered search range" in f["ci_note"]


if __name__ == "__main__":
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_")]
    for k, f in tests:
        f()
        print("ok", k)
    print(len(tests), "tests passed")
