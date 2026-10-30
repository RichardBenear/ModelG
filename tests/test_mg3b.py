"""Tests for src/mg3b.py and src/run_mg3b.py (MG3b driver logic). Run: python tests/test_mg3b.py (~1 min)."""
import os, sys, json, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import mg3, mg3b


def fake(rep, event=False, startup=False, exposure=430.0, wall=100.0):
    return dict(rep=rep, event=event, startup=startup, complete=True, exposure=0.0 if startup else exposure,
                expo_w1=0.0, expo_w2=0.0, in_w1=False, wall=wall)


def test_stop20_rule():
    # 16 runs with 20 events is impossible (one first event per run): 2 batches with 12 + 9 = 21 events -> stop at 32
    runs = [fake(i, event=(i < 12 or 16 <= i < 25)) for i in range(32)]
    assert mg3b.stop20_done(runs) == (True, 32)
    assert mg3b.stop20_done(runs[:16]) == (False, 32)
    # missing replicate inside batch 1 -> batch 1 is not complete
    assert mg3b.stop20_done([r for r in runs if r["rep"] != 3]) == (False, 16)
    # cap: 64 runs, few events -> done at 64
    few = [fake(i, event=(i % 10 == 0)) for i in range(64)]
    assert mg3b.stop20_done(few) == (True, 64)
    assert mg3b.stop20_done(few[:48]) == (False, 64)
    assert mg3b.stop20_done([]) == (False, 16)


def test_compare_and_branch():
    a = [fake(i, event=(i < 30), exposure=200.0) for i in range(64)]       # 30 events
    b = [fake(i, event=(i < 30), exposure=200.0) for i in range(64)]
    out, c = mg3b.compare("a", a, b)
    assert c["n_a"] == 30 and c["n_b"] == 30
    assert out in ("EQUIVALENT", "INCONCLUSIVE")
    big_a = [fake(i, event=True, exposure=10.0) for i in range(400)]
    big_b = [fake(i, event=True, exposure=10.0) for i in range(400)]
    assert mg3b.compare("a", big_a, big_b)[0] == "EQUIVALENT"
    tri = [fake(i, event=True, exposure=30.0) for i in range(400)]        # 3x lower rate -> DIFFERENT
    assert mg3b.compare("a", big_a, tri)[0] == "DIFFERENT"
    assert mg3b.compare(None, big_a, big_b)[0] == "INCONCLUSIVE"         # no eligible method
    assert mg3b.production_branch("EQUIVALENT") == ("tau", mg3b.TAU0, "P1_tau")
    assert mg3b.production_branch("INCONCLUSIVE") == ("tau2_inc", mg3b.TAU0 / 2, "P1_tau2")
    assert mg3b.production_branch("DIFFERENT")[0] == "part1b_needed"
    assert mg3b.production_branch("DIFFERENT", "EQUIVALENT") == ("tau2_1b", mg3b.TAU0 / 2, "P1_tau2")
    assert mg3b.production_branch("DIFFERENT", "INCONCLUSIVE")[0] == "stop"
    assert mg3b.production_branch("DIFFERENT", "DIFFERENT")[0] == "stop"


def test_budget_check():
    # prior projection: P1_tau 16 runs x 213.75 s x 1.25 = 1.1875 h; ceiling 6 h
    ok, info = mg3b.batch_allowed(0.0, [], 16, "P1_tau")
    assert ok and abs(info["projected_h"] - 16 * 3.8 * 3600 / 64 * 1.25 / 3600) < 1e-9
    ok, _ = mg3b.batch_allowed(5.0 * 3600, [], 16, "P1_tau")              # 5 + 1.19 > 6
    assert not ok
    ok, _ = mg3b.batch_allowed(5.0 * 3600, [], 16, "P1_tau", n_arm=96)     # ceiling 9 with 96 runs
    assert ok
    # completed-run mean replaces the prior
    ok, info = mg3b.batch_allowed(0.0, [3600.0], 16, "P2_tau")             # 16 x 1 h x 1.25 = 20 h > 12
    assert not ok and abs(info["projected_h"] - 20.0) < 1e-9
    # tau/4 ceiling: unscaled 40 h (PI decision); the measured cost still scales the prior projection
    assert mg3b.ceiling_h("P1_tau4", tau4_factor=0.54) == 40.0             # PI decision 2026-10-30: unscaled
    assert mg3b.ceiling_h("P2_tau2") == 30.0


def test_ledger_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "ledger.csv")
        mg3b.ledger_append(p, "P1_tau", "P1_tau/0", "run", 120.0, 9)
        mg3b.ledger_append(p, "P1_tau", "P1_tau/1", "run", 80.0, 9)
        mg3b.ledger_append(p, "timing", "TIMING_tau4/0", "timing", 50.0, 9)
        assert mg3b.ledger_used_s(p) == 250.0
        assert mg3b.ledger_used_s(p, "P1_tau") == 200.0
        rows = mg3b.ledger_rows(p)
        assert rows[0]["workers"] == "9" and rows[0]["numba"] != "" and len(rows) == 3


def test_summarize_b_keeps_every_confirmation_field():
    st = mg3.new_state(29.0, 100.0, mg3b.TAU0, (mg3b.GATE, 0, 3, 99, 0))   # test-only stage index 99
    st = mg3.sim_run(st, until=10.0)
    st["conf_cache"][5] = dict(cls="excursion", n_persistent=0, n_deep_end=0, min_phiX_end=0.0, maxdt_end=0.1)
    st["conf_cache"][12] = dict(cls="unresolved", n_persistent=0, n_deep_end=1, min_phiX_end=-0.1, maxdt_end=0.1)
    rec = mg3b.summarize_b(st)
    assert set(rec["conf_fields"]) == {5, 12}
    f = rec["conf_fields"][5]
    assert f["field"].shape == (3, mg3.N) and f["field"].dtype == np.float32 and f["cls"] == "excursion"
    assert np.array_equal(f["field"], np.asarray(st["phi"][5], np.float32))
    assert abs(f["t_conf"] - st["times"][5]) < 1e-12


def test_driver_refuses_without_decisions():
    with tempfile.TemporaryDirectory() as d:
        os.environ["MG3B_OUT"] = d                       # keep the test's log out of results/MG3b
        import run_mg3b
        try:
            run_mg3b.load_decisions(os.path.join(d, "missing.json"))
        except SystemExit as e:
            assert e.code == 5
        else:
            raise AssertionError("expected SystemExit")
    dec = run_mg3b.load_decisions(os.path.join(os.path.dirname(__file__), "..", "results", "MG3b", "calibration",
                                               "decisions.json"))
    assert dec["part1_runs_per_arm"] == 64 and dec["rate_ratio_ci"]["method"] == "a"


def test_stream_keys_distinct():
    from runner import seed_for
    keys = [(mg3b.GATE, 0, mg3b.om_idx(Om), mg3b.STAGE[s], r) for Om in mg3b.OMEGAS
            for s in ("P1_tau", "P1_tau2", "P1_tau4", "P2", "TIMING", "PILOT") for r in range(3)]
    seeds = {seed_for(*k) for k in keys}
    assert len(seeds) == len(keys)
    # MG3 used gate 3: no overlap with MG3's Stage-B keys
    assert seed_for(4, 0, 1, 0, 0) != seed_for(3, 0, 1, 0, 0)


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t(); print("ok", t.__name__)
    print(f"{len(tests)} passed")
