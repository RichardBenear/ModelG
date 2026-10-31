"""MG4c driver tests (proposal §7): the Part 2 rules as pure functions of constructed class tables, the threshold
eligibility rule, the high-precision balance recheck, and a small end-to-end smoke run (no registered graph is run).
    python tests/test_mg4c.py
"""
import os, sys, json, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import run_mg4c as rm

S, M, Dc, U = "single", "multiplied", "decayed", "unresolved"
BS = [29.0] + rm.GRID


def tab(rows):
    """rows: {b: list of 4 classes}; unspecified b are all multiplied (b = 29 is multiplied as in MG4)."""
    return {b: dict(D=rows.get(b, [M] * 4)) for b in BS}


def fake_run(truth):
    """truth(b) -> list of 4 classes; returns run_b and the list of visited b."""
    seen = []

    def run_b(b):
        seen.append(b)
        return dict(D=truth(b))
    return run_b, seen


def test_plan_case_A_interior_block():
    t = tab({28.25: [S] * 4, 28.0: [S] * 4, 27.75: [S] * 4, 27.0: [S] * 4})      # two blocks; the higher one is refined
    case, edges = rm.plan(t)
    assert case == "A" and edges == [(28.25, 28.5), (27.75, 27.5)]


def test_plan_case_A_unbracketed_lower():
    t = tab({b: [S] * 4 for b in rm.GRID[3:]})                                   # 28.0 ... 27.0 all single
    case, edges = rm.plan(t)
    assert case == "A" and edges == [(28.0, 28.25), "lower edge < 27.0"]


def test_plan_unresolved_is_not_single():
    t = tab({28.0: [S, S, U, S], 27.75: [S] * 4})
    case, edges = rm.plan(t)
    assert case == "A" and edges == [(27.75, 28.0), (27.75, 27.5)]


def test_plan_case_B_bracket_and_none():
    t = tab({28.0: [S, S, Dc, S], 27.75: [Dc] * 4})                              # no multiplied from 28.0 down
    assert rm.plan(t) == ("B", (28.0, 28.25))
    assert rm.plan(tab({})) == ("B", None)                                      # multiplied everywhere: no refinement
    t = tab({28.75: [Dc, S, S, S]})
    assert rm.plan(t) == ("B", (28.75, 29.0))


def test_bisect_case_A_edge():
    edge = 28.3                                                                   # truth: all single for b < 28.3
    run_b, seen = fake_run(lambda b: [S] * 4 if b < edge else [M] * 4)
    r = rm.bisect_edge(28.25, 28.5, rm.all_single, run_b)
    lo, hi = r["bracket"]
    assert len(seen) == 6 and lo < edge <= hi and abs(hi - lo - 0.25 / 64) < 1e-12
    # lower edge, bisected downward
    run_b, seen = fake_run(lambda b: [S] * 4 if b > 27.6 else [Dc] * 4)
    r = rm.bisect_edge(27.75, 27.5, rm.all_single, run_b)
    lo, hi = r["bracket"]
    assert lo <= 27.6 < hi and len(seen) == 6


def test_bisect_case_B_stops_on_all_single():
    # multiplied above 28.1, a narrow all-single band [28.05, 28.1), partial below
    truth = lambda b: [M] * 4 if b >= 28.1 else ([S] * 4 if b >= 28.05 else [S, Dc, S, S])
    run_b, seen = fake_run(truth)
    r = rm.bisect_edge(28.0, 28.25, lambda t: not rm.any_mult(t), run_b, stop_on=rm.all_single)
    assert r["bracket"] is None and 28.05 <= r["stopped_at"] < 28.1 and len(seen) < 6
    # no all-single point: runs the full 6 steps and returns the bracket of the "any multiplied" edge
    truth = lambda b: [M, S, S, S] if b >= 28.1 else [S, Dc, S, S]
    run_b, seen = fake_run(truth)
    r = rm.bisect_edge(28.0, 28.25, lambda t: not rm.any_mult(t), run_b, stop_on=rm.all_single)
    lo, hi = r["bracket"]
    assert len(seen) == 6 and lo < 28.1 <= hi


def test_threshold_eligibility():
    orig = rm.classify_run
    try:
        rm.classify_run = lambda rec, th: rec["cls"](th)
        T = lambda m, f: dict(min_smoothed_final=m, cls=f)
        r = rm.threshold_at(T(-8.0, lambda th: S))
        assert r["eligible"] and r["theta"] == -4.0
        r = rm.threshold_at(T(-1.5, lambda th: S))                               # too shallow: (ii) fails
        assert not r["eligible"] and r["theta"] == rm.THETA_29
        r = rm.threshold_at(T(-8.0, lambda th: Dc if th == rm.THETA_ELIG else S)) # (i) fails
        assert not r["eligible"] and r["theta"] == rm.THETA_29
        r = rm.threshold_at(T(-8.0, lambda th: S if th == rm.THETA_ELIG else M))  # T not single under theta_b
        assert not r["eligible"] and r["theta"] == rm.THETA_29
        assert rm.threshold_at(None)["theta"] == rm.THETA_29
    finally:
        rm.classify_run = orig


def test_keys():
    assert rm.job(28.75, "D", 1, 0)["key"] == "b28.750000_D1_p0"
    assert rm.job(29.0, "T", 0, 0, "cotan")["key"] == "w_b29.000000_T_p0"
    assert [j["key"] for j in rm.jobs_at(27.0)] == ["b27.000000_T_p0", "b27.000000_D0_p0", "b27.000000_D0_p1",
                                                    "b27.000000_D1_p0", "b27.000000_D1_p1"]


def test_smoke_small_graph_and_hp_recheck():
    """One short pipeline pass on a small T torus and a small Delaunay graph (not registered graphs): run_one, the
    record, the field file, the ledger, classification, and the high-precision recheck (unweighted and cotangent)."""
    import graphs as gr
    tmp = tempfile.mkdtemp()
    old = (rm.OUT, rm.graph_for)
    try:
        rm.OUT = tmp
        for d in ("runs", "fields"):
            os.makedirs(os.path.join(tmp, d))
        tor = gr.torus(0.2, 30, 34)
        gT = gr.triangular(tor)
        gD = gr.delaunay(gr.poisson_points(gT["N"], tor["Lx"], tor["Ly"], np.random.default_rng(5)), tor["Lx"], tor["Ly"])
        hT = float(np.sqrt(1.5) * tor["a"])
        hD = float(np.sqrt(gr.corrector(gD)["h2"]))
        c = np.array([tor["Lx"] / 2, tor["Ly"] / 2])
        rm.graph_for = lambda j, dec: (gT, hT, c) if j["kind"] == "T" else (gD, hD, c)
        recs = {}
        for j in (rm.job(29.0, "T"), rm.job(29.0, "D", 0, 0, "cotan")):
            r = rm.run_one(j)
            recs[j["key"]] = r
            assert os.path.exists(os.path.join(tmp, "fields", j["key"] + ".npz"))
            assert rm.classify_run(r, rm.THETA_29) in (S, M, Dc, U)
            hp = rm.hp_recheck(r)
            assert abs(hp["r1"]) <= rm.HP_TOL * (1 + hp["M1"]) and abs(hp["r2"]) <= rm.HP_TOL * (1 + hp["M2"]), hp
        ch, n = rm.spent()
        assert n == 2 and ch > 0
        # a forced flag is resolved by the recheck and recorded before any class is read
        k = "b29.000000_T_p0"
        recs[k]["balances"]["ok"] = False
        assert rm.resolve_flags(recs) == []
        assert json.load(open(os.path.join(tmp, "resolution.json")))[k]["resolved"]
    finally:
        rm.OUT, rm.graph_for = old


def test_evaluator_statuses():
    import evaluate_mg4c as ev
    t = lambda *rows: {b: dict(D=r) for b, r in zip(BS, rows)}
    allM = [[M] * 4] * 9
    assert ev.window_status(t(*allM), True, []) == "W2"
    assert ev.window_status(t(*allM), False, []) == "INCOMPLETE"
    assert ev.window_status(t(*allM), True, ["x"]) == "NOT SCORED"
    assert ev.window_status(t(*([[M] * 4] * 4 + [[S, S, S, S]] + [[M] * 4] * 4)), True, []) == "W1"
    assert ev.window_status(t(*([[M] * 4] * 4 + [[S, S, U, S]] + [[Dc] * 4] * 4)), True, []) == "W3"
    assert ev.window_status(t(*([[M] * 4] * 8 + [[Dc] * 4])), True, []) == "W2"   # decayed is not single
    ok = dict(T_control_maxdiff_X=1e-9)
    P = lambda d, T=S, **kw: {"w_b29.000000_D0_p0": d, "w_b29.000000_T_p0": T, **ok, **kw}
    assert ev.operator_status(True, P(S)) == "O1"
    assert ev.operator_status(True, P(M)) == "O2" and ev.operator_status(True, P(Dc)) == "O2"
    assert ev.operator_status(True, P(U)) == "O3"
    assert ev.operator_status(False, P(S)) == "O-invalid"
    assert ev.operator_status(True, P(S, T=M)) == "O-invalid"
    assert ev.operator_status(True, {**P(S), "T_control_maxdiff_X": 2e-6}) == "O-invalid"


def _synthetic_pipeline(truth, tmp):
    """Runs rm.main() and the evaluator on synthetic runs: a small T torus stands in for every graph, and each run's
    smoothed fields are wells whose number at t = 300 is truth(b, kind, slot, pl) (1 single, 2 multiplied, 0 decayed)."""
    import graphs as gr
    import evaluate_mg4c as ev
    tor = gr.torus(0.2, 62, 70)
    g = gr.triangular(tor)
    N = g["N"]
    c1, c2 = np.array([3.0, 3.0]), np.array([7.5, 6.0])
    well = lambda c: -8.0 * np.exp(-np.hypot(*gr.min_image(g["pos"] - c, g["Lx"], g["Ly"]).T) ** 2 / 0.5)

    def run_one(j):
        n = truth(j["b"], j["kind"], j["slot"], j["pl"])
        f = {0: np.zeros(N), 1: well(c1), 2: well(c1) + well(c2)}[n]
        sX = np.array([well(c1)] * 10 + [f])
        np.savez_compressed(os.path.join(rm.OUT, "fields", j["key"] + ".npz"), final=np.stack([f, f, f]),
                            checks_sX=sX.astype(np.float32))
        rec = dict(j, N=N, h=0.2, centre=c1.tolist(), wall_s=1.0, nfev=0, max_dphidt=0.0, stationary=True,
                   min_smoothed_final=float(sX[-1].min()), balances=dict(ok=True, r1=0, M1=0, r2=0, M2=0))
        json.dump(rec, open(os.path.join(rm.OUT, "runs", j["key"] + ".json"), "w"))
        rm.ledger(j["key"], 1.0)
        return rec

    np.save(os.path.join(tmp, "Tstate.npy"), np.stack([well(c1)] * 3).ravel())
    old = (rm.OUT, rm.graph_for, rm.run_one, rm.TSTATE, rm.part3_checks, ev.OUT)
    try:
        rm.OUT = ev.OUT = tmp
        rm.TSTATE = os.path.join(tmp, "Tstate.npy")
        rm.graph_for = lambda j, dec: (g, 0.2, c1)
        rm.run_one = run_one
        rm.part3_checks = lambda dec: dict(ok=True)
        rm._graphs.clear()
        log = rm.main()
        return log, ev.main()
    finally:
        rm.OUT, rm.graph_for, rm.run_one, rm.TSTATE, rm.part3_checks, ev.OUT = old
        rm._graphs.clear()


def test_pipeline_synthetic_case_A():
    if os.name == "nt":                     # the process pool needs fork to see the synthetic run_one
        return
    tmp = tempfile.mkdtemp()
    truth = lambda b, kind, slot, pl: 1 if kind == "T" or 27.6 <= b < 28.3 else 2
    log, ev = _synthetic_pipeline(truth, tmp)
    assert log["status"] == "complete" and log["case"] == "A"
    up, lo = log["edges"]
    assert up["bracket"][0] < 28.3 <= up["bracket"][1] and lo["bracket"][0] < 27.6 <= lo["bracket"][1]
    assert ev["window_status"] == "W1" and ev["operator_status"] == "O2"     # b = 29 multiplies in this truth
    assert log["runs"] == 40 + 2 + 2 * 6 * 5


def test_pipeline_synthetic_case_B_and_W2():
    if os.name == "nt":
        return
    tmp = tempfile.mkdtemp()
    truth = lambda b, kind, slot, pl: 1 if kind == "T" else (2 if b >= 28.1 else (1 if (slot, pl) != (1, 1) else 0))
    log, ev = _synthetic_pipeline(truth, tmp)
    assert log["case"] == "B" and log["case_B"]["bracket"][0] < 28.1 <= log["case_B"]["bracket"][1]
    assert ev["window_status"] == "W3"
    tmp = tempfile.mkdtemp()
    log, ev = _synthetic_pipeline(lambda b, kind, slot, pl: 1 if kind == "T" else 2, tmp)
    assert log["case_B"] == "multiplied at every tested b" and ev["window_status"] == "W2" and log["runs"] == 42


if __name__ == "__main__":
    import time
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            t0 = time.time(); f(); print(f"{name}: ok ({time.time() - t0:.1f}s)")
