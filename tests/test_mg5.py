"""MG5 tests (proposal §9). Small unregistered graphs and synthetic fields only; no registered graph is run.
    python tests/test_mg5.py
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import evolve as ev
import graphs as gr
import graph_ode as go
import measure_graph as mg
from lattice import P17

Lx, Ly = 6.0, 5.0


def pts(n, seed):
    return np.random.default_rng(seed).uniform((0, 0), (Lx, Ly), size=(n, 2))


def test_mesh_matches_graphs_delaunay():
    p = pts(400, 1)
    m = ev.Mesh(p, Lx, Ly)
    g = gr.delaunay(p, Lx, Ly)
    assert np.array_equal(m.g["I"], g["I"]) and np.array_equal(m.g["J"], g["J"])


def test_interpolation_exact_for_linear_fields_and_matches_independent():
    p = pts(500, 2)
    m = ev.Mesh(p, Lx, Ly)
    # a field linear in the minimal-image offset from (3, 2.5): exact wherever no triangle straddles the wrap
    c = np.array([3.0, 2.5])
    d = gr.min_image(p - c, Lx, Ly)
    F = np.vstack([1.0 + 0.3 * d[:, 0] - 0.2 * d[:, 1], -2 + d[:, 0], 0.5 * d[:, 1]])
    q = c + np.random.default_rng(3).uniform(-1, 1, size=(50, 2))
    val, rows, w = m.interpolate(F, q)
    dq = q - c
    exact = np.vstack([1.0 + 0.3 * dq[:, 0] - 0.2 * dq[:, 1], -2 + dq[:, 0], 0.5 * dq[:, 1]])
    assert np.abs(val - exact).max() < 1e-12
    assert np.allclose(w.sum(1), 1) and (w >= -1e-12).all()
    births = dict(q=q, rows=rows, w=w, val=val)
    worst, _ = ev.audit_births(m, F, births, n=50)
    assert worst < 1e-12, worst
    ex = ev.exact_births(m, F, births)
    assert ex["ok"], ex


def test_exact_recheck_catches_wrong_values():
    p = pts(300, 4)
    m = ev.Mesh(p, Lx, Ly)
    F = np.random.default_rng(5).normal(size=(3, 300))
    q = pts(5, 6)
    val, rows, w = m.interpolate(F, q)
    bad = val.copy()
    bad[1, 2] += 1e-6
    assert ev.exact_births(m, F, dict(q=q, rows=rows, w=w, val=val))["ok"]
    assert not ev.exact_births(m, F, dict(q=q, rows=rows, w=w, val=bad))["ok"]
    worst, _ = ev.audit_births(m, F, dict(q=q, rows=rows, w=w, val=bad), n=5)
    assert worst > ev.AUDIT_TOL


def test_homogeneous_state_unchanged_and_jump_accounting():
    N0 = 300
    p = pts(N0, 7)
    ids = np.arange(N0)
    m = ev.Mesh(p, Lx, Ly)
    E = ev.gen_events(ids, N0, 0.5, 2.0, Lx, Ly, 11)
    # homogeneous phi = 0: every node stays exactly 0
    F0 = np.zeros((3, N0))
    pos, ids2, F, acct, births = ev.apply_tick(m, F0, ids, E)
    assert np.all(F == 0) and acct["jump_ok"]
    # random field: identity I_after = I_before - R + B, ephemeral handling, row order
    F1 = np.random.default_rng(8).normal(size=(3, N0))
    pos, ids2, F, acct, births = ev.apply_tick(m, F1, ids, E)
    assert acct["jump_ok"], acct
    assert acct["n_eff_births"] + acct["n_ephemeral"] == acct["n_prop_births"]
    assert acct["n_eff_deaths"] + acct["n_ephemeral"] == acct["n_prop_deaths"]
    assert len(ids2) == N0 - acct["n_eff_deaths"] + acct["n_eff_births"]
    surv = ids2[:N0 - acct["n_eff_deaths"]]
    assert np.all(np.diff(surv) > 0)                       # survivors keep their relative order
    born = ids2[N0 - acct["n_eff_deaths"]:]
    assert np.all(np.diff(born) > 0)                       # births appended in event (= id) order
    eph = set(E["id"][E["type"] == 1]) & set(E["id"][E["type"] == 0])
    assert acct["n_ephemeral"] == len(eph) and not (eph & set(ids2.tolist()))


def test_event_stationarity_and_determinism():
    N0 = 2000
    E = ev.gen_events(np.arange(N0), N0, 0.2, 50.0, Lx, Ly, 123)
    n = N0 + np.cumsum(np.where(E["type"] == 1, 1, -1))
    assert abs(n.mean() - N0) < 4 * np.sqrt(N0)
    rate = len(E) / 50.0
    assert abs(rate - 2 * 0.2 * N0) < 0.05 * 2 * 0.2 * N0
    E2 = ev.gen_events(np.arange(N0), N0, 0.2, 50.0, Lx, Ly, 123)
    assert ev.events_hash(E) == ev.events_hash(E2)
    assert ev.regen_hash_fresh_process((np.arange(N0), N0, 0.2, 50.0, Lx, Ly, 123)) == ev.events_hash(E)
    E3 = ev.gen_events(np.arange(N0), N0, 0.2, 50.0, Lx, Ly, 124)
    assert ev.events_hash(E3) != ev.events_hash(E)
    # deaths only ever delete alive ids
    alive = set(range(N0))
    for t, typ, i in zip(E["t"], E["type"], E["id"]):
        if typ == 0:
            assert i in alive
            alive.remove(i)
        else:
            assert i not in alive
            alive.add(i)


def test_tick_slices():
    t = np.array([0.05, 0.2, 0.2000001, 0.35, 1.0])
    s = ev.tick_slices(t, 0.2, 1.0)
    assert s == [(0, 2), (2, 4), (4, 4), (4, 4), (4, 5)]
    for bad in (0.125,):
        try:
            ev.n_steps(bad, 0.05)
            assert False
        except ValueError:
            pass
    assert ev.n_steps(0.2, 0.05) == 4 and ev.n_steps(0.1, 0.025) == 4 and ev.n_steps(1.0, 0.05) == 20


def test_ros2_order_and_fixed_point():
    # order on the graph ODE (Jacobian frozen for 1 time unit, as in a frozen run), after the initial transient of the
    # continuum initial condition (in that transient the stiff order reduction is visible: ratios 2.5-3.3)
    p = pts(900, 9)
    g = gr.delaunay(p, Lx, Ly)
    h = 0.2
    P = {**P17, "b": 28.0}
    y0 = go.soliton_ic(g, (3.0, 2.5))
    s = ev.Ros2(g, h, P, 0.0125)
    for k in range(10):
        s.factor(y0)
        y0 = s.advance(y0, 80)

    def run(dt):
        s = ev.Ros2(g, h, P, dt)
        y = y0.copy()
        for k in range(2):
            s.factor(y)
            y = s.advance(y, ev.n_steps(1.0, dt))
        return y
    ref = run(0.003125)
    errs = [np.abs(run(dt) - ref).max() for dt in (0.1, 0.05, 0.025)]
    ratios = [errs[0] / errs[1], errs[1] / errs[2]]
    assert all(3.4 < r < 4.6 for r in ratios), (errs, ratios)
    # a stationary state is a fixed point of ROS2 (k1 = k2 = 0 when f(y) = 0): the homogeneous state
    s = ev.Ros2(g, h, P, 0.05)
    z = np.zeros(3 * g["N"])
    s.factor(z)
    assert np.abs(s.advance(z, 20)).max() == 0.0


def test_sparse_detection_matches_full():
    p = pts(1500, 10)
    rng = np.random.default_rng(11)
    for trial in range(4):
        centres = rng.uniform((0, 0), (Lx, Ly), size=(trial % 3 + 1, 2))
        s = np.zeros(len(p))
        for c in centres:
            d = gr.min_image(p - c, Lx, Ly)
            s += -8 * np.exp(-(d ** 2).sum(1) / 0.5)
        s[rng.integers(len(p))] = s.min()                  # an exact tie somewhere
        for thr in (-4.457582, -5.117354, -1.0):
            full = ev.detect_full(p, s, thr, Lx, Ly)
            rec = ev.sparse_record(p, np.arange(len(p)) + 1000, s, Lx, Ly)
            sp_ = ev.detect_sparse(rec, thr)
            assert len(full) == len(sp_)
            for a, b in zip(full, sp_):
                assert np.allclose(a["centre"], b["centre"], atol=1e-12) and a["depth"] == b["depth"]
                assert a["node"] == b["row"]


def test_classify_run_and_unwrapping():
    c0 = np.array([5.9, 2.0])
    mk = lambda *cs: [dict(centre=np.mod(np.array(c), (Lx, Ly))) for c in cs]
    # a particle drifting across the x boundary: 0.3 per check, 10 checks -> unwrapped displacement 3.0 (mobile)
    lists = [mk((5.9 + 0.3 * k, 2.0)) for k in range(1, 11)]
    r = ev.classify_run(lists, c0, Lx, Ly)
    assert r["cls"] == "single" and r["mobility"] == "mobile" and abs(r["max_d"] - 3.0) < 1e-9
    lists = [mk((5.9 + 0.01 * k, 2.0)) for k in range(1, 11)]
    assert ev.classify_run(lists, c0, Lx, Ly)["mobility"] == "pinned"
    # split vs remote; precedence multiplied > dissolved
    lists = [mk((5.9, 2.0)), mk((5.9, 2.0), (0.5, 2.0)), mk((5.9, 2.0), (0.6, 2.0)), []]
    r = ev.classify_run(lists, c0, Lx, Ly)
    assert r["cls"] == "multiplied" and r["sub"] == "split"
    big = [[dict(centre=np.array(c))] if len(c) == 2 and np.ndim(c[0]) == 0 else [dict(centre=np.array(x)) for x in c]
           for c in [(5.9, 2.0), [(5.9, 2.0), (12.0, 9.0)], [(5.9, 2.0), (12.0, 9.0)]]]
    assert ev.classify_run(big, c0, 20.0, 20.0)["sub"] == "remote"        # second particle 9.2 away (> 4.0)
    assert ev.classify_run([mk((5.9, 2.0)), [], []], c0, Lx, Ly)["cls"] == "dissolved"
    assert ev.classify_run([mk((5.9, 2.0)), [], mk((5.9, 2.0))], c0, Lx, Ly)["cls"] == "unresolved"
    assert ev.classify_run([mk((5.9, 2.0)), mk((5.9, 2.0), (1, 1)), mk((5.9, 2.0))], c0, Lx, Ly)["cls"] == "unresolved"
    assert ev.classify_run([mk((5.9, 2.0)), mk((2.5, 2.0))], c0, Lx, Ly)["cls"] == "unresolved"   # step > 2


def test_end_to_end_small_graph_and_scoring():
    """run_ticks on a small unregistered graph with fast rewriting, committed-record round trip, evaluator re-detection
    from sparse records, the independent recompute's detector, and score() on constructed tables."""
    import tempfile
    import run_mg5 as R
    import evaluate_mg5 as E
    import recompute_mg5 as RC
    L = 10.0
    p = np.random.default_rng(21).uniform((0, 0), (L, L), size=(3450, 2))
    F0 = go.soliton_ic(ev.Mesh(p, L, L).g, (5.0, 5.0)).reshape(3, -1)
    P = {**P17, "b": 28.0}
    events = ev.gen_events(np.arange(len(p)), len(p), 0.2, 20.0, L, L, 5)
    res, checks, ticks = R.run_ticks(p, F0, 0.2, P, R.THETA_28, 0.2, events, 0.2, 0.05, 20.0, L, L, seed_audit=3)
    assert res["status"] == R.COMPLETE, res
    assert res["max_jump_resid"] < ev.JUMP_TOL and res["max_audit"] < ev.AUDIT_TOL, res
    assert len(ticks["t"]) == 99 and len(checks) == 3
    n_eff = sum(ticks["n_eff_births"]) - sum(ticks["n_eff_deaths"])
    assert res["N_final"] == len(p) + n_eff
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "x.npz")
        R.save_checks(path, checks, ticks)
        recs, raw = R.load_checks(path)
        assert len(recs) == 3 and raw["tick_t"].shape == (99,)
        mine = E.classes_from_checks(path, R.THETA_28)
        summ = R.summarize(checks, R.THETA_28, L, L, 20.0)
        assert mine["counts"] == summ["counts"] and mine["cls"]["cls"] == summ["cls"]["cls"]
        for r, c in zip(recs, checks):
            a = RC.detect(r["pos"], r["s"], r["rows"], R.THETA_28, L, L)
            b = c["parts"]
            assert len(a) == len(b)
            for x, y in zip(a, b):
                assert np.allclose(x["centre"], y["centre"], atol=1e-12)
    # scoring rules on constructed tables
    S = lambda d=0.05, mob="pinned": dict(cls="single", max_d=d, mobility=mob)
    M, Dd = dict(cls="multiplied"), dict(cls="dissolved")
    def table(A, B):
        t = {}
        for nu, cl in A.items():
            for i, c in enumerate(cl):
                t[f"A{nu}{i}"] = dict(part="A", nu=nu, status=R.COMPLETE, cls=c)
        for nu, cl in B.items():
            for i, c in enumerate(cl):
                t[f"B{nu}{i}"] = dict(part="B", nu=nu, status=R.COMPLETE, cls=c)
        return t
    ok_ctrl = [S()] * 4
    r = E.score(table({0.0: ok_ctrl, 0.002: [S()] * 4, 0.02: [S(2.0, "mobile")] * 3 + [M], 0.2: [Dd] * 4},
                      {0.0: [M] * 4, 0.2: [S()] * 3 + [M]}))
    assert r["gate"] == "PARTIAL" and r["largest_surviving_rate"] == 0.02
    assert r["partA"]["rates"]["0.02"]["mobility"] == "MOBILE" and r["partB"]["outcome"] == "A1"
    r = E.score(table({0.0: [S(0.2)] + [S()] * 3, 0.002: [S()] * 4, 0.02: [S()] * 4, 0.2: [S()] * 4},
                      {0.0: [S()] * 4, 0.2: [Dd] * 4}))
    assert r["gate"] == "INVALID-CONTROL" and r["partB"]["outcome"] == "A-void"
    r = E.score(table({0.0: ok_ctrl, 0.002: [S()] * 4, 0.02: [S()] * 4, 0.2: [S()] * 4},
                      {0.0: [M] * 4, 0.2: [Dd] * 4}))
    assert r["gate"] == "PASS" and r["partB"]["outcome"] == "A3"
    t = table({0.0: ok_ctrl, 0.002: [S()] * 4, 0.02: [S()] * 4, 0.2: [S()] * 4}, {0.0: [M] * 4, 0.2: [M] * 4})
    t["A0.20"]["status"] = R.FAILED
    r = E.score(t)
    assert r["gate"] == "INCOMPLETE" and r["partB"]["outcome"] == "A2"
    r = E.score(table({0.0: ok_ctrl, 0.002: [M] * 4, 0.02: [Dd] * 4, 0.2: [M, M, S(), S()]}, {}))
    assert r["gate"] == "FAIL" and r["partB"]["outcome"] == "NOT SCORED"


def test_registered_start_states():
    """The start mesh equals MG4c's graph; the Part B shift is -0.9/14 on phi_Y only (no dynamics are run)."""
    import run_mg5 as R
    g, h, F, shift = R.start_state(0, 0, 28.0)
    assert np.all(shift == 0) and abs(h - 0.20007382779836555) < 1e-15 and F.shape == (3, 18901)
    m = ev.Mesh(g["pos"], g["Lx"], g["Ly"])
    assert np.array_equal(m.g["I"], g["I"]) and np.array_equal(m.g["J"], g["J"])
    g, h, F29, shift = R.start_state(1, 1, 29.0)
    assert shift[0] == 0 and shift[1] == 0 and abs(shift[2] + 0.9 / 14) < 1e-15
    assert len(R.run_specs()) == 24 and len({s["key"] for s in R.run_specs()}) == 24


if __name__ == "__main__":
    n = 0
    for k, f in list(globals().items()):
        if k.startswith("test_") and callable(f):
            t = time.time()
            f()
            n += 1
            print(f"ok  {k}  ({time.time() - t:.1f} s)", flush=True)
    print(f"{n} tests passed")
