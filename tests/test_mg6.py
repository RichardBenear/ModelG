"""MG6 tests (proposal §7, §8). Small unregistered graphs and synthetic fields only; no registered graph or stream
is run.
    python tests/test_mg6.py
"""
import os, sys, time, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import evolve as ev
import graphs as gr
import graph_ode as go
import mg6
import quad6_check as qc
from lattice import P17

L = 10.0


def small(seed=21, n=3450):
    p = np.random.default_rng(seed).uniform((0, 0), (L, L), size=(n, 2))
    F0 = go.soliton_ic(ev.Mesh(p, L, L).g, (5.0, 5.0)).reshape(3, -1)
    return p, F0


def test_candidates_deterministic_hash_and_rate():
    a = mg6.gen_candidates(1000, 0.2, 50.0, L, L, 7)
    b = mg6.gen_candidates(1000, 0.2, 50.0, L, L, 7)
    assert mg6.candidates_hash(a) == mg6.candidates_hash(b)
    assert mg6.regen_hash_fresh_process((1000, 0.2, 50.0, L, L, 7)) == mg6.candidates_hash(a)
    exp = 0.2 * 1000 * mg6.FBOUND * 50.0
    assert abs(len(a["tb"]) - exp) < 5 * np.sqrt(exp)
    assert np.all(np.diff(a["tb"]) >= 0) and a["tb"].max() <= 50.0
    assert abs(a["life0"].mean() - 5.0) < 0.5 and abs(a["ell"].mean() - 5.0) < 0.2
    c = mg6.gen_candidates(1000, 0.2, 50.0, L, L, 8)
    assert mg6.candidates_hash(c) != mg6.candidates_hash(a)


def test_shift_and_phi_lim():
    for s in range(20):
        v = mg6.draw_shift(s, 19.9, 20.1)
        assert np.hypot(*gr.min_image(v, 19.9, 20.1)) >= 6.0
    assert mg6.phi_lim(0) == np.inf and abs(mg6.phi_lim(4.0) - 0.1) < 1e-15 and abs(mg6.phi_lim(2.0) - 0.2) < 1e-15
    assert abs(np.exp(4.0 * mg6.phi_lim(4.0)) - mg6.FBOUND) < 1e-12


def test_quadrature_tiles_torus_discs_and_dec_zero():
    Lx, Ly = 19.922516574636518, 20.08183258569795
    g = mg6.Grid(Lx, Ly)
    assert abs(g.dA * g.nx * g.ny - Lx * Ly) <= 1e-12 * Lx * Ly
    rng = np.random.default_rng(1)
    for _ in range(5):
        c = rng.random(2) * (Lx, Ly)
        for R in (0.5,) + mg6.RADII:
            e = abs(g.disc_weights(c, R).sum() - np.pi * R * R) / (np.pi * R * R)
            assert e <= (0.005 if R == 0.5 else 0.001), (R, e)
    p, F0 = small()
    m = ev.Mesh(p, L, L)
    gs = mg6.Grid(L, L, 200, 200)
    out = mg6.measure(m, F0, (5.0, 5.0), 0.0, None, gs, len(p) / L ** 2)
    assert np.all(out["P0"] == 0) and np.all(out["Pc"] == 0)
    assert np.array_equal(out["D0"], out["N_R"] - len(p) / L ** 2 * np.pi * np.square(mg6.RADII))


def test_cross_check_agrees_with_measure():
    p, F0 = small()
    m = ev.Mesh(p, L, L)
    gs = mg6.Grid(L, L, 200, 200)
    rho0 = len(p) / L ** 2
    for kap, sh in ((4.0, None), (2.0, np.array([5.0, 4.0]))):
        a = mg6.measure(m, F0, (5.0, 5.0), kap, sh, gs, rho0)
        x = qc.cross_check(p, F0[0], L, L, (5.0, 5.0), kap, sh, rho0)
        for k, v in (("P0", a["P0"][1]), ("Pc", a["Pc"][1]), ("M_minus", a["M_minus"]), ("M_plus", a["M_plus"])):
            assert abs(v - x[k]) <= 0.01 * max(abs(v), abs(x[k])) + (0.5 if k in ("P0", "Pc") else 0.0), \
                (kap, k, v, x[k])


def test_measure_matches_brute_force_counts_and_annuli():
    p, F0 = small()
    m = ev.Mesh(p, L, L)
    gs = mg6.Grid(L, L, 200, 200)
    c = np.array([1.0, 9.5])
    out = mg6.measure(m, F0, c, 2.0, None, gs, len(p) / L ** 2)
    r = np.array([np.hypot(*gr.min_image(q - c, L, L)) for q in p])
    assert out["N_R"].tolist() == [int((r < R).sum()) for R in mg6.RADII]
    ann = [int(((r >= a) & (r < b)).sum()) for a, b in zip(mg6.ANN[:-1], mg6.ANN[1:])]
    assert out["ann_count"].tolist() == ann


def test_coupled_kappa0_rate_ephemeral_and_bitmap_rebuild():
    """kappa = 0 accepts each evaluated candidate with probability 1/F; ephemeral candidates are never evaluated or
    accepted; the node set at every boundary equals the rebuild from lifetimes, candidates and the accepted bitmap."""
    import run_mg6 as R
    p, F0 = small()
    cand = mg6.gen_candidates(len(p), 0.2, 10.0, L, L, 3)
    co = mg6.Coupled(cand, 0.0, None)
    gs = mg6.Grid(L, L, 100, 100)
    seen = {}

    def snap(t, pos, ids, F):
        seen[round(t, 6)] = ids.copy()
    ff = dict(G=lambda q: np.zeros(len(q)), centre=np.array([5.0, 5.0]))
    res, checks, ticks = R.run_ticks6(p, None, 0.0, None, R.THETA, co, 0.2, 0.05, 10.0, L, L, 0.0, None, gs,
                                      len(p) / L ** 2, measure_every=1e9, cross=False, fixed_field=ff, snap=snap)
    assert res["status"] == R.COMPLETE, res
    k = np.ceil(cand["tb"] / 0.2 - 1e-12)
    eph = cand["tb"] + cand["ell"] <= k * 0.2
    assert not np.any(co.evaluated & eph) and not np.any(co.accepted & eph)
    assert np.all(co.evaluated[(cand["tb"] <= 9.8) & ~eph])
    frac = co.accepted[co.evaluated].mean()
    n = co.evaluated.sum()
    assert abs(frac - 1 / mg6.FBOUND) < 5 * np.sqrt(frac * (1 - frac) / n)
    # snapshots at t = 10 (the final boundary) show M_(K-1): edits up to t = 9.8
    tp = 9.8
    N0 = len(p)
    alive0 = np.flatnonzero(cand["life0"] > tp)
    born = np.flatnonzero(co.accepted & (cand["tb"] <= tp) & (cand["tb"] + cand["ell"] > tp)) + N0
    assert sorted(seen[10.0].tolist()) == sorted(np.concatenate([alive0, born]).tolist())


def test_acceptance_follows_field_kappa4():
    import run_mg6 as R
    p, _ = small()
    cand = mg6.gen_candidates(len(p), 0.2, 10.0, L, L, 4)
    co = mg6.Coupled(cand, 4.0, None, keep_p=True)
    ff = dict(G=lambda q: -0.25 * np.exp(-(gr.min_image(np.asarray(q) - 5.0, L, L) ** 2).sum(1) / 2.0),
              centre=np.array([5.0, 5.0]))
    gs = mg6.Grid(L, L, 100, 100)
    res, _, _ = R.run_ticks6(p, None, 0.0, None, R.THETA, co, 0.2, 0.05, 10.0, L, L, 4.0, None, gs, len(p) / L ** 2,
                             measure_every=1e9, cross=False, fixed_field=ff)
    J = co.evaluated
    pp = co.p[J]
    z = (co.accepted.sum() - pp.sum()) / np.sqrt((pp * (1 - pp)).sum())
    assert abs(z) < 4, z
    assert pp.min() < np.exp(-1.0) / mg6.FBOUND * 1.05 and pp.max() <= 1 / mg6.FBOUND + 1e-12
    co2 = mg6.Coupled(cand, 4.0, None)
    R.run_ticks6(p, None, 0.0, None, R.THETA, co2, 0.2, 0.05, 10.0, L, L, 4.0, None, gs, len(p) / L ** 2,
                 measure_every=1e9, cross=False, fixed_field=ff)
    assert np.array_equal(co.accepted, co2.accepted)


def test_bound_breach_stops_run():
    import run_mg6 as R
    p, _ = small()
    cand = mg6.gen_candidates(len(p), 0.2, 2.0, L, L, 5)
    co = mg6.Coupled(cand, 4.0, None)
    ff = dict(G=lambda q: np.full(len(q), 0.11), centre=np.array([5.0, 5.0]))
    gs = mg6.Grid(L, L, 50, 50)
    res, _, _ = R.run_ticks6(p, None, 0.0, None, R.THETA, co, 0.2, 0.05, 2.0, L, L, 4.0, None, gs, len(p) / L ** 2,
                             measure_every=1e9, cross=False, fixed_field=ff)
    assert res["status"] == R.INVALID_BOUND and res["flags"][0]["kind"] == "bound"
    co = mg6.Coupled(cand, 2.0, None)
    res, _, _ = R.run_ticks6(p, None, 0.0, None, R.THETA, co, 0.2, 0.05, 2.0, L, L, 2.0, None, gs, len(p) / L ** 2,
                             measure_every=1e9, cross=False, fixed_field=ff)
    assert res["status"] == R.COMPLETE                       # 0.11 <= 0.4 / 2


def test_end_to_end_small_graph_records():
    """Coupled ODE run on a small unregistered graph: complete, accounting, record round trip and re-detection."""
    import run_mg6 as R
    import evaluate_mg6 as E
    p, F0 = small()
    cand = mg6.gen_candidates(len(p), 0.2, 20.0, L, L, 9)
    co = mg6.Coupled(cand, 4.0, None)
    gs = mg6.Grid(L, L, 200, 200)
    res, checks, ticks = R.run_ticks6(p, F0, 0.2, {**P17, "b": 28.0}, R.THETA, co, 0.2, 0.05, 20.0, L, L, 4.0, None,
                                      gs, len(p) / L ** 2, seed_audit=2)
    assert res["status"] == R.COMPLETE, res
    assert len(checks) == 3 and all(c.get("x_ok") for c in checks if "m" in c)
    assert res["N_final"] == len(p) + sum(ticks["n_eff_births"]) - sum(ticks["n_eff_deaths"])
    assert sum(ticks["n_acc"]) == co.accepted.sum() - 0                    # all accepted births were applied
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "x.npz")
        R.save_checks6(path, checks, ticks, co.accepted)
        with np.load(path) as z:
            dd = {k: z[k] for k in z.files}
        assert np.array_equal(np.unpackbits(dd["accepted_bits"])[:int(dd["n_candidates"])].astype(bool), co.accepted)
        lists = E.redetect(dd, R.THETA)
        assert [len(x) for x in lists] == dd["n_parts"].tolist()
        assert dd["m_D0"].shape == (3, len(mg6.RADII))


def test_run_order_and_scoring():
    import run_mg6 as R
    import evaluate_mg6 as E
    sp = R.run_specs()
    assert [s["arm"] for s in sp] == ["DEC"] * 4 + ["SH2"] * 4 + ["SH4"] * 4 + ["CP2"] * 4 + ["CP4"] * 4
    assert len({s["key"] for s in sp}) == 20

    def table(cp2, sh2, dec, single=True, cls="single"):
        t = {}
        for a, ys in (("CP2", cp2), ("SH2", sh2), ("DEC", dec), ("CP4", cp2), ("SH4", sh2)):
            for (sl, pl), y in zip(R.DRUNS, ys):
                t[f"{a}_D{sl}_p{pl}"] = dict(status=R.COMPLETE, cls=dict(cls=cls), single_in_window=single, Y_D0=y,
                                             Y_Dc=y, Y_P0=-30.0, obs_over_pred=y / -30.0, ann_dn0=[0.0] * 10,
                                             ann_phiG=[-0.2, -0.1, 0.0, 0.03, 0.02, 0.01, 0, 0, 0, 0],
                                             fine_dn0=[0.0] * 12, fine_phiG=[0.0] * 12, fine_P0=[0.0] * 12)
        return t
    assert E.score(table([-30, -32, -28, -35], [1, -2, 3, 0], [2, 0, -1, 1]))["gate"] == "DEMONSTRATED"
    assert E.score(table([-30, 5, -28, -35], [1, -2, 3, 0], [2, 0, -1, 1]))["gate"] == "NOT DEMONSTRATED"
    assert E.score(table([-30] * 4, [0] * 4, [0] * 4, single=False))["gate"] == "NOT SCORED (feedback)"
    s = E.score(table([-30, -32, -28, -35], [1, -2, 3, 0], [2, 0, -1, 1], cls="multiplied", single=False))
    assert s["feedback"]["CP2"]["outcome"] == "CHANGED"
    t = table([-30, -32, -28, -35], [1, -2, 3, 0], [2, 0, -1, 1])
    t["DEC_D0_p0"]["status"] = R.INVALID_BOUND
    s = E.score(t)
    assert s["gate"] == "INCOMPLETE" and "campaign" in s
    # falloff on a start-like profile: first zero crossing below 1.5, so no eligible annulus -> UNRESOLVED
    f = E.falloff([[0.0] * 10] * 4, [[0.0] * 10] * 4, [-0.22, -0.10, -0.0004, 0.033, 0.027, 0.009, -0.002, 0, 0, 0],
                  mg6.ANN)
    assert f["result"] == "UNRESOLVED" and f["eligible"] == [] and 1.0 < f["first_zero_crossing"] < 1.5


def test_mg021_amended_checks_from_committed_records():
    """MG021: waiver condition, checks 5-6 and final-state provenance from the committed Part 0 records."""
    import json, tempfile
    import run_mg6_part0 as P
    path = os.path.join(P.P0, "decisions.json")
    if not os.path.exists(path):
        return
    old = json.load(open(path))
    with tempfile.NamedTemporaryFile(delete=False, suffix=".npz") as f:
        f.write(b"not the final state")
    for man in (None, "0" * 64):                                   # missing row, wrong hash: both unverified
        o = P.amended_checks(old, f.name, man)
        assert o["final_state"]["verified"] is False
    assert P.amended_checks(old, f.name, P._sha256(f.name))["final_state"]["verified"] is True
    assert P.amended_checks(old, f.name + ".missing", "0" * 64)["final_state"]["verified"] is False
    os.unlink(f.name)
    w = o["waiver"]
    assert w["ok"] and [c["check"] for c in w["centre_failures"]] == [6] and w["max_depth_rel"] <= 0.03
    bad = json.loads(json.dumps(old))
    bad["tick_check"][1]["part_a"]["failures"].append(dict(check=3, why="count", na=1, nb=2))
    assert not P.waiver(bad)["ok"]
    assert set(o["bound_margin"]["max_phiG_by_run"]) == set(P.AMEND_RUNS) and abs(o["bound_margin"]["limit"] - 0.07) < 1e-12
    assert o["transient"]["run"] == "d0.1_dt0.05" and o["transient"]["P0_bar"] == old["tick_check"][1]["P0_bar"]
    assert o["transient"]["window"] in ([100.0, 200.0], [150.0, 250.0])


def test_mg021_authorise_guards():
    """MG021 review items B and C: incomplete settle or timing runs, a bound-margin breach in timing, or a projection
    above the cutoff never authorise the registered runs."""
    import run_mg6 as R
    import run_mg6_part0 as P
    good_s = dict(status=R.COMPLETE, t_reached=300.0)
    good_t = [dict(status=R.COMPLETE, t_reached=20.0, max_phiG=0.05)] * 8
    assert P.authorise({}, good_s, good_t, 60.0) == (True, None)
    assert not P.authorise({}, dict(status=R.INCOMPLETE, t_reached=150.0), good_t, 60.0)[0]
    assert not P.authorise({}, dict(status=R.COMPLETE, t_reached=150.0), good_t, 60.0)[0]
    assert not P.authorise({}, good_s, good_t[:7] + [dict(status=R.FAILED, t_reached=3.0, max_phiG=0.0)], 60.0)[0]
    assert not P.authorise({}, good_s, good_t[:7] + [dict(status=R.COMPLETE, t_reached=10.0, max_phiG=0.0)], 60.0)[0]
    assert not P.authorise({}, good_s, good_t[:7] + [dict(status=R.COMPLETE, t_reached=20.0, max_phiG=0.08)], 60.0)[0]
    ok, why = P.authorise({}, good_s, good_t, 111.11)
    assert not ok and "PI decision" in why
    assert not P.authorise({}, good_s, good_t, float("nan"))[0]


def test_mg021_budget_decision_authorises_driver():
    """MG021d: the PI's recorded budget decision lifts only the budget stop, and only up to its cutoff."""
    import run_mg6 as R
    d = os.path.join(R.OUT, "part0", "decisions_mg021.json")
    if not os.path.exists(d) or not os.path.exists(os.path.join(R.OUT, "part0", "budget_decision_mg021.json")):
        return
    p = R.timing_params()
    assert p["proceed"] and p["cutoff_core_h"] == 110.0 and p["delta"] == 0.1 and p["t_end"] == 200.0

if __name__ == "__main__":
    n = 0
    for k, f in list(globals().items()):
        if k.startswith("test_") and callable(f):
            t = time.time()
            f()
            n += 1
            print(f"ok  {k}  ({time.time() - t:.1f} s)", flush=True)
    print(f"{n} tests passed")
