"""MG8 tests (spec/MG8_proposal.md draft 4; §4.1, §5). Run: python tests/test_mg8.py
Only unregistered streams are used (seed_for(8, 9, ...) or small literal seeds); no production run is performed."""
import os, sys, json, tempfile, itertools
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
import mg8
from mg8 import J, R, L, P, TAU0, H
from lattice import tau_leap_euler_1d, euler_mean_map
from runner import seed_for
from measure import fourier_shift, analyze_run

PM = mg8.params()
REF, C = mg8.reference_start()
ARGS = (PM["a"], PM["b"], PM["g"], PM["p"], PM["q"], PM["dx"], PM["dy"])
S9 = seed_for(8, 9, 0, 0, 0)


def counts(Om):
    return [np.rint(Om * C[k]).astype(np.int64) for k in range(3)]


# ------------------------------------------------------------------ §4.1 implementation checks
def test_J_regression_stochastic():
    a, b = counts(1000.0), counts(1000.0)
    r1 = tau_leap_euler_1d(a[0], a[1], a[2], 1000.0, *ARGS, H, TAU0, 400, 0.0, 0.0, REF["x"], 100, S9)
    r2 = mg8.run_stoch(b[0], b[1], b[2], J, 1000.0, *ARGS, H, TAU0, 400, 100, S9)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    assert all(np.array_equal(r1[i], r2[i]) for i in range(3))


def test_J_segmented_equals_single_call():
    a, b = counts(1000.0), counts(1000.0)
    mg8.run_stoch(a[0], a[1], a[2], R, 1000.0, *ARGS, H, TAU0, 400, 400, S9)
    mg8.run_stoch(b[0], b[1], b[2], R, 1000.0, *ARGS, H, TAU0, 200, 200, S9)
    mg8.run_stoch(b[0], b[1], b[2], R, 1000.0, *ARGS, H, TAU0, 200, 200, -1)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))


def test_det_J_equals_euler_mean_map():
    G, X, Y = (C[k].copy() for k in range(3))
    X[10] += 0.3; Y[40] -= 0.2
    G2, X2, Y2 = G.copy(), X.copy(), Y.copy()
    clips = mg8.run_det(G, X, Y, J, *ARGS, H, TAU0, 1500, 1)
    G2, X2, Y2 = euler_mean_map(PM, H, TAU0, 1500 * TAU0, G2, X2, Y2)
    assert clips.sum() == 0
    assert max(np.abs(G - G2).max(), np.abs(X - X2).max(), np.abs(Y - Y2).max()) <= 1e-12


def test_mirror_stochastic_and_deterministic():
    s = counts(1000.0); m = [v[::-1].copy() for v in s]
    mg8.run_stoch(s[0], s[1], s[2], R, 1000.0, *ARGS, H, TAU0, 300, 300, S9)
    mg8.run_stoch(m[0], m[1], m[2], L, 1000.0, *ARGS, H, TAU0, 300, 300, S9)
    assert all(np.array_equal(x[::-1], y) for x, y in zip(s, m))
    d = [C[k].copy() for k in range(3)]; e = [C[k][::-1].copy() for k in range(3)]
    mg8.run_det(d[0], d[1], d[2], R, *ARGS, H, TAU0, 3000, S9)
    mg8.run_det(e[0], e[1], e[2], L, *ARGS, H, TAU0, 3000, S9)
    assert max(np.abs(x[::-1] - y).max() / np.abs(x).max() for x, y in zip(d, e)) <= 1e-12


def test_conservation_without_reactions():
    """No source or sink (a = b = p = q = 0, g = 1: only X->G, Y->X conversions and hops): the total token count is
    conserved by every scheduler."""
    for sch in (J, R, L, P):
        s = counts(1000.0); tot = sum(v.sum() for v in s)
        mg8.run_stoch(s[0], s[1], s[2], sch, 1000.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 12.0, H, TAU0, 200, 200, S9)
        assert sum(v.sum() for v in s) == tot


# ------------------------------------------------------------------ §5.1 operator benchmark
def test_conditional_matrices_two_constructions():
    rng = np.random.default_rng(1)
    for r in (0.01, 0.05, 0.2):
        for s in range(12):
            o = np.array([(s + k) % 12 for k in range(12)])
            assert np.abs(mg8.map_matrix(o, r) - mg8.chain_matrix(o, r)).max() <= 1e-12
        for _ in range(30):
            o = rng.permutation(6)
            assert np.abs(mg8.map_matrix(o, r) - mg8.chain_matrix(o, r)).max() <= 1e-12


def _avg_R(N, r):
    return np.mean([mg8.map_matrix(np.array([(s + k) % N for k in range(N)]), r) for s in range(N)], axis=0)


def _avg_P(N, r):
    return np.mean([mg8.map_matrix(np.array(o), r) for o in itertools.permutations(range(N))], axis=0)


def test_refinement_ratios():
    """§5.1(c): N = 12 (R) and 6 (P), t = 1, tau in {0.2, ..., 0.0125}; Frobenius norm; last two ratios in [0.4, 0.6]."""
    for avg, N in ((_avg_R, 12), (_avg_P, 6)):
        D = []
        for tau in (0.2, 0.1, 0.05, 0.025, 0.0125):
            n = int(round(1 / tau))
            TJ = mg8.map_matrix(np.arange(N), tau, sync=True)
            D.append(np.linalg.norm(np.linalg.matrix_power(avg(N, tau), n) - np.linalg.matrix_power(TJ, n)))
        rat = [D[i + 1] / D[i] for i in range(4)]
        assert all(x < 1 for x in rat)
        assert all(0.4 <= x <= 0.6 for x in rat[-2:]), rat


def test_displacement_distribution_moments():
    for r in (0.05, 0.2):
        for s in range(12):
            o = [(s + k) % 12 for k in range(12)]
            for j in (0, 5):
                A = mg8.displacement_dist(o, r, j, 12)
                d = np.arange(-12, 13)
                assert abs(A.sum() - 1) < 1e-12 and abs((d * A).sum()) < 1e-12
    # infinite-ring formulas (A4) vs a seam far from the token: variance 2r + 2r^2/(1-r)
    r = 0.2
    A = mg8.displacement_dist(list(range(40)), r, 20, 40)      # token far from the seam
    d = np.arange(-40, 41)
    assert abs((d ** 2 * A).sum() - (2 * r + 2 * r * r / (1 - r))) < 1e-9
    assert abs((d ** 3 * A).sum() - 6 * r * r / (1 - r) ** 2) < 1e-9


def test_stochastic_benchmark_conditional():
    """§5.1(b) in miniature: per-replicate chi-square against the conditional-order distribution."""
    from scipy import stats
    from part0_mg8 import pooled_chisq
    N, n0, dmax = 12, 10 ** 5, 12
    for sch in (R, P):
        ps = []
        for rep in range(60):
            cnt, order = mg8.bench_one_step(n0, N, 0.2, sch, seed_for(8, 9, 1, sch, rep), dmax)
            assert cnt.sum() == n0
            A = mg8.displacement_dist(list(order), 0.2, 0, dmax)
            ps.append(pooled_chisq(cnt, n0 * A))
        assert stats.kstest(ps, "uniform").pvalue >= 1e-3 and min(ps) > 1e-6


# ------------------------------------------------------------------ §5.4 synthetic validation
def test_shape_error_shift_and_deformation():
    phi = REF["phi"]
    for s in (0.037, -0.081, 0.6):
        E, sg, _ = mg8.shape_error(fourier_shift(phi, H, -s), phi, H)
        assert E < 1e-8 and abs(sg - s) < 1e-6
    x = REF["x"]
    for amp in (1e-3, 1e-2):
        odd = x * np.exp(-x ** 2 / 4)
        A = phi + amp * np.outer(np.abs(phi).max(axis=1), odd / np.abs(odd).max())
        A = fourier_shift(A, H, -0.05)
        E, sg, _ = mg8.shape_error(A, phi, H)
        nB = np.sqrt((phi ** 2).sum(axis=1))
        grid = 0.05 + np.arange(-0.2, 0.2, 1e-5 * H * 20)          # brute force (coarser than spec; spec grid in Part 0)
        bf = min((float(np.max(np.sqrt(((fourier_shift(A, H, g) - phi) ** 2).sum(axis=1)) / nB)), g) for g in grid)
        assert E <= bf[0] * (1 + 1e-6) and abs(E - bf[0]) <= 0.01 * bf[0] and abs(sg - bf[1]) <= 1e-3


def test_centre_velocity_with_wrap():
    x = REF["x"]
    for v in (1e-3, 1e-5):
        t = np.arange(0, 400.0)
        start = 49.9 if v > 1e-4 else 0.0
        c = [mg8.centre(np.sin(0) + fourier_shift(-np.exp(-x ** 2), H, -(start + v * ti) % 100), x, H) for ti in t[::20]]
        vv = mg8.slope(t[::20], mg8.unwrap(c, 100.0))
        assert abs(vv - v) <= 0.02 * v


def _synthetic_counts(positions_per_snap, Om=1e6):
    """Counts with the reference particle at each listed position (list per snapshot)."""
    S0 = REF["S0"]
    out = []
    for ps in positions_per_snap:
        f = np.repeat(S0[:, None], C.shape[1], axis=1).astype(float)
        for p in ps:
            f = f + fourier_shift(REF["phi"], H, -p)
        out.append(np.rint(Om * f).astype(np.int64))
    a = np.array(out)
    return a[:, 0], a[:, 1], a[:, 2]


def test_tracking_and_fate_synthetic():
    import run_mg8 as D
    from measure import smooth_periodic, detect
    lam = REF["lam"]
    t = np.arange(1, 121, dtype=float)
    cases = {
        "moving_wrap": [[((0.9 * ti + 50) % 100) - 50] for ti in t],          # wraps at t ~ 56
        "loss": [[0.0] if ti < 40 else [] for ti in t],
        "replication": [[0.0] if ti < 40 else [0.0, 1.5 * lam] for ti in t],
        "jump": [[0.0] if ti < 40 else [2.5] for ti in t],
    }
    for name, pos in cases.items():
        G, X, Y = _synthetic_counts(pos)
        a = analyze_run(G, X, Y, t, 1e6, H, REF["S0"], REF["A_s"], lam, t_discard=10.0, K=2, T=120.0, min_snap_per_block=1)
        ps = smooth_periodic(X / 1e6 - REF["S0"][1], H)
        dets = [detect(ps[s], REF["x"], H, REF["thr"], lam) for s in range(len(t))]
        tf, _ = D.track_fail(dets, np.array(a["positions"], float))
        kinds = {e[0] for e in a["events"]}
        if name == "moving_wrap":
            assert not tf and not (kinds & {"loss", "replication", "merge"})
        elif name == "loss":
            assert "loss" in kinds and not tf
        elif name == "replication":
            assert "replication" in kinds
        else:
            assert tf


def test_labels_production():
    f = mg8.label_production
    assert f(False, 0, 0, 0, 0, True)["fate"] == "CHANGED"
    assert f(True, 0.02, 1e-9, 0, 1e-8, True)["summary"] == "SENSITIVE"
    assert f(True, 5e-3, 1e-9, 1e-8, 1e-8, True)["summary"] == "MINOR"
    assert f(True, 1e-4, 1e-9, 1e-8, 1e-8, True)["summary"] == "PRACTICALLY INSENSITIVE"
    assert f(True, 1e-4, 1e-9, 2e-6, 1e-8, True)["summary"] == "SENSITIVE"
    assert f(True, 0.5, 2e-4, 0, 1e-8, True)["shape"] == "INSUFFICIENT RESOLUTION"
    assert f(True, 1e-4, 1e-9, 0, 1e-6, True)["drift"] == "INSUFFICIENT RESOLUTION"
    assert f(True, 1e-4, 1e-9, 0, 1e-8, False)["summary"] == "UNRESOLVED"
    lab = [f(True, 1e-4, 1e-9, 0, 1e-8, True)] * 2 + [f(True, 0.02, 1e-9, 0, 1e-8, True)]
    assert mg8.combine_replicates(lab) == "SEAM-DEPENDENT"


def _ladder(E, floor=1e-10, S=0.0, single=True):
    return {k: dict(Ek=E[k], Efloor=floor, S=S, all_single=single, all_settled=True) for k in range(7)}


def test_labels_refinement():
    hz = lambda E: dict(single=True, settled=True, E300=E, E600=E, Efloor=1e-10)
    lin = [1e-2 / 2 ** k for k in range(7)]
    assert mg8.label_refinement(_ladder(lin), hz(lin[6]))[0] == "CONVERGES TO J"
    flat = [1e-2] * 7
    assert mg8.label_refinement(_ladder(flat), hz(1e-2))[0] == "NOT RESOLVED BY THIS LADDER"
    slow = [1e-2 / 2 ** (0.3 * k) for k in range(7)]
    assert mg8.label_refinement(_ladder(slow), hz(slow[6]))[0] == "NOT RESOLVED BY THIS LADDER"
    tiny = [1e-6, 3e-7, 8e-8, 2e-8, 5e-9, 1e-9, 5e-10]
    assert mg8.label_refinement(_ladder(tiny, floor=1e-9), dict(single=True, settled=True, E300=5e-10, E600=5e-10,
                                                                Efloor=1e-9))[0] == "BELOW MEASUREMENT RESOLUTION"
    assert mg8.label_refinement(_ladder(tiny, floor=2e-5), hz(5e-10))[0] == "NOT RESOLVED BY THIS LADDER"
    assert mg8.label_refinement(_ladder(lin, S=0.5), hz(lin[6]))[0] == "NOT RESOLVED BY THIS LADDER"
    bad_h = dict(single=True, settled=True, E300=lin[6], E600=2 * lin[6], Efloor=1e-10)
    assert mg8.label_refinement(_ladder(lin), bad_h)[0] == "NOT RESOLVED BY THIS LADDER"


def test_b2_statistics_synthetic():
    import evaluate_mg8 as EV
    prof = np.array(REF["phi"][:, np.abs(REF["x"]) < 10])
    rng = np.random.default_rng(3)
    runs = {}
    for a in ("J", "R", "P"):
        for i in range(8):
            runs[f"sto_{a}_k0_r{i}"] = dict(events=[] if i else [("loss", 300.0)], eligible=i > 0, track_fail=False,
                                           v=float(rng.normal(0, 1e-4)), profile=(prof + rng.normal(0, 1e-4, prof.shape)).tolist(),
                                           positions=list(rng.normal(0, 0.01, 500)), farvar=[1, 2, 3], clip_frac=[0, 0, 0])
    res = EV.evaluate_b2(runs)
    assert res["arms"]["R"]["n_eligible"] == 7 and res["arms"]["R"]["runs_with_event"]["loss"] == 1
    assert res["R_vs_J"]["drift_diff"]["verdict"] == "equivalent within margin"
    assert res["R_vs_J"]["E_stoch"]["hi"] < 1e-2 and res["J_split_half"] is not None


def test_driver_smoke():
    """det_run and sto_run on a temp dir with unregistered seeds and tiny horizons (not production runs)."""
    import run_mg8 as D
    tmp = tempfile.mkdtemp()
    D.RUNS, D.FIELDS, D.SNAP, D.FULL = (os.path.join(tmp, d) for d in ("runs", "fields", "snap", "full"))
    for d in (D.RUNS, D.FIELDS, D.SNAP, D.FULL):
        os.makedirs(d)
    D.seed = lambda job: seed_for(8, 9, 7, job["sched"], job["rep"])
    r = D.det_run(dict(kind="det", sched=R, k=1, rep=0, T=50))
    assert len(r["t"]) == 50 and r["finite"] and os.path.exists(os.path.join(D.FIELDS, r["key"] + ".npz"))
    s = D.sto_run(dict(kind="sto", sched=P, k=0, rep=0, T=20))
    assert "eligible" in s and len(s["dets"]) == 20
    assert len(D.jobs()) == 55       # 57 registered runs; the two continuations are the k = 6 J and R r0 records run to 600


def test_pooled_chisq_tail_pooling():
    """Every pooled cell has expected >= 5; a tail excess of 4 tokens over an expectation of 0.3 is not flagged
    below 1e-6 (attempt-1 failure mode)."""
    from part0_mg8 import pooled_chisq
    e = np.array([0.01, 0.2, 50000.0, 902500.0, 45125.0, 2256.25, 112.8, 5.64, 0.28, 0.014])
    o = np.array([0, 0, 50122, 902415, 45120, 2221, 113, 5, 4, 0])
    assert pooled_chisq(o, e) > 1e-3
    assert pooled_chisq(np.rint(e), e) > 0.5


def test_evaluate_b1_synthetic():
    """Evaluator B1 path on fabricated records: R = J shifted plus an odd deformation proportional to tau."""
    import evaluate_mg8 as EV
    tmp = tempfile.mkdtemp()
    os.makedirs(os.path.join(tmp, "runs")); os.makedirs(os.path.join(tmp, "fields"))
    phi, x = REF["phi"], REF["x"]
    odd = np.outer(np.abs(phi).max(axis=1), x * np.exp(-x ** 2 / 4))
    def write(key, T, fields_fn, centre_v=0.0):
        t = list(np.arange(1, T + 1) * 1.00005)
        rec = dict(t=t, centre=[centre_v * ti for ti in t], n_dets=[1] * T, finite=True)
        json.dump(rec, open(os.path.join(tmp, "runs", key + ".json"), "w"))
        np.savez(os.path.join(tmp, "fields", key + ".npz"), **{f"t{u}": fields_fn(u) for u in range(50, T + 1, 50)})
    for k in range(7):
        T = 600 if k == 6 else 300
        amp = 1e-2 / 2 ** k
        write(f"det_J_k{k}_r0", T, lambda u: phi)
        for j in range(3):
            write(f"det_R_k{k}_r{j}", T if j == 0 else 300,
                  lambda u, a=amp: fourier_shift(phi + a * odd, H, -0.03))
    ev = EV.evaluate_b1(tmp, EV.load(tmp))
    assert ev["levels"][0]["summary"] == "SENSITIVE" and ev["levels"][0]["labels"][0]["shape"] == "MATERIAL"
    assert ev["refinement"] == "CONVERGES TO J", ev.get("refinement_info")


if __name__ == "__main__":
    names = [n for n in dir() if n.startswith("test_")]
    for n in names:
        globals()[n]()
        print("ok", n)
    print(f"{len(names)} tests passed")
