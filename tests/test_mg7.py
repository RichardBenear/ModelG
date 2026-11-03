"""MG7 tests (proposal §3, §5, §7, §8). Static checks on the registered start files (no dynamics), Part 0 decision
rules on synthetic records, and the evaluator's scoring on MG6's committed records relabelled (no MG7 run is made).
    python tests/test_mg7.py
"""
import os, sys, time, json, tempfile, subprocess, shutil
import numpy as np
HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import run_mg5 as R5
import run_mg7 as R
import run_mg7_part0 as P
import evaluate_mg7 as E

REPO = os.path.join(HERE, "..")
PINNED = ("src/mg6.py", "src/run_mg6.py", "src/run_mg5.py", "src/evolve.py", "src/quad6_check.py", "src/graph_ode.py",
          "src/lattice.py", "src/graphs.py", "src/runner.py", "src/run_mg4c.py", "src/run_mg6_part0.py")


def _git(*a):
    return subprocess.check_output(["git", *a], cwd=REPO, stderr=subprocess.DEVNULL).decode().strip()


def test_reused_code_is_mg6_production_revision():
    try:
        _git("cat-file", "-e", R.MG6_REVISION)
    except Exception:
        print("    (git or revision c795770 not available; skipped)")
        return
    for f in PINNED:
        assert _git("hash-object", f) == _git("rev-parse", f"{R.MG6_REVISION}:{f}"), f"{f} differs from c795770"


def test_endpoint_files_match_manifest_and_load_without_shift():
    ok, rows = R.verify_endpoint_files()
    assert ok and len(rows) == 8
    for b in R.B_END:
        g, h, F, sha = R.start_state7(1, 0, b)
        with np.load(os.path.join(R.ROOT, R.field_relpath(b, 1, 0))) as z:
            assert np.array_equal(F, z["final"]) and F.dtype == np.float64
        g5, h5, F5, shift = R5.start_state(1, 0, b)          # MG6's loader: the b = 28 field, shifted
        assert np.array_equal(g["pos"], g5["pos"]) and h == h5
        assert not np.allclose(F, F5), "MG6's loader would have started from a different field"
    try:
        R.start_state7(0, 0, 28.0)
        raise AssertionError("b = 28 accepted")
    except ValueError:
        pass
    man = R.read_manifest()
    rp = R.field_relpath(27.5, 0, 0)
    man[rp] = dict(man[rp], sha256="0" * 64)
    try:
        R.start_state7(0, 0, 27.5, manifest=man)
        raise AssertionError("hash mismatch accepted")
    except RuntimeError:
        pass
    assert not R.verify_endpoint_files(man)[0]


def test_run_specs_keys_streams_and_params():
    sp = R.run_specs()
    keys = [s["key"] for s in sp]
    assert len(sp) == 16 and len(set(keys)) == 16
    assert keys[:4] == ["SH2_b27.50_D0_p0", "SH2_b27.50_D0_p1", "SH2_b27.50_D1_p0", "SH2_b27.50_D1_p1"]
    assert all(s["arm"] == "SH2" for s in sp[:8]) and all(s["arm"] == "CP2" for s in sp[8:])
    assert all(f"_b{s['b']:.2f}_" in s["key"] for s in sp)
    seeds = {(s["slot"], s["pl"]): set() for s in sp}
    for s in sp:
        seeds[(s["slot"], s["pl"])].add((R.cand_seed(s["slot"], s["pl"]), R.shift_seed(s["slot"], s["pl"])))
    assert all(len(v) == 1 for v in seeds.values()) and len({next(iter(v)) for v in seeds.values()}) == 4
    for b in R.B_END:
        assert R.params_for(b)["b"] == b
    try:
        R.params_for(28.0)
        raise AssertionError
    except AssertionError as e:
        assert "endpoint" in str(e)
    assert R.MARGIN == 0.14 and R.CUTOFF == 90.0 and R.RUN_CAP == 22 and (R.DELTA, R.DT) == (0.1, 0.05)


def _run(counts, depths, D0, centres=None, status="COMPLETE"):
    t = [10.0 * i for i in range(len(counts))]
    centres = centres or [(5.0, 5.0)] * len(counts)
    lists = [[dict(centre=list(c), depth=d)] * n for n, d, c in zip(counts, depths, centres)]
    return dict(res=dict(status=status), lists=lists, t=t, D0_15=D0)


def test_part0_compare_rules():
    n = 11
    base = _run([1] * n, [-8.0] * n, [-60.0] * n)
    same = _run([1] * n, [-8.0] * n, [-60.0] * n, centres=[(5.5, 5.0)] * n)        # centre 0.5 apart: not scored
    assert P.compare(base, same, 20, 20, P.DEPTH_TICK, P.D0_TICK)["ok"]
    assert P.compare(base, same, 20, 20, P.DEPTH_TICK, P.D0_TICK)["max_centre_diff_recorded"] > 0.49
    cnt = _run([1] * 5 + [2] + [1] * 5, [-8.0] * n, [-60.0] * n)
    assert not P.compare(base, cnt, 20, 20, P.DEPTH_TICK, P.D0_TICK)["ok"]
    dep = _run([1] * n, [-8.0] * 5 + [-8.3] + [-8.0] * 5, [-60.0] * n)               # 3.75%
    assert P.compare(base, dep, 20, 20, 0.04, P.D0_TICK)["ok"]
    assert not P.compare(base, dep, 20, 20, P.DEPTH_TICK, P.D0_TICK)["ok"]
    d5 = _run([1] * n, [-8.0] * n, [-60.0] * 5 + [-65.0] * 6)                         # window t = 50-100 shifts by 5
    assert P.compare(base, d5, 20, 20, P.DEPTH_TICK, P.D0_TICK)["ok"]                 # tick: <= 6 nodes
    assert not P.compare(base, d5, 20, 20, P.DEPTH_STEP, P.D0_STEP)["ok"]             # step: <= 2 nodes
    inc = _run([1] * n, [-8.0] * n, [-60.0] * n, status="INCOMPLETE")
    assert not P.compare(base, inc, 20, 20, P.DEPTH_TICK, P.D0_TICK)["ok"]


def _good_dec():
    runs = {f"r{i}": dict(status="COMPLETE") for i in range(6)}
    return dict(inputs=dict(ok=True),
                settle={R.bkey(b): dict(status="COMPLETE", n_parts_end=1) for b in R.B_END}, runs=runs,
                step_check={R.bkey(b): dict(ok=True) for b in R.B_END},
                tick_check={R.bkey(b): dict(ok=True) for b in R.B_END},
                bound_margin=dict(ok=True), timing=dict(projected_core_h_16_runs=75.0))


def test_part0_authorise_guards():
    assert P.authorise(_good_dec()) == (True, None)
    breaks = [lambda d: d["inputs"].update(ok=False),
              lambda d: d["settle"]["b27.50"].update(status="INCOMPLETE"),
              lambda d: d["settle"]["b28.50"].update(n_parts_end=2),
              lambda d: d["runs"]["r3"].update(status="INCOMPLETE"),
              lambda d: d["step_check"]["b28.50"].update(ok=False),
              lambda d: d["tick_check"]["b27.50"].update(ok=False),
              lambda d: d["bound_margin"].update(ok=False),
              lambda d: d["timing"].update(projected_core_h_16_runs=95.0),
              lambda d: d["timing"].update(projected_core_h_16_runs=float("nan"))]
    for f in breaks:
        d = json.loads(json.dumps(_good_dec()))
        f(d)
        ok, why = P.authorise(d)
        assert not ok and why
    d = _good_dec()
    d["timing"]["projected_core_h_16_runs"] = 95.0
    assert "PI decision required" in P.authorise(d)[1]


def test_part0_never_overwrites_and_driver_needs_proceed():
    tmp = tempfile.mkdtemp()
    try:
        f = os.path.join(tmp, "x.npz")
        assert P._fresh(f) == f
        open(f, "w").write("x")
        try:
            P._fresh(f)
            raise AssertionError("overwrite allowed")
        except P.Exists:
            pass
        dp = os.path.join(tmp, "decisions.json")
        json.dump(dict(proceed=False, stop="projected cost 95.0 core-h exceeds ...",
                       timing=dict(projected_core_h_16_runs=95.0)), open(dp, "w"))
        try:
            R.part0_decision(dp)
            raise AssertionError("launched without proceed")
        except SystemExit:
            pass
        json.dump(dict(pi_approved=True, cutoff_core_h=100.0), open(os.path.join(tmp, "budget_decision.json"), "w"))
        try:
            R.part0_decision(dp)                                # proceed lifted, but Delta/dt missing
            raise AssertionError
        except SystemExit as e:
            assert "Delta" in str(e)
        json.dump(dict(proceed=False, stop="projected cost 95.0", delta=0.1, dt=0.05,
                       timing=dict(projected_core_h_16_runs=95.0)), open(dp, "w"))
        assert R.part0_decision(dp)["cutoff_core_h"] == 100.0
        json.dump(dict(proceed=False, stop="check 2 (step) failed", delta=0.1, dt=0.05,
                       timing=dict(projected_core_h_16_runs=50.0)), open(dp, "w"))
        try:
            R.part0_decision(dp)                                # a budget decision never lifts a check failure
            raise AssertionError
        except SystemExit:
            pass
    finally:
        shutil.rmtree(tmp)


def _mg6_rows():
    rows = {}
    for arm in ("CP2", "SH2", "CP4", "SH4"):
        for st in E.STARTS:
            rec = json.load(open(os.path.join(E.MG6_OUT, "runs", f"{arm}_{st}.json")))
            r, bad = E.run_row(rec, E.load(os.path.join(E.MG6_OUT, "checks", f"{arm}_{st}.npz")))
            assert not bad
            rows[(arm, st)] = r
    return rows


def test_evaluator_scoring_on_relabelled_mg6_records():
    """Scoring check only: MG6's kappa = 2 pair plays b = 27.5 and its kappa = 4 pair plays b = 28.5."""
    rows = _mg6_rows()
    lab = {("CP2", 27.5): "CP2", ("SH2", 27.5): "SH2", ("CP2", 28.5): "CP4", ("SH2", 28.5): "SH4"}
    table = {f"{a}_b{b:.2f}_{st}": dict(rows[(src, st)]) for (a, b), src in lab.items() for st in E.STARTS}
    mid = {f"{a}_b28.00_{st}": rows[(a, st)] for a in ("CP2", "SH2") for st in E.STARTS}
    out = E.score(table, mid)
    c2 = np.array([rows[("CP2", s)]["Y_D0"] - rows[("SH2", s)]["Y_D0"] for s in E.STARTS])
    c4 = np.array([rows[("CP4", s)]["Y_D0"] - rows[("SH4", s)]["Y_D0"] for s in E.STARTS])
    D = c4 - c2
    t = D.mean() / (D.std(ddof=1) / 2)
    assert np.allclose(out["primary"]["values"], D) and abs(out["primary"]["t"] - t) < 1e-9
    exp = "LARGER AT HIGHER b" if (np.all(D < 0) and t <= -2.353) else "NOT DEMONSTRATED"
    assert out["gate"] == exp
    mg6ev = json.load(open(os.path.join(E.MG6_OUT, "evaluation.json")))
    assert np.allclose(out["contrasts"]["b28.00"]["c"], mg6ev["primary"]["CP2_minus_SH2"]["c"])
    assert np.allclose(out["contrasts"]["b27.50"]["c"], c2)
    assert set(out["residuals"]) == {"start_field", "in_run"} and out["falloff"]["b27.50"]["result"] == "UNRESOLVED"
    assert out["feedback"]["b28.50"]["outcome"] == "SURVIVES"
    # INCOMPLETE precedes the feedback rule
    t2 = {k: dict(v) for k, v in table.items()}
    t2["SH2_b28.50_D1_p1"]["status"] = "INCOMPLETE"
    t2["CP2_b27.50_D0_p0"]["single_in_window"] = False
    assert E.score(t2, mid)["gate"] == "INCOMPLETE"
    t2["SH2_b28.50_D1_p1"]["status"] = "COMPLETE"
    assert E.score(t2, mid)["gate"] == "NOT SCORED (feedback)"


def test_fieller():
    a = np.array([-10.0, -12.0, -9.0, -11.0])
    f = E.fieller(a, a / 2)                                   # exact ratio 2: a degenerate interval at 2
    assert f["bounded"] and abs(f["ratio"] - 2) < 1e-12 and np.allclose(f["ci95"], [2, 2])
    f = E.fieller(a, np.array([-1.0, 1.0, -0.5, 0.6]))        # denominator indistinguishable from zero
    assert not f["bounded"] and f["ci95"] is None
    rng = np.random.default_rng(3)
    b = np.array([-8.0, -9.0, -7.0, -8.5])
    f = E.fieller(b * 1.5 + rng.normal(0, 0.3, 4), b)
    assert f["bounded"] and f["ci95"][0] < f["ratio"] < f["ci95"][1]


if __name__ == "__main__":
    n = 0
    for k, f in list(globals().items()):
        if k.startswith("test_") and callable(f):
            t = time.time()
            f()
            n += 1
            print(f"ok  {k}  ({time.time() - t:.1f} s)", flush=True)
    print(f"{n} tests passed")
