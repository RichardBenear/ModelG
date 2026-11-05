"""MG9d tests (spec/MG9d_proposal.md draft 3 §5). Run: python tests/test_mg9d.py
No MG9d integration is performed: the end-to-end test uses a fake ODE (exact translation), and the wrapper test
integrates MG003's b = 29 reference for 5 time units only (not an MG9d state, not ledgered)."""
import os, sys, json, shutil, tempfile
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import run_mg9 as R
import run_mg9d as D

VF = 5e-9                     # the floor validation's value (delta ~1e-9 -> minimum 5e-9)


# ------------------------------------------------------------------ solver wrapper and maps
def test_wrapper_bitwise_matches_mg9_ode():
    ref = np.load(os.path.join(D.ROOT, "results", "MG003", "ref_profiles.npz"))["h0.2_b29.0"]
    y0 = np.concatenate([ref[1], ref[2], ref[3]])
    te = np.arange(1.0, 5.0 + 1e-9, 1.0)
    a = R.mg9_ode(29.0, 500, 0.2, 5.0, y0, seed_scale=0.0, t_eval=te)
    b = D.mg9d_ode(29.0, 500, 0.2, 5.0, y0, seed_scale=0.0, t_eval=te, max_step=0.5)
    assert np.array_equal(a["y"], b["y"]) and a["residual"] == b["residual"] and a["t_reached"] == b["t_reached"]


def test_exact_e_pat_and_source_hash():
    assert D.e_pat() == 6.60086115019687
    D.source_state()                                      # raises Stop on a hash mismatch


def test_nyquist_filter():
    src = D.source_state()
    y0 = D.nyquist_filter(src, D.N02)
    for c in range(3):
        assert abs(np.fft.rfft(y0[c * 500:(c + 1) * 500])[-1]) < 1e-9
    assert np.abs(y0 - src).max() <= D.Y0_MAX_CHANGE
    assert np.abs(D.nyquist_filter(y0, D.N02) - y0).max() < 1e-13


def test_shift_reflection_and_d01_maps():
    y0 = D.common_start()
    for s in D.SHIFTS.values():
        assert np.abs(R.fourier_shift(R.fourier_shift(y0, 500, 0.2, s), 500, 0.2, -s) - y0).max() <= 1e-12
    assert np.array_equal(D.reflect(D.reflect(y0, 500), 500), y0)
    x = D.xgrid(500, 0.2)
    g = np.exp(-(x - 1.0) ** 2)
    gr = D.reflect(np.concatenate([g, g, g]), 500)[:500]
    assert np.abs(gr - np.exp(-(x + 1.0) ** 2)).max() < 1e-12          # x -> -x about node x = 0
    y1 = D.to_d01(y0)
    assert len(y1) == 3000 and np.abs(D.from_d01_nodes(y1) - y0).max() <= D.R4_BACK_TOL


def test_transfer_checks_pass():
    out = D.transfer_checks()
    assert out["ok"], out
    assert set(out["start_sha256"]) == set(D.ORDER)


# ------------------------------------------------------------------ measurement
def test_align_recovers_shift_and_flags_bound():
    y0 = D.common_start()
    f0 = y0[500:1000]
    for s in (0.003, -0.08, 0.15):
        d, S, amb = D.align(f0, D.shift_field(f0, 0.2, s), 0.2, 0.0)
        assert abs(d - s) < 1e-8 and S < 1e-8 and not amb
    d, S, amb = D.align(f0, D.shift_field(f0, 0.2, D.LAM / 4 + 0.05), 0.2, 0.0)
    assert amb                                            # optimum pinned at the bound -> ALIGNMENT-AMBIGUOUS


def test_floor_validation_small_and_vf_rule():
    out = D.floor_validation()
    assert out["delta"] < 1e-8 and out["v_f"] == 5e-9
    assert np.isclose(D.vf_from_delta(1e-5), 5e-9) and np.isclose(D.vf_from_delta(1e-4), 5e-8)
    for c in out["cases"].values():
        assert c["matched"] == c["of"] and not c["ambiguous"]


# ------------------------------------------------------------------ decision rules
def test_compare_bands():
    assert D.compare(-4.1e-7, -4e-7, VF) == "AGREE"
    assert D.compare(-3e-7, -4e-7, VF) == "DIFFER"
    assert D.compare(-3.5e-7, -4e-7, VF) == "INCONCLUSIVE"
    assert D.compare(1e-10, -1e-10, VF) == "AGREE"          # a tiny sign change near zero is not a difference


def test_c1_truth_table():
    s0 = dict(a=-4e-7, b=-3.5e-7, v=-3.75e-7)
    assert D.c1(s0, dict(s0), True, True, VF) == "NUMERICS-AGREE"
    assert D.c1(s0, dict(s0, b=-2e-7), True, True, VF) == "NUMERICS-SENSITIVE"
    assert D.c1(s0, dict(s0, b=-3.1e-7), True, True, VF) == "INCONCLUSIVE"
    assert D.c1(s0, dict(s0), True, False, VF) == "INCONCLUSIVE"


def test_c3_reflection_sign_applied_once_and_not_evaluated():
    assert D.c3(-4e-7, +4e-7, True, True, VF) == "SYMMETRIC"            # raw R3 moves right = left in R0 frame
    assert D.c3(-4e-7, -4e-7, True, True, VF) == "ASYMMETRY-FLAG"       # negating twice would pass this
    assert D.c3(-4e-7, +4e-7, True, False, VF) == "NOT EVALUATED"
    assert D.reflect_adjust(3e-7) == -3e-7


def test_c2():
    arms = [(-4e-7, True)] * 3
    assert D.c2(-4e-7, True, arms, VF) == "ALIGNMENT-INDEPENDENT"
    assert D.c2(-4e-7, True, [(-4e-7, True), (-2e-7, True), (-4e-7, False)], VF) == "ALIGNMENT-DEPENDENT"
    assert D.c2(-4e-7, True, [(-4e-7, True), (-2e-7, False), (-4e-7, True)], VF) == "INCONCLUSIVE"
    assert D.c2(-4e-7, False, arms, VF) == "INCONCLUSIVE"


def test_c4_all_branches_and_worked_cases():
    L = lambda a, b: D.c4_label(a, b, VF)
    assert L(1e-9, -1e-9) == "NO RESOLVED MOTION"
    assert L(1e-9, -4e-7) == "ONSET"
    assert L(+4e-7, -4e-7) == "REVERSAL"                   # worked case 1
    assert L(-4e-7, +1e-9) == "DECREASED BELOW RESOLUTION" # worked case 3: no ratio formed
    assert L(-4e-7, -3e-7) == "DECREASING"
    assert L(-4e-7, -4e-7) == "PERSISTENT"
    assert L(-4e-7, -4.8e-7) == "INCREASING"
    assert L(-4e-7, -3.7e-7) == "INCONCLUSIVE"             # rho 0.925, in the gap
    s0, s1 = dict(a=-4e-7, b=-3e-7, v=-3.5e-7), dict(a=-4e-7, b=-4e-7, v=-4e-7)
    assert D.c4(s0, s0, True, True, "NUMERICS-AGREE", VF) == "DECREASING"
    assert D.c4(s0, s1, True, True, "NUMERICS-AGREE", VF) == "INCONCLUSIVE"   # R1 must give the same label


def test_c5_ordered_and_exclusive():
    assert D.c5(-4e-7, -4e-7, False, VF) == "R4 INVALID / LOST-IDENTITY"
    assert D.c5(1e-10, -1e-10, True, VF) == "NOT EVALUATED"                   # the draft-1 overlap case
    assert D.c5(-4e-7, +1e-7, True, VF) == "REVERSED"
    assert D.c5(-4e-7, -3.5e-7, True, VF) == "SIMILAR"
    assert D.c5(-4e-7, -1e-7, True, VF) == "REDUCED"
    assert D.c5(-4e-7, -2.5e-7, True, VF) == "OTHER"
    assert D.c5(-4e-7, -9e-7, True, VF) == "OTHER"


def _lab(**kw):
    lab = dict(C1="NUMERICS-AGREE", C2="ALIGNMENT-INDEPENDENT", C3="SYMMETRIC", C4="DECREASING", C5="SIMILAR")
    lab.update(kw)
    return lab


def test_overall_order_and_worked_cases():
    sp = dict(a=-4e-7, b=-3e-7, v=-3.5e-7)
    O = lambda lab, s0=sp, s1=sp, v0=True, v1=True: D.overall(lab, s0, s1, v0, v1, VF)
    rev = dict(a=4e-7, b=-4e-7, v=0.0)
    assert O(_lab(C4="REVERSAL"), rev, rev) == "INCONCLUSIVE (C4 REVERSAL)"        # worked case 1
    assert O(_lab(C3="NOT EVALUATED")) == "INCONCLUSIVE"                            # worked case 2
    assert O(_lab(C3="ASYMMETRY-FLAG")) == "INCONCLUSIVE"
    assert O(_lab(), v0=False) == "INCONCLUSIVE"
    assert O(_lab(C1="NUMERICS-SENSITIVE")) == "NUMERICAL ACCURACY UNRESOLVED"
    z = dict(a=1e-9, b=-1e-9, v=0.0)
    assert O(_lab(C4="NO RESOLVED MOTION", C5="NOT EVALUATED"), z, z) == "NO RESOLVED DRIFT"
    assert O(_lab(C4="NO RESOLVED MOTION", C5="NOT EVALUATED"), z, dict(z, b=-4e-7)) != "NO RESOLVED DRIFT"
    assert O(_lab(C2="ALIGNMENT-DEPENDENT")) == "MESH EFFECT SUPPORTED"
    assert O(_lab(C5="REDUCED")) == "MESH EFFECT SUPPORTED"
    assert O(_lab(C4="REVERSAL", C2="ALIGNMENT-DEPENDENT"), rev, rev) == "MESH EFFECT SUPPORTED"  # not hidden
    assert O(_lab()) == "RELAXATION SUPPORTED [C4 DECREASING]"
    assert O(_lab(C4="DECREASED BELOW RESOLUTION")) == "RELAXATION SUPPORTED [C4 DECREASED BELOW RESOLUTION]"
    assert O(_lab(C5="R4 INVALID / LOST-IDENTITY")).endswith("(on D02; the finer-mesh comparison is unavailable)")
    assert O(_lab(C4="PERSISTENT")) == "MOVING STATE SUPPORTED"
    assert O(_lab(C5="OTHER")) == "INCONCLUSIVE (C4 DECREASING)"
    assert O(_lab(C2="INCONCLUSIVE")) == "INCONCLUSIVE (C4 DECREASING)"


# ------------------------------------------------------------------ end to end with a fake ODE (no integration)
def _fake_ode(v):
    def ode(b, N, h, t_end, y, *, seed_scale, t_eval, t_start, **kw):
        ys = np.stack([R.fourier_shift(y, N, h, v * (t - t_start)) for t in t_eval], axis=1)
        return dict(t=np.asarray(t_eval), y=ys, y_final=ys[:, -1], residual=0.0, success=True,
                    t_reached=float(t_eval[-1]), complete=True, retry=None)
    return ode


def test_end_to_end_fake_translation_and_registration_guard():
    tmp = tempfile.mkdtemp()
    old = (D.OUT, D.RUNS, D.SNAP, D.REGISTERED_FLAG)
    try:
        D.OUT, D.RUNS, D.SNAP = tmp, os.path.join(tmp, "runs"), os.path.join(tmp, "snap")
        D.REGISTERED_FLAG = os.path.join(tmp, "REGISTERED")
        try:
            D.stage_runs()
            raise AssertionError("runs started without registration")
        except D.Stop:
            pass
        st = D.start_states()
        rec = D.one_run("R0", st["R0"], D.e_pat(), ode=_fake_ode(-4e-7))
        assert len(rec["diag"]) == 201 and len(rec["track"]) == 2001 and len(rec["snapshots"]) == 11
        sp = D.speeds(rec["diag"])
        assert all(abs(sp[k] + 4e-7) < 1e-11 for k in "abv"), sp
        V = D.validity(rec)
        assert V["valid"], V
        assert D.c4_label(sp["a"], sp["b"], VF) == "PERSISTENT"
        lg = D.ledger()
        assert len(lg["integrations"]) == 1 and lg["integrations"][0]["wall_s"] is not None
        # the ledger counts at start: a 10th start is the last allowed
        lg["integrations"] = [dict(key=f"k{i}", wall_s=1.0) for i in range(10)]
        D._ledger_write(lg)
        try:
            D._ledger_start("x")
            raise AssertionError("cap not enforced")
        except D.Stop:
            pass
        # validity failure: a shape change beyond S_MAX after t_ref is caught (item 5)
        bad = json.loads(json.dumps(rec))
        bad["diag"][-1]["S_ref"] = 2e-3
        assert 5 in D.validity(bad)["failed"]
        bad = json.loads(json.dumps(rec))
        bad["diag"][50]["ambiguous"] = True
        assert 4 in D.validity(bad)["failed"]
    finally:
        D.OUT, D.RUNS, D.SNAP, D.REGISTERED_FLAG = old
        shutil.rmtree(tmp)


if __name__ == "__main__":
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_")]
    for k, f in tests:
        f()
        print("ok", k)
    print(len(tests), "tests passed")
