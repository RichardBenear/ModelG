"""Regression tests for measure.track / analyze_homogeneous (MG005, ChatGPT finding P1: first homogeneous birth omitted).
Run: python tests/test_track.py   (or pytest)
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from measure import track, analyze_homogeneous
from runner import reference

L = 100.0
TIMES = np.arange(0.0, 100.0 + 1e-9, 0.5)


def dets_from(spans):
    """spans: list of (x, t_on, t_off) -> per-snapshot detection arrays."""
    return [np.array(sorted(x for (x, a, b) in spans if a - 1e-9 <= t <= b + 1e-9)) for t in TIMES]


def kinds(ev):
    return sorted((e[0], round(float(e[1]), 2)) for e in ev)


def test_homogeneous_single_birth():
    # ChatGPT's regression input: one detection at x = 0 from t = 0.5 to 30
    tr = track(dets_from([(0.0, 0.5, 30.0)]), TIMES, L, x0=None)
    assert kinds(tr["events"]) == [("birth", 0.5)], tr["events"]
    assert np.all(tr["orig_hit"] == -1)


def test_homogeneous_two_births():
    tr = track(dets_from([(0.0, 5.0, 60.0), (30.0, 40.0, 100.0)]), TIMES, L, x0=None)
    assert kinds(tr["events"]) == [("birth", 5.0), ("birth", 40.0)], tr["events"]


def test_homogeneous_short_track_is_not_birth():
    tr = track(dets_from([(0.0, 5.0, 20.0)]), TIMES, L, x0=None)   # duration 15 < 20
    assert tr["events"] == [], tr["events"]


def test_seeded_original_no_events():
    tr = track(dets_from([(0.0, 0.0, 100.0)]), TIMES, L, x0=0.0)
    assert tr["events"] == [] and np.all(tr["orig_hit"] == 0)


def test_seeded_loss():
    tr = track(dets_from([(0.0, 0.0, 40.0)]), TIMES, L, x0=0.0)
    assert kinds(tr["events"]) == [("loss", 40.0)], tr["events"]


def test_seeded_plus_birth():
    tr = track(dets_from([(0.0, 0.0, 100.0), (25.0, 30.0, 60.0)]), TIMES, L, x0=0.0)
    assert kinds(tr["events"]) == [("birth", 30.0)], tr["events"]


def _field_with_dips(spans, Om=1e8):
    """Counts arrays (snap, N) with the b = 29 reference soliton placed at each active span."""
    ref = reference(29.0, 0.2)
    x, S0, phi = ref["x"], ref["S0"], ref["phi"]
    N = x.size
    F = np.zeros((len(TIMES), 3, N))
    for s, t in enumerate(TIMES):
        for (x0, a, b) in spans:
            if a - 1e-9 <= t <= b + 1e-9:
                F[s] += np.roll(phi, int(round(x0 / 0.2)), axis=1)
    n = [Om * (S0[k] + F[:, k, :]) for k in range(3)]
    return n, ref


def test_analyze_homogeneous_first_birth_and_burn_in():
    (nG, nX, nY), ref = _field_with_dips([(0.0, 10.0, 60.0), (40.0, 70.0, 100.0)])
    a = analyze_homogeneous(nG, nX, nY, TIMES, 1e8, 0.2, ref["S0"], ref["A_s"], ref["lam"], t_start=0.0)
    assert a["births"] == [10.0, 70.0], a["births"]
    assert a["first_birth"] == 10.0
    b = analyze_homogeneous(nG, nX, nY, TIMES, 1e8, 0.2, ref["S0"], ref["A_s"], ref["lam"], t_start=50.0)
    assert b["first_birth"] == 70.0, b["first_birth"]          # burn-in excludes the t = 10 birth


def test_analyze_homogeneous_null():
    (nG, nX, nY), ref = _field_with_dips([])
    a = analyze_homogeneous(nG, nX, nY, TIMES, 1e8, 0.2, ref["S0"], ref["A_s"], ref["lam"])
    assert a["births"] == [] and a["first_birth"] is None and max(a["n_dets"]) == 0


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for f in fns:
        f(); print("ok", f.__name__)
    print(f"{len(fns)} tests passed")
