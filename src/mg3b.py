"""MG3b library (GATES.md MG3b; spec/MG3b_followup_proposal.md, registered at 202c00f; Part 3 decisions in
results/MG3b/calibration/decisions.json).

Simulation, detection, tracking and confirmation are MG3's (`mg3.sim_run`); this module adds what MG3b changes:
- summarize_b: the MG3 run record plus EVERY confirmation field (float32, 3 x N) with its class and time;
- the stop20 rule for Part 2 cells, the Part 1 / Part 1b decisions and the production branch;
- the budget ledger and the per-batch ceiling check (enforcement rules 1-3).
"""
import os, csv, time, platform
import numpy as np
import mg3
from calib_mg3b import part1_decision

TAU0 = 6.666666666666667e-4
GATE = 4
OMEGAS = (130.0, 120.0, 110.0, 100.0)            # om_idx 0..3
STAGE = dict(P1_tau=0, P1_tau2=1, P1_tau4=2, P2=3, TIMING=6, PILOT=7)
BATCH, TARGET, CAP = 16, 20, 64
OUTCOME = {0: "EQUIVALENT", 1: "DIFFERENT", 2: "INCONCLUSIVE"}

# Registered core-hour ceilings (GATES.md MG3b) and prior projections per run (core-s; proposal tables / run counts).
CEIL_H = dict(P1_tau=6.0, P1_tau2=20.0, P1_tau4=40.0, P2_tau=12.0, P2_tau2=30.0)
CEIL_H_96 = dict(P1_tau=9.0, P1_tau2=30.0)
PRIOR_S = dict(P1_tau=3.8 * 3600 / 64, P1_tau2=13.0 * 3600 / 64, P1_tau4=26.7 * 3600 / 64,
               P2_tau=7.2 * 3600 / 117, P2_tau2=19.4 * 3600 / 117)
TAU4_ASSUMED = 5.0                                # core-s per time unit assumed for tau/4 in the registration
SAFETY = 1.25


def om_idx(Om):
    return OMEGAS.index(float(Om))


# ------------------------------------------------------------------------------------------------ run records
def summarize_b(st):
    """MG3 record + all confirmation fields (protocol: every persistent / excursion / unresolved field is saved)."""
    rec = mg3.summarize(st)
    times = np.array(st["times"])
    rec["conf_fields"] = {int(sc): dict(field=np.asarray(st["phi"][sc], np.float32), t_conf=float(times[sc]),
                                        **{k: v for k, v in cc.items()})
                          for sc, cc in st["conf_cache"].items()}
    return rec


def stop20_done(runs):
    """Part 2 cell: complete when, at the end of a batch, counted events >= 20, or at the cap of 64 runs.
    runs: finished records (with 'rep'). A batch counts only when all its replicates are finished.
    Returns (done, n_runs_needed)."""
    have = {r["rep"]: r for r in runs}
    for b in range(1, CAP // BATCH + 1):
        reps = range(b * BATCH)
        if not all(i in have for i in reps):
            return False, b * BATCH
        ev = sum(have[i]["event"] for i in reps if have[i]["complete"])
        if ev >= TARGET or b * BATCH >= CAP:
            return True, b * BATCH
    return True, CAP


def arm_counts(runs):
    r = mg3.cell_rate(runs)
    return r["n"], r["exposure"]


def compare(method, runs_a, runs_b):
    """Rate ratio k(b)/k(a) decision (Part 1: a = tau, b = tau/2; Part 1b: a = tau/2, b = tau/4)."""
    n1, E1 = arm_counts(runs_a)
    n2, E2 = arm_counts(runs_b)
    d = part1_decision(method, np.array([n1]), np.array([float(E1)]), np.array([n2]), np.array([float(E2)]))[0]
    return OUTCOME[int(d)], dict(n_a=n1, E_a=E1, n_b=n2, E_b=E2)


def production_branch(out1, out1b=None):
    """Returns (branch, tau, omega120_cell) per GATES.md MG3b; branch 'stop' means MG3b stops (time-step grounds)."""
    if out1 == "EQUIVALENT":
        return "tau", TAU0, "P1_tau"
    if out1 == "INCONCLUSIVE":
        return "tau2_inc", TAU0 / 2, "P1_tau2"
    if out1b is None:
        return "part1b_needed", None, None
    if out1b == "EQUIVALENT":
        return "tau2_1b", TAU0 / 2, "P1_tau2"
    return "stop", None, None


# ------------------------------------------------------------------------------------------------ ledger and budget
LEDGER_COLS = ["item", "key", "kind", "wall_s", "workers", "machine", "os", "python", "numpy", "scipy", "numba",
               "utc_time"]


def _versions():
    import numpy, scipy
    try:
        import numba
        nbv = numba.__version__
    except ImportError:
        nbv = "n/a"
    return dict(machine=platform.node(), os=platform.platform(), python=platform.python_version(),
                numpy=numpy.__version__, scipy=scipy.__version__, numba=nbv)


def ledger_append(path, item, key, kind, wall_s, workers):
    new = not os.path.exists(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    v = _versions()
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(LEDGER_COLS)
        w.writerow([item, key, kind, f"{wall_s:.1f}", workers, v["machine"], v["os"], v["python"], v["numpy"],
                    v["scipy"], v["numba"], time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())])


def ledger_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def ledger_used_s(path, item=None, kinds=None):
    return sum(float(r["wall_s"]) for r in ledger_rows(path)
               if (item is None or r["item"] == item) and (kinds is None or r["kind"] in kinds))


def ceiling_h(item, n_arm=64, tau4_factor=1.0):
    # P1_tau4: the registered rule scaled 40 h by the measured tau/4 cost / 5.0. PI decision 2026-10-30 (before the
    # Part 1 result; GATES.md MG3b): the timing ran 2 runs alone and underestimates loaded cost, so the ceiling is the
    # unscaled 40 core-hours. tau4_factor still scales the prior per-run projection.
    c = CEIL_H_96.get(item, CEIL_H[item]) if n_arm == 96 else CEIL_H[item]
    return c


def batch_allowed(used_s, completed_walls, batch_size, item, n_arm=64, tau4_factor=1.0):
    """Enforcement rule 2: submit only if used + projected batch cost <= the item's ceiling.
    Projection = batch size x (mean cost per completed run in the item, or the prior projection) x 1.25."""
    per = float(np.mean(completed_walls)) if len(completed_walls) else PRIOR_S[item] * (
        tau4_factor if item == "P1_tau4" else 1.0)
    proj = batch_size * per * SAFETY
    ceil_s = ceiling_h(item, n_arm, tau4_factor) * 3600
    return used_s + proj <= ceil_s, dict(used_h=used_s / 3600, projected_h=proj / 3600, ceiling_h=ceil_s / 3600)
