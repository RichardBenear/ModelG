"""MG9c driver (spec/MG9c_proposal.md draft 2, registered at 372336b). Deterministic only. Resumable at run level.

Stages (proposal §8 order of work):
    python src/run_mg9c.py q1        main fronts at 29.25, 29.30, 29.35 (replay of MG9 [0, 1800]) and 29.22
    python src/run_mg9c.py robust    robustness (a)-(c) of every SETTLED Q1 point
    python src/run_mg9c.py q2        upper settling bracket (only if 29.35 is SETTLED-ROBUST)
    python src/run_mg9c.py desc      descriptive width trials (d)-(e)
Frozen MG9 machinery is imported from run_mg9 unchanged (proposal §3). Prepared states come from the committed
results/MG9/states.npz. Every integration is recorded in results/MG9c/ledger.json; cap 60, cutoff 2 core-h.
"""
import os, sys, json, time, hashlib
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import landau_mg9 as L
import run_mg9 as R

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results", "MG9c")
RUNS = os.path.join(OUT, "runs")
SNAP = os.path.join(OUT, "snapshots")              # one .npy per snapshot; packed into snapshots.npz at the end
MG9_STATES = os.path.join(ROOT, "results", "MG9", "states.npz")
MG9_RUNS = os.path.join(ROOT, "results", "MG9", "runs")
RUN_CAP, CUTOFF_CORE_S = 60, 2 * 3600.0
LAM = R.LAM
RES_TOL = R.ARREST_RES_TOL                          # 1e-5
REPLAY_TOL = 1e-6
Q1_B = (29.25, 29.30, 29.35)
DESC_B = 29.22
B_INVADES = 29.40                                   # MG9 grid: INVADES
BR = (29.395313, 29.396094)                         # MG9 rechecked b_r bracket
Q2_TOL = 0.01

# Segment boundaries. Main / D01 / D02L runs: MG9's [0, 600], [600, 1800] replayed exactly, then the window
# boundaries; extension 12000 -> 15000 -> 30000. Perturbation runs: fixed 6000, no extension.
MAIN_SEGS = (0.0, 600.0, 1800.0, 3600.0, 6000.0, 12000.0)
EXT_SEGS = (12000.0, 15000.0, 30000.0)
MAIN_W = ((1800.0, 3600.0), (3600.0, 6000.0), (6000.0, 12000.0))
EXT_W = ((1800.0, 3600.0), (6000.0, 15000.0), (15000.0, 30000.0))
PERT_SEGS = (0.0, 1500.0, 3000.0, 6000.0)
PERT_W = ((0.0, 1500.0), (1500.0, 3000.0), (3000.0, 6000.0))
SNAP_FROM = 1800.0                                  # main runs: snapshots at every boundary >= 1800


class Stop(Exception):
    pass


# ------------------------------------------------------------------ ledger
def _ledger_path():
    return os.path.join(OUT, "ledger.json")


def ledger():
    try:
        with open(_ledger_path()) as f:
            return json.load(f)
    except FileNotFoundError:
        return dict(integrations=[], cap=RUN_CAP, cutoff_core_s=CUTOFF_CORE_S)


def _ledger_add(key, wall):
    lg = ledger()
    lg["integrations"].append(dict(key=key, wall_s=wall, at=time.strftime("%Y-%m-%dT%H:%M:%S")))
    os.makedirs(OUT, exist_ok=True)
    with open(_ledger_path(), "w") as f:
        json.dump(lg, f, indent=1)


def _ledger_check():
    lg = ledger()
    n, used = len(lg["integrations"]), sum(r["wall_s"] for r in lg["integrations"])
    if n + 1 > RUN_CAP:
        raise Stop(f"run cap {RUN_CAP} reached")
    if used >= CUTOFF_CORE_S:
        raise Stop(f"cutoff {CUTOFF_CORE_S / 3600:.1f} core-h reached")


# ------------------------------------------------------------------ records
def rec_path(key):
    return os.path.join(RUNS, key + ".json")


def load_rec(key):
    try:
        with open(rec_path(key)) as f:
            return json.load(f)
    except FileNotFoundError:
        return None


def save_rec(key, rec):
    os.makedirs(RUNS, exist_ok=True)
    with open(rec_path(key), "w") as f:
        json.dump(rec, f, indent=1, default=float)


def save_snap(name, y):
    os.makedirs(SNAP, exist_ok=True)
    np.save(os.path.join(SNAP, name + ".npy"), np.asarray(y, float))


def load_snap(name):
    return np.load(os.path.join(SNAP, name + ".npy"))


def sha(y):
    return hashlib.sha256(np.ascontiguousarray(y, dtype=np.float64).tobytes()).hexdigest()


def mg9_state(key):
    with np.load(MG9_STATES) as z:
        return np.array(z[key], float)


def mg9_front_track(b):
    with open(os.path.join(MG9_RUNS, f"front_grid_{b:.4f}.json")) as f:
        return json.load(f)


def platform():
    import platform as pf, scipy
    return dict(python=pf.python_version(), numpy=np.__version__, scipy=scipy.__version__, system=pf.platform())


# ------------------------------------------------------------------ §4 classification
def _interp(track, t, key):
    ts = np.array([r["t"] for r in track])
    xs = np.array([r[key] for r in track], float)
    return float(np.interp(t, ts, xs))


def late_classify(track, T, residual, complete, windows, lam=LAM, t_struct=100.0):
    """Proposal §4, rules 1-5; rule 6 returns 'SLOW' (the caller extends or maps it to UNRESOLVED-SLOW).
    track: list of dict(t, status, xL, xR, n, breakup, remote) every 5 time units from t = 0."""
    if not complete or not track or track[0]["status"] != "ok":
        return "UNRESOLVED-NUMERICAL", {}
    st = [r["status"] for r in track]
    first = next((i for i, s in enumerate(st) if s in ("fill", "loss")), None)       # rule 2: from t = 0
    if first is not None:
        ev = st[first]
        info = dict(event=ev, t_event=track[first]["t"])
        if all(s == ev for s in st[first:]):
            return ("INVADES" if ev == "fill" else "RETREATS"), info
        return "UNRESOLVED-TEMPORARY", info
    if any(s == "untracked" for s in st):                                              # rule 3 (as MG9)
        return "UNRESOLVED-INTERNAL", dict(reason="untracked")
    w = [r for r in track if r["t"] >= t_struct - 1e-9]
    if any(r["breakup"] or r["remote"] for r in w):
        return "UNRESOLVED-INTERNAL", dict(reason="breakup/remote")
    oL = -(np.array([r["xL"] for r in w]) - w[0]["xL"])
    oR = np.array([r["xR"] for r in w]) - w[0]["xR"]
    if R.has_reversal(oL, lam / 2) or R.has_reversal(oR, lam / 2):
        return "UNRESOLVED-REVERSAL", {}
    # windows: outward displacement D and mean speed s per interface
    D = {"L": [], "R": []}
    for a, b in windows:
        D["L"].append(-(_interp(track, b, "xL") - _interp(track, a, "xL")))
        D["R"].append(_interp(track, b, "xR") - _interp(track, a, "xR"))
    s = {k: [d / (b - a) for d, (a, b) in zip(D[k], windows)] for k in D}
    dL, dR = D["L"][2], D["R"][2]
    info = dict(D=D, s=s, residual=residual)
    if dL >= lam / 4 and dR >= lam / 4:                                                # rule 4
        return "ADVANCES", info
    if dL <= -lam / 4 and dR <= -lam / 4:
        return "RECEDES", info
    if max(abs(dL), abs(dR)) >= lam / 4:
        return "UNRESOLVED-ASYMMETRIC", info
    t = np.array([r["t"] for r in track])
    last = t >= T - 1000.0 - 1e-9
    xl, xr = np.array([r["xL"] for r in track]), np.array([r["xR"] for r in track])
    vL = float(np.polyfit(t[last], -xl[last], 1)[0])
    vR = float(np.polyfit(t[last], xr[last], 1)[0])
    n = np.array([r["n"] for r in track if r["t"] >= windows[0][0] - 1e-9])
    dec = {k: (abs(s[k][2]) < abs(s[k][1]) < abs(s[k][0])) or abs(D[k][2]) <= lam / 400 for k in D}
    info.update(vL=vL, vR=vR, n_range=int(n.max() - n.min()), decreasing=dec)
    if (abs(dL) <= lam / 40 and abs(dR) <= lam / 40 and abs(vL) <= lam / 40000 and abs(vR) <= lam / 40000
            and dec["L"] and dec["R"] and residual <= RES_TOL and n.max() - n.min() <= 1):     # rule 5
        return "SETTLED", info
    return "SLOW", info


# ------------------------------------------------------------------ long front runs
def long_front(key, b, N, h, y0, E_pat, *, kind="main", snap_tag=None):
    """Integrate y0 through the segment boundaries of `kind` ('main' or 'pert'), sampling the tracker every 5 units,
    saving full-field snapshots at the registered boundaries; classify by §4; extend main runs once on rule 6."""
    rec = load_rec(key)
    if rec is not None:
        return rec
    _ledger_check()
    t0w = time.time()
    tr = R.Tracker(N, h, E_pat)
    track = [dict(t=0.0, **tr.sample(y0[N:2 * N]))]
    snaps = []
    snap_tag = snap_tag or key

    def run_segs(segs, y, snap_from):
        r = None
        for ta, tb in zip(segs, segs[1:]):
            te = np.arange(ta + 5.0, tb + 1e-9, 5.0)
            r = R.mg9_ode(b, N, h, tb, y, seed_scale=0.0, t_eval=te, t_start=ta)
            for j, tj in enumerate(r["t"]):
                track.append(dict(t=float(tj), **tr.sample(r["y"][N:2 * N, j])))
            if not r["complete"]:
                return r
            y = r["y_final"]
            if tb >= snap_from - 1e-9:
                name = f"{snap_tag}_t{int(round(tb))}"
                save_snap(name, y)
                snaps.append(dict(name=name, t=tb, sha256=sha(y)))
        return r

    if kind == "main":
        r = run_segs(MAIN_SEGS, y0, SNAP_FROM)
        T, extended = MAIN_SEGS[-1], False
        label, info = late_classify(track, T, r["residual"], r["complete"], MAIN_W)
        if r["complete"] and label == "SLOW":
            r = run_segs(EXT_SEGS, r["y_final"], SNAP_FROM)
            T, extended = EXT_SEGS[-1], True
            label, info = late_classify(track, T, r["residual"], r["complete"], EXT_W)
    elif kind == "pert":
        r = run_segs(PERT_SEGS, y0, PERT_SEGS[1])
        T, extended = PERT_SEGS[-1], False
        label, info = late_classify(track, T, r["residual"], r["complete"], PERT_W)
    else:
        raise ValueError(kind)
    if label == "SLOW":
        label = "UNRESOLVED-SLOW"
    _ledger_add(key, time.time() - t0w)
    rec = dict(key=key, b=b, N=N, h=h, kind=kind, E_pat=E_pat, init_sha256=sha(y0), T=T, extended=extended,
               label=label, info=info, complete=r["complete"], success=r["success"], t_reached=r["t_reached"],
               residual=r["residual"], retry=r["retry"], seed_scale=0.0, snapshots=snaps, track=track,
               platform=platform())
    save_rec(key, rec)
    print(key, label, flush=True)
    return rec


def replay_check(track, mg9_track, tol=REPLAY_TOL):
    """§3: the new track over [0, 1800] must reproduce MG9's to tol in xL and xR (and the same statuses)."""
    old = {round(r["t"], 6): r for r in mg9_track}
    worst, status_ok, n = 0.0, True, 0
    for r in track:
        if r["t"] > 1800.0 + 1e-9:
            break
        o = old.get(round(r["t"], 6))
        if o is None:
            status_ok = False
            continue
        n += 1
        if r["status"] != o["status"]:
            status_ok = False
        if r["xL"] is not None and o["xL"] is not None:
            worst = max(worst, abs(r["xL"] - o["xL"]), abs(r["xR"] - o["xR"]))
    ok = status_ok and n == len(mg9_track) and worst <= tol
    return dict(ok=bool(ok), max_abs_diff=worst, samples=n, mg9_samples=len(mg9_track), status_ok=status_ok)


# ------------------------------------------------------------------ preparations
def d02_prepared(b):
    """Q1/desc points: MG9's prepared state, as used by the MG9 grid run (or the sweep state for 29.22)."""
    if abs(b - DESC_B) < 1e-9:
        return mg9_state("sweep_29.22"), "sweep_29.22"
    prov = mg9_front_track(b)["prep"]
    return mg9_state(prov), prov


def d02_prepare_mid(b):
    """Q2 midpoints: MG9 §6.1, carry from the nearest sweep state above, settle to 1000 (extend to 3000)."""
    key = f"prep_{b:.5f}"
    rec = load_rec(key)
    d = L.DOMAINS["D02"]
    N, h = d["N"], d["h"]
    if rec is None:
        above = min(s for s in R.SWEEP_B if s > b + 1e-9)
        src = f"sweep_{above:.2f}"
        y0 = R.carry_concentration(mg9_state(src), N, above, b)
        _ledger_check()
        t0w = time.time()
        r = R.mg9_ode(b, N, h, 1000.0, y0, seed_scale=0.0)
        ext = False
        if r["complete"] and r["residual"] > R.STAT_TOL:
            r = R.mg9_ode(b, N, h, 3000.0, r["y_final"], seed_scale=0.0, t_start=1000.0)
            ext = True
        _ledger_add(key, time.time() - t0w)
        status = R.prep_classify(r, N, h, None)
        rec = dict(key=key, b=b, source=src, status=status, extended=ext, complete=r["complete"],
                   residual=r["residual"], carry="concentration", seed_scale=0.0)
        if r["complete"]:
            save_snap(key, r["y_final"])
            rec["sha256"] = sha(r["y_final"])
        save_rec(key, rec)
    if rec["status"] != "PERSISTS":
        return None, rec
    return load_snap(key), rec


def d01_prepare(b, y, tag):
    """Robustness (a): Fourier-interpolate the D02 state to D01, settle to 1000 (extension 3000), PERSISTS check."""
    d, d1 = L.DOMAINS["D02"], L.DOMAINS["D01"]
    key = f"prep01_{tag}"
    rec = load_rec(key)
    if rec is None:
        y1 = R.fourier_resample(y, d["N"], d1["N"])
        _ledger_check()
        t0w = time.time()
        r = R.mg9_ode(b, d1["N"], d1["h"], 1000.0, y1, seed_scale=0.0)
        ext = False
        if r["complete"] and r["residual"] > R.STAT_TOL:
            r = R.mg9_ode(b, d1["N"], d1["h"], 3000.0, r["y_final"], seed_scale=0.0, t_start=1000.0)
            ext = True
        _ledger_add(key, time.time() - t0w)
        status = R.prep_classify(r, d1["N"], d1["h"], None)
        rec = dict(key=key, b=b, status=status, extended=ext, complete=r["complete"], residual=r["residual"])
        if r["complete"]:
            save_snap(key, r["y_final"])
            rec["sha256"] = sha(r["y_final"])
        save_rec(key, rec)
    return (load_snap(key) if rec["status"] == "PERSISTS" else None), rec


def d02l_prepare(b, y, tag):
    """Robustness (b): tile twice, settle on D02L to 300; PERSISTS = complete, residual <= 1e-6, 68 cells (MG9)."""
    d, dl = L.DOMAINS["D02"], L.DOMAINS["D02L"]
    N, NL, h = d["N"], dl["N"], dl["h"]
    key = f"prepL_{tag}"
    rec = load_rec(key)
    if rec is None:
        tiled = np.concatenate([np.tile(y[c * N:(c + 1) * N], 2) for c in range(3)])
        _ledger_check()
        t0w = time.time()
        r = R.mg9_ode(b, NL, h, 300.0, tiled, seed_scale=0.0)
        _ledger_add(key, time.time() - t0w)
        ok = r["complete"] and r["residual"] <= R.STAT_TOL
        n = None
        if r["complete"]:
            pX = r["y_final"][NL:2 * NL]
            n = R.cell_count(pX, h, float(R.envelope(pX, h).mean())) if np.abs(pX).max() >= 0.1 else 0
        status = ("PERSISTS" if ok and n == 68 else "PREP-UNRESOLVED" if not ok
                  else "PREP-LOST" if n == 0 else "PREP-CHANGED")
        rec = dict(key=key, b=b, status=status, complete=r["complete"], residual=r["residual"], cells=n)
        if r["complete"]:
            save_snap(key, r["y_final"])
            rec["sha256"] = sha(r["y_final"])
        save_rec(key, rec)
    return (load_snap(key) if rec["status"] == "PERSISTS" else None), rec


def e_pat(y, N, h):
    return float(R.envelope(y[N:2 * N], h).mean())


def patch(y, N, h, half_width):
    return R.splice_patch(R.phase_convention(y, N, h), N, h, half_width)


# ------------------------------------------------------------------ §5 trial and point status
def trial_status(kind, rec, ref_edges=None, lam=LAM):
    """PASS / FAIL / UNRESOLVED for one robustness trial (§5). rec is a front record or a prep-failure record."""
    if rec is None:
        return "UNRESOLVED"
    lab = rec.get("label", "UNRESOLVED-PREP")
    if lab in ("ADVANCES", "RECEDES", "INVADES", "RETREATS"):
        return "FAIL"
    if lab != "SETTLED":
        return "UNRESOLVED"
    if kind == "pert":
        xl, xr = rec["track"][-1]["xL"], rec["track"][-1]["xR"]
        if max(abs(xl - ref_edges[0]), abs(xr - ref_edges[1])) > lam / 10:
            return "FAIL"
    return "PASS"


def point_status(statuses):
    """statuses: dict trial -> PASS/FAIL/UNRESOLVED for (a), (b), (c+), (c-)."""
    v = list(statuses.values())
    if all(s == "PASS" for s in v):
        return "SETTLED-ROBUST"
    if any(s == "FAIL" for s in v):
        return "NOT-ROBUST"
    return "ROBUSTNESS-UNRESOLVED"


# ------------------------------------------------------------------ §6 bracket comparison, §7 outcome
def q2_compare(b1, b2, br=BR, tol=Q2_TOL):
    """Interval rule: difference in [b1 - r_hi, b2 - r_lo]."""
    lo, hi = round(b1 - br[1], 9), round(b2 - br[0], 9)            # rounded so the inclusive limits are exact
    if lo >= -tol and hi <= tol:
        return "NEAR b_r"
    if hi < -tol:
        return "BELOW b_r"
    if lo > tol:
        return "ABOVE b_r"
    return "UNRESOLVED RELATIVE TO THE TOLERANCE"


def q1_outcome(points):
    """points: dict b -> dict(label, replay_ok (None if exempt), status (point status or None))."""
    for p in points.values():
        if p["label"].startswith("UNRESOLVED") or p.get("replay_ok") is False:
            return "UNRESOLVED"
        if p["label"] == "SETTLED" and p.get("status") == "ROBUSTNESS-UNRESOLVED":
            return "UNRESOLVED"
    if all(p["label"] in ("ADVANCES", "INVADES") for p in points.values()):
        return "SLOW INVASION"
    if all(p["label"] == "SETTLED" and p.get("status") == "SETTLED-ROBUST" for p in points.values()):
        return "BOUNDED"
    return "MIXED"


# ------------------------------------------------------------------ stages
def main_key(b):
    return f"front_D02_{b:.5f}"


def stage_q1():
    d = L.DOMAINS["D02"]
    N, h = d["N"], d["h"]
    for b in Q1_B + (DESC_B,):
        y, prov = d02_prepared(b)
        E = e_pat(y, N, h)
        rec = long_front(main_key(b), b, N, h, patch(y, N, h, 25.0), E)
        if "replay" not in rec:
            if b in Q1_B:
                old = mg9_front_track(b)
                rec["replay"] = replay_check(rec["track"], old["track"])
                rec["replay"]["E_pat_diff"] = abs(E - old["E_pat"])
            else:
                rec["replay"] = dict(ok=None, note="exempt: no MG9 front track at 29.22")
            rec["prep"] = prov
            save_rec(main_key(b), rec)
            print(main_key(b), "replay", rec["replay"], flush=True)


def robustness(b, tag=None):
    """Trials (a)-(c) for a SETTLED main run at b; returns dict of statuses and records' keys."""
    d, d1, dl = L.DOMAINS["D02"], L.DOMAINS["D01"], L.DOMAINS["D02L"]
    N, h = d["N"], d["h"]
    tag = tag or f"{b:.5f}"
    main = load_rec(main_key(b))
    y = d02_prepared(b)[0] if b in Q1_B + (DESC_B,) else d02_prepare_mid(b)[0]
    out = {}
    # (a) D01
    y1, prec = d01_prepare(b, y, tag)
    if y1 is None:
        out["a"] = dict(status="UNRESOLVED", prep=prec["status"])
    else:
        r = long_front(f"front_D01_{tag}", b, d1["N"], d1["h"], patch(y1, d1["N"], d1["h"], 25.0),
                       e_pat(y1, d1["N"], d1["h"]))
        out["a"] = dict(status=trial_status("a", r), label=r["label"])
    # (b) D02L, half-width 25
    yl, prec = d02l_prepare(b, y, tag)
    if yl is None:
        out["b"] = dict(status="UNRESOLVED", prep=prec["status"])
    else:
        r = long_front(f"front_D02L25_{tag}", b, dl["N"], dl["h"], patch(yl, dl["N"], dl["h"], 25.0),
                       e_pat(yl, dl["N"], dl["h"]))
        out["b"] = dict(status=trial_status("b", r), label=r["label"])
    # (c) perturbations of the settled final field of the main run
    yT = load_snap(f"{main['key']}_t{int(round(main['T']))}")
    ref = (main["track"][-1]["xL"], main["track"][-1]["xR"])
    for name, f in (("c+", 1.05), ("c-", 0.95)):
        r = long_front(f"pert_{'p' if f > 1 else 'm'}05_{tag}", b, N, h, f * yT, main["E_pat"], kind="pert")
        out[name] = dict(status=trial_status("pert", r, ref), label=r["label"],
                         edge_shift=[r["track"][-1]["xL"] - ref[0], r["track"][-1]["xR"] - ref[1]])
    status = point_status({k: v["status"] for k, v in out.items()})
    res = dict(b=b, trials=out, status=status)
    save_rec(f"robust_{tag}", res)
    print(f"robust {tag}", status, {k: v["status"] for k, v in out.items()}, flush=True)
    return res


def stage_robust():
    for b in Q1_B + (DESC_B,):
        rec = load_rec(main_key(b))
        if rec is None:
            raise RuntimeError("run stage q1 first")
        if rec["label"] == "SETTLED" and (b == DESC_B or rec["replay"]["ok"]):
            if load_rec(f"robust_{b:.5f}") is None:
                robustness(b)


def stage_q2():
    r35 = load_rec(f"robust_{29.35:.5f}")
    if r35 is None or r35["status"] != "SETTLED-ROBUST":
        save_rec("q2", dict(eligible=False, reason="29.35 not SETTLED-ROBUST"))
        print("Q2 not eligible", flush=True)
        return
    d = L.DOMAINS["D02"]
    N, h = d["N"], d["h"]
    lo, hi = 29.35, B_INVADES
    samples, status = [], "resolved"
    for _ in range(3):
        mid = round((lo + hi) / 2, 6)
        y, prec = d02_prepare_mid(mid)
        if y is None:
            samples.append(dict(b=mid, prep=prec["status"], label=None))
            status = "unresolved (preparation)"
            break
        rec = long_front(main_key(mid), mid, N, h, patch(y, N, h, 25.0), e_pat(y, N, h))
        samples.append(dict(b=mid, prep="PERSISTS", label=rec["label"]))
        if rec["label"] == "SETTLED":
            lo = mid
        elif rec["label"] in ("INVADES", "ADVANCES"):
            hi = mid
        else:
            status = f"unresolved ({rec['label']})"
            break
    final = dict(eligible=True, samples=samples, bracket=[lo, hi], status=status)
    if status == "resolved" and lo != 29.35:
        rb = robustness(lo, tag=f"{lo:.5f}")
        final["final_point_robustness"] = rb["status"]
        if rb["status"] != "SETTLED-ROBUST":
            final["status"] = "unresolved (final point not SETTLED-ROBUST)"
    final["comparison"] = q2_compare(lo, hi) if final["status"] == "resolved" else "NOT EVALUATED"
    save_rec("q2", final)
    print("Q2", final, flush=True)


def stage_desc():
    d, dl = L.DOMAINS["D02"], L.DOMAINS["D02L"]
    N, h = d["N"], d["h"]
    for b in Q1_B + (DESC_B,):
        main = load_rec(main_key(b))
        if main is None or main["label"] != "SETTLED":
            continue
        tag = f"{b:.5f}"
        yl, prec = d02l_prepare(b, d02_prepared(b)[0], tag)
        if yl is not None:
            long_front(f"front_D02L50_{tag}", b, dl["N"], dl["h"], patch(yl, dl["N"], dl["h"], 50.0),
                       e_pat(yl, dl["N"], dl["h"]))
        y = d02_prepared(b)[0]
        long_front(f"front_D02w_{tag}", b, N, h, patch(y, N, h, 25.0 + LAM / 4), e_pat(y, N, h))


def pack_snapshots():
    """Collect every snapshot into results/MG9c/snapshots.npz (committed) with a SHA-256 manifest."""
    if not os.path.isdir(SNAP):
        return
    names = sorted(f[:-4] for f in os.listdir(SNAP) if f.endswith(".npy"))
    arrs = {n: load_snap(n) for n in names}
    np.savez_compressed(os.path.join(OUT, "snapshots.npz"), **arrs)
    with open(os.path.join(OUT, "snapshots_manifest.json"), "w") as f:
        json.dump({n: sha(a) for n, a in arrs.items()}, f, indent=1)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = dict(q1=stage_q1, robust=stage_robust, q2=stage_q2, desc=stage_desc, pack=pack_snapshots).get(stage)
    if fn is None:
        print(__doc__)
        sys.exit(1)
    try:
        fn()
    except Stop as e:
        print("STOP:", e, flush=True)
        sys.exit(2)
