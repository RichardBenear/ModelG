"""MG5 driver (spec/MG5_proposal.md draft 3, frozen at e2a25e3; GATES.md MG5 and its pre-run implementation
specification). A settled particle on an evolving periodic Delaunay graph.

    $env:MODELG_WORKERS=9; python src/run_mg5.py part0       (PowerShell, PI's PC): Part 0 checks (unregistered graph)
    $env:MODELG_WORKERS=9; python src/run_mg5.py             Parts A and B (needs results/MG5/part0/decisions.json)

Resumable at run level. Committed outputs (results/MG5/): runs/<key>.json, checks/<key>.npz, budget_ledger.csv,
local_manifest.csv, driver_log.json. PC-only outputs (results/MG5/local/, gitignored): event lists, full-field snapshots
and final states, listed with SHA-256 in local_manifest.csv.
"""
import os, sys, json, time, platform, csv, hashlib, subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(__file__))
import graphs as gr
import graph_ode as go
import evolve as ev
from lattice import P17, steady_state
from runner import seed_for

ROOT = os.path.join(os.path.dirname(__file__), "..", "results")
OUT = os.environ.get("MG5_OUT") or os.path.join(ROOT, "MG5")
LOCAL = os.path.join(OUT, "local")
DEC4 = os.path.join(ROOT, "MG4", "part0", "decisions.json")
MG4C_FIELDS = os.path.join(ROOT, "MG4c", "fields")
P0DEC = os.path.join(OUT, "part0", "decisions.json")
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))

THETA_28, THETA_29 = -4.457581675528346, -5.117354421620525      # MG4c T reference at b = 28 (1/2 x -8.915163); MG4
N0 = 18901
RUN_CAP, CUTOFF = 24, 40.0
CHECK_EVERY, SNAP_EVERY = 10.0, 50.0
DRUNS = [(0, 0), (0, 1), (1, 0), (1, 1)]
RATES = [0.0, 0.002, 0.02, 0.2]
FAILED, COMPLETE, INCOMPLETE = "FAILED-NUMERICAL", "COMPLETE", "INCOMPLETE"


def code_revision():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=os.path.dirname(__file__),
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def platform_info():
    import scipy
    return dict(machine=platform.node(), os=platform.platform(), python=platform.python_version(),
                numpy=np.__version__, scipy=scipy.__version__)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def is_multiple(t, step):
    return abs(t / step - round(t / step)) < 1e-9


# ------------------------------------------------------------------ the simulation engine (registered and Part 0)
def run_ticks(pos0, F0, h, P, theta, nu, events, delta, dt, t_end, Lx, Ly, seed_audit=0, snap=None, log=None):
    """The registered tick loop (proposal §3.2): for k = 1..K, integrate on M_(k-1) to t_k; measure a scheduled check on
    M_(k-1) and the boundary field; apply tick k's events (births interpolated on M_(k-1)); rebuild M_k; the next
    Jacobian is formed at the edited state. Returns (result, checks, ticks)."""
    t_wall = time.time()
    mesh = ev.Mesh(pos0, Lx, Ly)
    ids = np.arange(len(pos0), dtype=np.int64)
    F = np.asarray(F0, float).copy()
    stepper = ev.Ros2(mesh.g, h, P, dt)
    K = ev.n_steps(t_end, delta)
    nsub = ev.n_steps(delta, dt)
    slices = ev.tick_slices(events["t"], delta, t_end) if events is not None else None
    rng = np.random.default_rng(seed_audit)
    s0 = ev.smoothed_X(mesh.g, F)
    checks = [dict(t=0.0, rec=ev.sparse_record(mesh.pos, ids, s0, Lx, Ly),
                   parts=ev.detect_full(mesh.pos, s0, theta, Lx, Ly))]
    ticks = {k: [] for k in ("t", "n_prop_births", "n_prop_deaths", "n_eff_births", "n_eff_deaths", "n_ephemeral",
                             "jump_resid", "audit_max", "N", "I_before", "R", "B", "I_after")}
    flags, status, note = [], COMPLETE, None
    n_factor = 0
    y = F.ravel()
    k_done = 0
    for k in range(1, K + 1):
        tk = k * delta
        stepper.factor(y)
        n_factor += 1
        y = stepper.advance(y, nsub)
        if not np.all(np.isfinite(y)) or np.abs(y).max() > 1e3:
            status, note = FAILED, f"non-finite or |phi| > 1e3 at t = {tk:.4f}"
            break
        N = len(ids)
        bal = go.balances(y, stepper.rhs(0.0, y), N, P)
        if not bal["ok"]:
            hp = ev.hp_balance(mesh.g, h, P, y)
            flags.append(dict(kind="balance", t=tk, r1=bal["r1"], r2=bal["r2"], recheck=hp))
            if not hp["resolved"]:
                status, note = FAILED, f"unresolved balance flag at t = {tk:.4f}"
                break
        F = y.reshape(3, N)
        if is_multiple(tk, CHECK_EVERY):                     # check on M_(k-1) and the boundary field, before edits
            s = ev.smoothed_X(mesh.g, F)
            checks.append(dict(t=tk, rec=ev.sparse_record(mesh.pos, ids, s, Lx, Ly),
                               parts=ev.detect_full(mesh.pos, s, theta, Lx, Ly)))
            if log:
                log(f"t = {tk:g}: {len(checks[-1]['parts'])} particle(s), N = {N}")
        if snap is not None and (is_multiple(tk, SNAP_EVERY) or k == K):
            snap(tk, mesh.pos, ids, F)
        k_done = k
        if events is None or k == K:                         # no edits after the final boundary
            continue
        a, b = slices[k - 1]                                 # tick k's events, (t_(k-1), t_k]
        E = events[a:b]
        pos_n, ids_n, F_n, acct, births = ev.apply_tick(mesh, F, ids, E)
        if not acct["jump_ok"]:
            dead_rows = np.flatnonzero(~np.isin(ids, ids_n))
            ex = ev.exact_jump(F, dead_rows, births["val"], F_n)
            flags.append(dict(kind="jump", t=tk, resid=acct["jump_resid"], recheck=ex))
            if not ex["resolved"]:
                status, note = FAILED, f"unresolved jump-accounting flag at t = {tk:.4f}"
                break
        worst, sel = ev.audit_births(mesh, F, births, rng)
        if worst > ev.AUDIT_TOL:
            exr = ev.exact_births(mesh, F, births)
            flags.append(dict(kind="interpolation", t=tk, audit_max=worst, audited=sel.tolist(), recheck=exr))
            if not exr["ok"]:
                status, note = FAILED, f"unresolved interpolation-audit flag at t = {tk:.4f}"
                break
        for key in ("n_prop_births", "n_prop_deaths", "n_eff_births", "n_eff_deaths", "n_ephemeral", "jump_resid",
                    "I_before", "R", "B", "I_after"):
            ticks[key].append(acct[key])
        ticks["t"].append(tk); ticks["audit_max"].append(worst); ticks["N"].append(len(ids_n))
        mesh = ev.Mesh(pos_n, Lx, Ly)
        ids, F = ids_n, F_n
        stepper = ev.Ros2(mesh.g, h, P, dt)
        y = F.ravel()
    res = dict(status=status, note=note, flags=flags, t_reached=k_done * delta, n_factor=n_factor,
               wall_s=time.time() - t_wall, N_final=int(len(ids)),
               max_jump_resid=float(max(ticks["jump_resid"], default=0.0)),
               max_audit=float(max(ticks["audit_max"], default=0.0)))
    if status == COMPLETE and k_done < K:
        res["status"] = INCOMPLETE
    return res, checks, ticks



def save_checks(path, checks, ticks):
    """Committed per-run record: sparse check records (concatenated with offsets) and the per-tick table."""
    off = np.cumsum([0] + [len(c["rec"]["s"]) for c in checks])
    cat = lambda f, dt_: (np.concatenate([np.asarray(c["rec"][f], dt_) for c in checks]) if off[-1]
                          else np.zeros((0,) + np.asarray(checks[0]["rec"][f]).shape[1:], dt_))
    out = dict(t=np.array([c["t"] for c in checks]), off=off,
               rows=cat("rows", np.int64), pos=cat("pos", np.float64).reshape(-1, 2), ids=cat("ids", np.int64),
               s=cat("s", np.float64), N=np.array([c["rec"]["N"] for c in checks]),
               smin=np.array([c["rec"]["smin"] for c in checks]),
               smin_pos=np.array([c["rec"]["smin_pos"] for c in checks]).reshape(-1, 2),
               Lx=checks[0]["rec"]["Lx"], Ly=checks[0]["rec"]["Ly"])
    for k, v in ticks.items():
        out["tick_" + k] = np.asarray(v, float if k in ("t", "jump_resid", "audit_max", "I_before", "R", "B",
                                                        "I_after") else np.int64)
    np.savez_compressed(path, **out)


def load_checks(path):
    with np.load(path) as z:                                 # copy out and close (Windows keeps open files locked)
        d = {k: z[k] for k in z.files}
    recs = []
    for i in range(len(d["t"])):
        a, b = d["off"][i], d["off"][i + 1]
        recs.append(dict(t=float(d["t"][i]), rows=d["rows"][a:b], pos=d["pos"][a:b], ids=d["ids"][a:b],
                         s=d["s"][a:b], N=int(d["N"][i]), smin=float(d["smin"][i]), smin_pos=d["smin_pos"][i],
                         Lx=float(d["Lx"]), Ly=float(d["Ly"])))
    return recs, d


def summarize(checks, theta, Lx, Ly, t_end):
    """Class and mobility from the driver's full-field detections (the evaluator recomputes them from the sparse
    records)."""
    c0 = checks[0]["parts"]
    lists = [c["parts"] for c in checks[1:]]
    out = dict(n0=len(c0), counts=[len(c) for c in lists],
               depths=[[p["depth"] for p in c] for c in lists],
               centres=[[np.asarray(p["centre"]).tolist() for p in c] for c in lists])
    if len(c0) != 1:
        out["cls"] = dict(cls="unresolved", sub="start not single")
        return out
    out["c0"] = np.asarray(c0[0]["centre"]).tolist()
    r = ev.classify_run(lists, c0[0]["centre"], Lx, Ly)
    if r.get("cls") == "single":
        r["D_p"] = ev.msd_dp(r["track"])
    out["cls"] = r
    return out


# ------------------------------------------------------------------ registered runs
def run_specs():
    """The 24 registered runs (proposal §5). Part A: b = 28, t = 500, nu in RATES; Part B: b = 29, t = 300, nu 0, 0.2."""
    sp_ = []
    for part, b, t_end, rates in (("A", 28.0, 500.0, RATES), ("B", 29.0, 300.0, [0.0, 0.2])):
        for nu in rates:
            for slot, pl in DRUNS:
                sp_.append(dict(key=f"{part}_nu{nu:g}_D{slot}_p{pl}", part=part, b=b, t_end=t_end, nu=nu, slot=slot,
                                pl=pl, seed=dict(gate=ev.GATE, b_idx=0 if part == "A" else 1,
                                                 om_idx=RATES.index(nu), stage=1, replicate=2 * slot + pl)))
    return sp_


def start_state(slot, pl, b):
    """MG4c's stored float64 final state at b = 28 on panel graph D<slot> placement p<pl>, with the panel's h. At b = 29
    the potentials are shifted by S0(28) - S0(29) (only phi_Y changes: -0.9/14)."""
    import run_mg4c as R4
    dec = json.load(open(DEC4))
    g, h, c = R4.graph_for(R4.job(28.0, "D", slot, pl), dec)
    F = np.load(os.path.join(MG4C_FIELDS, f"b28.000000_D{slot}_p{pl}.npz"))["final"].astype(np.float64)
    shift = np.array(steady_state({**P17, "b": 28.0})) - np.array(steady_state({**P17, "b": b}))
    F = F + shift[:, None]
    return g, float(h), F, shift


def timing_params():
    d = json.load(open(P0DEC))
    if not d.get("proceed"):
        raise SystemExit("Part 0 did not authorise Part A (see results/MG5/part0/decisions.json)")
    return float(d["dt"]), float(d["delta_fast"])


def run_one(spec):
    path = os.path.join(OUT, "runs", spec["key"] + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    dt, delta_fast = timing_params()
    g, h, F0, shift = start_state(spec["slot"], spec["pl"], spec["b"])
    Lx, Ly = g["Lx"], g["Ly"]
    m0 = ev.Mesh(g["pos"], Lx, Ly)
    assert np.array_equal(m0.g["I"], g["I"]) and np.array_equal(m0.g["J"], g["J"]), "start mesh differs from MG4c graph"
    P = {**P17, "b": spec["b"]}
    theta = THETA_28 if spec["part"] == "A" else THETA_29
    nu = spec["nu"]
    delta = 1.0 if nu <= 0.04 else delta_fast
    seed = seed_for(**spec["seed"])
    rec = dict(key=spec["key"], spec=spec, b=spec["b"], nu=nu, h=h, dt=dt, delta=delta, theta=theta, Lx=Lx, Ly=Ly,
               t_end=spec["t_end"], phi_shift=shift.tolist(), platform=platform_info(), code_revision=code_revision(),
               registered="GATES.md MG5; spec/MG5_proposal.md draft 3 (e2a25e3)")
    events = None
    if nu > 0:
        args = (np.arange(N0), N0, nu, spec["t_end"], Lx, Ly, seed)
        events = ev.gen_events(*args)
        hsh = ev.events_hash(events)
        hsh2 = ev.regen_hash_fresh_process(args)
        rec["replay"] = dict(seed_for=spec["seed"], seed=seed, bit_generator="PCG64 (numpy.random.default_rng)",
                             numpy=np.__version__, code_revision=rec["code_revision"],
                             schema="t <f8, type i1 (0 death, 1 birth), id <i8, x <f8, y <f8; time order; "
                                    "SHA-256 over the five columns concatenated in that order",
                             n_events=int(len(events)), sha256=hsh, fresh_process_sha256=hsh2)
        if hsh != hsh2:
            rec.update(status=INCOMPLETE, note="event-list hash mismatch in a fresh process")
            json.dump(rec, open(path, "w"), indent=1)
            return rec
        os.makedirs(os.path.join(LOCAL, "events"), exist_ok=True)
        np.savez(os.path.join(LOCAL, "events", spec["key"] + ".npz"), events=events)
    os.makedirs(os.path.join(LOCAL, "fields", spec["key"]), exist_ok=True)

    def snap(t, pos, ids, F):
        np.savez_compressed(os.path.join(LOCAL, "fields", spec["key"], f"t{t:06.1f}.npz"),
                            pos=pos, ids=ids, F=F if is_multiple(t, spec["t_end"]) else F.astype(np.float32),
                            sX=ev.smoothed_X(ev.Mesh(pos, Lx, Ly).g, F).astype(np.float32))

    res, checks, ticks = run_ticks(g["pos"], F0, h, P, theta, nu, events, delta, dt, spec["t_end"], Lx, Ly,
                                   seed_audit=seed, snap=snap)
    os.makedirs(os.path.join(OUT, "checks"), exist_ok=True)
    save_checks(os.path.join(OUT, "checks", spec["key"] + ".npz"), checks, ticks)
    rec.update(res)
    if res["status"] == COMPLETE:
        rec["summary"] = summarize(checks, theta, Lx, Ly, spec["t_end"])
    rec["wall_s"] = res["wall_s"]
    json.dump(rec, open(path, "w"), indent=1, default=float)
    return rec


def ledger(key, wall):
    path = os.path.join(OUT, "budget_ledger.csv")
    new = not os.path.exists(path)
    pi = platform_info()
    with open(path, "a", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        if new:
            wr.writerow(["key", "wall_s", "workers", "machine", "os", "python", "numpy", "utc_time"])
        wr.writerow([key, f"{wall:.2f}", WORKERS, pi["machine"], pi["os"], pi["python"], pi["numpy"],
                     time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())])


def spent():
    path = os.path.join(OUT, "budget_ledger.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8"))) if os.path.exists(path) else []
    return sum(float(r["wall_s"]) for r in rows) / 3600, {r["key"] for r in rows}


def write_manifest():
    rows = []
    for dp, _, fs in os.walk(LOCAL):
        for f in sorted(fs):
            p = os.path.join(dp, f)
            rows.append([os.path.relpath(p, OUT).replace(os.sep, "/"), os.path.getsize(p), sha256(p)])
    with open(os.path.join(OUT, "local_manifest.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["path", "bytes", "sha256"])
        wr.writerows(sorted(rows))


def _worker(spec):
    return run_one(spec)


def main():
    os.makedirs(os.path.join(OUT, "runs"), exist_ok=True)
    timing_params()
    go.continuum2d()
    specs = run_specs()
    todo = [s for s in specs if not os.path.exists(os.path.join(OUT, "runs", s["key"] + ".json"))]
    # the fast runs first, so the longest jobs do not finish last
    todo.sort(key=lambda s: -s["nu"])
    stopped = None
    with ProcessPoolExecutor(WORKERS) as ex:
        live = {}
        while todo or live:
            while todo and len(live) < WORKERS:
                ch, done = spent()
                n = len(done | set(live.values()))
                if ch >= CUTOFF or n >= RUN_CAP:
                    stopped = f"ledger {ch:.2f} core-h, {n} runs started"
                    todo.clear()
                    break
                s = todo.pop(0)
                live[ex.submit(_worker, s)] = s["key"]
            if not live:
                break
            fut = next(as_completed(live))
            key = live.pop(fut)
            r = fut.result()
            ledger(key, r.get("wall_s", 0.0))
            cls = r.get("summary", {}).get("cls", {})
            print(f"{key}: {r['status']}, wall {r.get('wall_s', 0):.0f}s, "
                  f"class {cls.get('cls')} {cls.get('mobility', '')}", flush=True)
    write_manifest()
    ch, done = spent()
    log = dict(status=f"stopped: {stopped}" if stopped else "complete", core_h=ch, runs=len(done),
               platform=platform_info(), code_revision=code_revision())
    json.dump(log, open(os.path.join(OUT, "driver_log.json"), "w"), indent=1)
    print("MG5 driver finished:", log["status"], f"{ch:.2f} core-h, {len(done)} runs", flush=True)
    print("Evaluate with: python src/evaluate_mg5.py", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "part0":
        import run_mg5_part0
        run_mg5_part0.main()
    else:
        main()
