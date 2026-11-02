"""MG6 driver (spec/MG6_proposal.md draft 3, frozen at 106b3e2; GATES.md MG6 and its pre-run implementation
specification). G-coupled birth-death rewriting of the Delaunay graph around a settled b = 28 particle.

    $env:MODELG_WORKERS=5; python src/run_mg6.py part0       (PowerShell, PI's PC): Part 0 checks (unregistered graph)
    $env:MODELG_WORKERS=5; python src/run_mg6.py             the 20 registered runs (needs results/MG6/part0/decisions.json)

Resumable at run level. Committed outputs (results/MG6/): runs/<key>.json, checks/<key>.npz, accepted/<key>.npz,
budget_ledger.csv, local_manifest.csv, driver_log.json. PC-only (results/MG6/local/, gitignored): candidate lists,
full-field snapshots and final states, listed with SHA-256 in local_manifest.csv.
"""
import os, sys
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):   # one BLAS thread per worker (MG017)
    os.environ.setdefault(_v, "1")
import json, time, csv
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(__file__))
import graph_ode as go
import evolve as ev
import mg6
import quad6_check as qc
import run_mg5 as R5
from lattice import P17
from runner import seed_for

ROOT = os.path.join(os.path.dirname(__file__), "..", "results")
OUT = os.environ.get("MG6_OUT") or os.path.join(ROOT, "MG6")
LOCAL = os.path.join(OUT, "local")
P0DEC = os.path.join(OUT, "part0", "decisions.json")
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))

B, NU, T_END, CHECK_EVERY, SNAP_EVERY = 28.0, 0.2, 200.0, 10.0, 50.0
THETA = R5.THETA_28
N0 = R5.N0
RUN_CAP, CUTOFF = 24, 40.0
DRUNS = R5.DRUNS
ARMS = [("DEC", 0.0, False), ("SH2", 2.0, True), ("SH4", 4.0, True), ("CP2", 2.0, False), ("CP4", 4.0, False)]
QTOL_REL, QTOL_ABS = 0.01, 0.5               # in-run quadrature cross-check (pre-run spec): 1% relative, plus an
                                             # absolute floor of 0.5 nodes for P0 and Pc only
FAILED, INVALID_BOUND, COMPLETE, INCOMPLETE = "FAILED-NUMERICAL", "INVALID-BOUND", "COMPLETE", "INCOMPLETE"


def cand_seed(slot, pl):
    return seed_for(mg6.GATE, 0, 0, 1, 2 * slot + pl)


def shift_seed(slot, pl):
    return seed_for(mg6.GATE, 0, 1, 1, 2 * slot + pl)


def qclose(a, b, floor=0.0):
    return abs(a - b) <= QTOL_REL * max(abs(a), abs(b)) + floor


# ------------------------------------------------------------------ the engine (registered runs and Part 0)
def run_ticks6(pos0, F0, h, P, theta, coupled, delta, dt, t_end, Lx, Ly, kappa, shift, grid, rho0, seed_audit=0,
               snap=None, log=None, measure_every=CHECK_EVERY, cross=True, fixed_field=None):
    """Proposal §3: for k = 1..K integrate on M_(k-1) to t_k (skipped if fixed_field: Part 0 graph-process checks);
    full-field bound check (max nodal phi_G <= 0.4/kappa) before thinning; scheduled check on M_(k-1) and the boundary
    field (detection, sparse record, §4.2 measurements, quadrature cross-check); tick k's edits (deaths, then accepted
    births interpolated on M_(k-1)); rebuild M_k. Returns (result, checks, ticks)."""
    t_wall = time.time()
    mesh = ev.Mesh(pos0, Lx, Ly)
    ids = np.arange(len(pos0), dtype=np.int64)
    F = None if fixed_field is not None else np.asarray(F0, float).copy()
    stepper = None if fixed_field is not None else ev.Ros2(mesh.g, h, P, dt)
    K = ev.n_steps(t_end, delta)
    nsub = ev.n_steps(delta, dt)
    rng = np.random.default_rng(seed_audit)
    lim = mg6.phi_lim(kappa)
    checks, flags = [], []
    ticks = {k: [] for k in ("t", "n_cand", "n_eph", "n_eval", "n_acc", "n_dead", "n_eff_births", "n_eff_deaths",
                             "jump_resid", "audit_max", "N", "max_phiG")}
    status, note = COMPLETE, None

    def do_check(tk, F):
        rec = dict(t=tk)
        if fixed_field is None:
            s = ev.smoothed_X(mesh.g, F)
            rec["rec"] = ev.sparse_record(mesh.pos, ids, s, Lx, Ly)
            parts = ev.detect_full(mesh.pos, s, theta, Lx, Ly)
            rec["parts"] = parts
            c = parts[0]["centre"] if parts else None
        else:
            rec["parts"] = []
            c = fixed_field["centre"]
        if c is not None:
            m = mg6.measure(mesh, F, c, kappa, shift, grid, rho0)
            rec["m"] = m
            if cross:
                x = qc.cross_check(mesh.pos, F[0], Lx, Ly, c, kappa, shift, rho0)
                rec["x"] = x
                ok = (qclose(m["P0"][1], x["P0"], QTOL_ABS) and qclose(m["Pc"][1], x["Pc"], QTOL_ABS)
                      and qclose(m["M_minus"], x["M_minus"]) and qclose(m["M_plus"], x["M_plus"]))
                rec["x_ok"] = bool(ok)
        return rec

    def field_of(pos):
        return np.vstack([fixed_field["G"](pos), np.zeros(len(pos)), np.zeros(len(pos))])

    if fixed_field is not None:
        F = field_of(mesh.pos)
    checks.append(do_check(0.0, F))
    y = F.ravel()
    k_done = 0
    for k in range(1, K + 1):
        tk = k * delta
        N = len(ids)
        if fixed_field is None:
            stepper.factor(y)
            y = stepper.advance(y, nsub)
            if not np.all(np.isfinite(y)) or np.abs(y).max() > 1e3:
                status, note = FAILED, f"non-finite or |phi| > 1e3 at t = {tk:.4f}"
                break
            bal = go.balances(y, stepper.rhs(0.0, y), N, P)
            if not bal["ok"]:
                hp = ev.hp_balance(mesh.g, h, P, y)
                flags.append(dict(kind="balance", t=tk, r1=bal["r1"], r2=bal["r2"], recheck=hp))
                if not hp["resolved"]:
                    status, note = FAILED, f"unresolved balance flag at t = {tk:.4f}"
                    break
            F = y.reshape(3, N)
        mx = float(F[0].max())
        if mx > lim:                                          # §3.2: full-field bound, before thinning
            status, note = INVALID_BOUND, f"max nodal phi_G {mx:.6f} > {lim:.4f} at t = {tk:.4f}"
            flags.append(dict(kind="bound", t=tk, max_phiG=mx, node=int(np.argmax(F[0])),
                              pos=mesh.pos[int(np.argmax(F[0]))].tolist()))
            break
        if R5.is_multiple(tk, measure_every):
            checks.append(do_check(tk, F))
            if checks[-1].get("x_ok") is False:
                flags.append(dict(kind="quadrature", t=tk, m=[float(checks[-1]["m"]["P0"][1]),
                                                              float(checks[-1]["m"]["Pc"][1]),
                                                              checks[-1]["m"]["M_minus"], checks[-1]["m"]["M_plus"]],
                                  x=checks[-1]["x"]))
                re = do_check(tk, F)                          # re-evaluation (pre-run spec)
                if re.get("x_ok") is False:
                    status, note = FAILED, f"quadrature cross-check failed twice at t = {tk:.4f}"
                    break
            if log:
                log(f"t = {tk:g}: {len(checks[-1]['parts'])} particle(s), N = {N}")
        if snap is not None and (R5.is_multiple(tk, SNAP_EVERY) or k == K):
            snap(tk, mesh.pos, ids, F)
        k_done = k
        if k == K:                                            # no edits after the final boundary
            continue
        E, cnt = coupled.tick_edits(mesh, F, ids, tk - delta, tk)
        if coupled.max_p > 1.0 + 1e-15:                       # candidate-level assertion (§3.2)
            status, note = INVALID_BOUND, f"acceptance probability {coupled.max_p} > 1 at t = {tk:.4f}"
            break
        pos_n, ids_n, F_n, acct, births = ev.apply_tick(mesh, F, ids, E)
        if not acct["jump_ok"]:
            dead_rows = np.flatnonzero(~np.isin(ids, ids_n))
            exr = ev.exact_jump(F, dead_rows, births["val"], F_n)
            flags.append(dict(kind="jump", t=tk, resid=acct["jump_resid"], recheck=exr))
            if not exr["resolved"]:
                status, note = FAILED, f"unresolved jump-accounting flag at t = {tk:.4f}"
                break
        worst, sel = ev.audit_births(mesh, F, births, rng)
        if worst > ev.AUDIT_TOL:
            exr = ev.exact_births(mesh, F, births)
            flags.append(dict(kind="interpolation", t=tk, audit_max=worst, audited=sel.tolist(), recheck=exr))
            if not exr["ok"]:
                status, note = FAILED, f"unresolved interpolation-audit flag at t = {tk:.4f}"
                break
        for key in ("n_cand", "n_eph", "n_eval", "n_acc", "n_dead"):
            ticks[key].append(cnt[key])
        ticks["n_eff_births"].append(acct["n_eff_births"]); ticks["n_eff_deaths"].append(acct["n_eff_deaths"])
        ticks["jump_resid"].append(acct["jump_resid"]); ticks["audit_max"].append(worst)
        ticks["t"].append(tk); ticks["N"].append(len(ids_n)); ticks["max_phiG"].append(mx)
        mesh = ev.Mesh(pos_n, Lx, Ly)
        ids = ids_n
        F = field_of(mesh.pos) if fixed_field is not None else F_n
        if fixed_field is None:
            stepper = ev.Ros2(mesh.g, h, P, dt)
        y = F.ravel()
    res = dict(status=status, note=note, flags=flags, t_reached=k_done * delta, wall_s=time.time() - t_wall,
               N_final=int(len(ids)), max_jump_resid=float(max(ticks["jump_resid"], default=0.0)),
               max_audit=float(max(ticks["audit_max"], default=0.0)),
               max_phiG=float(max(ticks["max_phiG"], default=-np.inf)), max_p=coupled.max_p)
    if status == COMPLETE and k_done < K:
        res["status"] = INCOMPLETE
    return res, checks, ticks


M_ARRAYS = ("N_R", "D0", "Dc", "P0", "Pc", "quad_area", "ann_count", "ann_dn0", "ann_dnc", "ann_P0", "ann_Pc",
            "ann_phiG", "fine_count", "fine_dn0", "fine_dnc", "fine_P0", "fine_Pc", "fine_phiG")
M_SCALARS = ("N", "rho_bar", "rho0", "M_minus", "M_plus", "M")
X_SCALARS = ("P0", "Pc", "M_minus", "M_plus")


def save_checks6(path, checks, ticks, accepted=None):
    """Committed per-run record: MG5's sparse records (when the run has fields), the §4.2 measurements per check (NaN
    where no particle), the cross-check values, the per-tick table, and the accepted-candidate bitmap."""
    T = len(checks)
    out = dict(t=np.array([c["t"] for c in checks]))
    if all("rec" in c for c in checks):
        off = np.cumsum([0] + [len(c["rec"]["s"]) for c in checks])
        cat = lambda f, dt_, shp=(): (np.concatenate([np.asarray(c["rec"][f], dt_).reshape((-1,) + shp)
                                                      for c in checks]))
        out.update(off=off, rows=cat("rows", np.int64), pos=cat("pos", np.float64, (2,)), ids=cat("ids", np.int64),
                   s=cat("s", np.float64), smin=np.array([c["rec"]["smin"] for c in checks]),
                   smin_pos=np.array([c["rec"]["smin_pos"] for c in checks]).reshape(-1, 2),
                   Lx=checks[0]["rec"]["Lx"], Ly=checks[0]["rec"]["Ly"],
                   n_parts=np.array([len(c["parts"]) for c in checks]),
                   centre=np.array([c["parts"][0]["centre"] if c["parts"] else (np.nan, np.nan) for c in checks]))
    ref = next((c["m"] for c in checks if "m" in c), None)
    if ref is not None:
        for a in M_ARRAYS:
            out["m_" + a] = np.array([c["m"][a] if "m" in c else np.full(len(ref[a]), np.nan) for c in checks], float)
        for a in M_SCALARS:
            out["m_" + a] = np.array([c["m"][a] if "m" in c else np.nan for c in checks], float)
        for a in X_SCALARS:
            out["x_" + a] = np.array([c["x"][a] if "x" in c else np.nan for c in checks], float)
    for k, v in ticks.items():
        out["tick_" + k] = np.asarray(v, float if k in ("t", "jump_resid", "audit_max", "max_phiG") else np.int64)
    if accepted is not None:
        out["accepted_bits"] = np.packbits(accepted)
        out["n_candidates"] = len(accepted)
    np.savez_compressed(path, **out)


def summarize6(checks, Lx, Ly):
    c0 = checks[0]["parts"]
    lists = [c["parts"] for c in checks[1:]]
    out = dict(n0=len(c0), counts=[len(c) for c in lists], depths=[[p["depth"] for p in c] for c in lists])
    if len(c0) != 1:
        out["cls"] = dict(cls="unresolved", sub="start not single")
        return out
    r = ev.classify_run(lists, c0[0]["centre"], Lx, Ly)
    if r.get("cls") == "single":
        r["D_p"] = ev.msd_dp(r["track"])
    out["cls"] = r
    return out


# ------------------------------------------------------------------ registered runs
def run_specs():
    """The 20 registered runs (proposal §3.3), in launch order: DEC, SH2, SH4, CP2, CP4 (nulls first, §6)."""
    return [dict(key=f"{arm}_D{slot}_p{pl}", arm=arm, kappa=kap, shifted=sh, slot=slot, pl=pl)
            for arm, kap, sh in ARMS for slot, pl in DRUNS]


def timing_params():
    d = json.load(open(P0DEC))
    if not d.get("proceed"):
        raise SystemExit("Part 0 did not authorise the registered runs (decisions.json: proceed is not true)")
    return d


def run_one(spec):
    path = os.path.join(OUT, "runs", spec["key"] + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    dec = timing_params()
    dt, delta, t_end = dec["dt"], dec["delta"], dec.get("t_end", T_END)
    g, h, F0, _ = R5.start_state(spec["slot"], spec["pl"], B)
    Lx, Ly = g["Lx"], g["Ly"]
    m0 = ev.Mesh(g["pos"], Lx, Ly)
    assert np.array_equal(m0.g["I"], g["I"]) and np.array_equal(m0.g["J"], g["J"]), "start mesh differs from MG4c graph"
    P = {**P17, "b": B}
    grid = mg6.Grid(Lx, Ly, dec.get("nx", mg6.NX), dec.get("ny", mg6.NY))
    rho0 = N0 / grid.A
    seed = cand_seed(spec["slot"], spec["pl"])
    args = (N0, NU, t_end, Lx, Ly, seed)
    cand = mg6.gen_candidates(*args)
    hsh, hsh2 = mg6.candidates_hash(cand), mg6.regen_hash_fresh_process(args)
    shift = mg6.draw_shift(shift_seed(spec["slot"], spec["pl"]), Lx, Ly) if spec["shifted"] else None
    rec = dict(key=spec["key"], spec=spec, b=B, nu=NU, kappa=spec["kappa"], h=h, dt=dt, delta=delta, t_end=t_end,
               theta=THETA, Lx=Lx, Ly=Ly, rho0=rho0, F=mg6.FBOUND, phi_lim=mg6.phi_lim(spec["kappa"]),
               shift=None if shift is None else shift.tolist(), grid=[grid.nx, grid.ny],
               platform=R5.platform_info(), code_revision=R5.code_revision(),
               registered="GATES.md MG6; spec/MG6_proposal.md draft 3 (106b3e2)",
               replay=dict(seed_for=[mg6.GATE, 0, 0, 1, 2 * spec["slot"] + spec["pl"]], seed=seed,
                           bit_generator="PCG64 (numpy.random.default_rng)", numpy=np.__version__,
                           schema="N0 <i8, then life0, tb, x, y, u, ell as <f8 arrays; SHA-256 over that byte stream",
                           n_candidates=int(len(cand["tb"])), sha256=hsh, fresh_process_sha256=hsh2))
    if hsh != hsh2:
        rec.update(status=INCOMPLETE, note="candidate-list hash mismatch in a fresh process", wall_s=0.0)
        json.dump(rec, open(path, "w"), indent=1, default=float)
        return rec
    os.makedirs(os.path.join(LOCAL, "candidates"), exist_ok=True)
    np.savez(os.path.join(LOCAL, "candidates", spec["key"] + ".npz"), **{k: cand[k] for k in mg6.CAND_FIELDS})
    os.makedirs(os.path.join(LOCAL, "fields", spec["key"]), exist_ok=True)

    def snap(t, pos, ids, F):
        np.savez_compressed(os.path.join(LOCAL, "fields", spec["key"], f"t{t:06.1f}.npz"),
                            pos=pos, ids=ids, F=F if R5.is_multiple(t, t_end) else F.astype(np.float32))

    coupled = mg6.Coupled(cand, spec["kappa"], shift)
    res, checks, ticks = run_ticks6(g["pos"], F0, h, P, THETA, coupled, delta, dt, t_end, Lx, Ly, spec["kappa"], shift,
                                    grid, rho0, seed_audit=seed, snap=snap)
    for d in ("checks",):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    save_checks6(os.path.join(OUT, "checks", spec["key"] + ".npz"), checks, ticks, coupled.accepted)
    rec.update(res)
    rec["n_accepted"] = int(coupled.accepted.sum())
    if res["status"] == COMPLETE:
        rec["summary"] = summarize6(checks, Lx, Ly)
    json.dump(rec, open(path, "w"), indent=1, default=float)
    return rec


def ledger(key, wall):
    path = os.path.join(OUT, "budget_ledger.csv")
    new = not os.path.exists(path)
    pi = R5.platform_info()
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
            rows.append([os.path.relpath(p, OUT).replace(os.sep, "/"), os.path.getsize(p), R5.sha256(p)])
    with open(os.path.join(OUT, "local_manifest.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["path", "bytes", "sha256"])
        wr.writerows(sorted(rows))


def main():
    os.makedirs(os.path.join(OUT, "runs"), exist_ok=True)
    dec = timing_params()
    cutoff = dec.get("cutoff_core_h", CUTOFF)
    todo = [s for s in run_specs() if not os.path.exists(os.path.join(OUT, "runs", s["key"] + ".json"))]
    stopped = None
    with ProcessPoolExecutor(WORKERS) as ex:
        live = {}
        while todo or live:
            while todo and len(live) < WORKERS:
                ch, done = spent()
                n = len(done | set(live.values()))
                if ch >= cutoff or n >= RUN_CAP:
                    stopped = f"ledger {ch:.2f} core-h, {n} runs started"
                    todo.clear()
                    break
                s = todo.pop(0)
                live[ex.submit(run_one, s)] = s["key"]
            if not live:
                break
            fut = next(as_completed(live))
            key = live.pop(fut)
            r = fut.result()
            ledger(key, r.get("wall_s", 0.0))
            cls = r.get("summary", {}).get("cls", {})
            print(f"{key}: {r['status']}, wall {r.get('wall_s', 0):.0f}s, class {cls.get('cls')} "
                  f"{cls.get('mobility', '')}", flush=True)
    write_manifest()
    ch, done = spent()
    log = dict(status=f"stopped: {stopped}" if stopped else "complete", core_h=ch, runs=len(done),
               platform=R5.platform_info(), code_revision=R5.code_revision())
    json.dump(log, open(os.path.join(OUT, "driver_log.json"), "w"), indent=1)
    print("MG6 driver finished:", log["status"], f"{ch:.2f} core-h, {len(done)} runs", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "part0":
        import run_mg6_part0
        run_mg6_part0.main()
    else:
        main()
