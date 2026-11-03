"""MG7 driver (spec/MG7_proposal.md draft 2, frozen at 5b44cde; GATES.md MG7 and its pre-run implementation
specification). Mass scaling of the coupled density deficit: CP2 and SH2 at b = 27.5 and 28.5 from MG4c's own settled
states at those b. Everything not stated is MG6 as registered and amended (MG021), with MG6's code reused unchanged
(src/mg6.py and src/run_mg6.py are byte-identical to the MG6 production revision c795770; tests/test_mg7.py checks it).

    $env:MODELG_WORKERS=8; python src/run_mg7.py part0     (PowerShell, PI's PC): Part 0 (unregistered inputs)
    $env:MODELG_WORKERS=8; python src/run_mg7.py           the 16 registered runs (needs results/MG7/part0/decisions.json)

Resumable at run level. Committed outputs (results/MG7/): runs/<key>.json, checks/<key>.npz, budget_ledger.csv,
local_manifest.csv, driver_log.json, endpoint_manifest.csv (pre-run). PC-only (results/MG7/local/, gitignored):
candidate lists and full-field snapshots, listed with SHA-256 in local_manifest.csv.
"""
import os, sys
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):   # one BLAS thread per worker (MG017)
    os.environ.setdefault(_v, "1")
import json, time, csv, shutil
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(__file__))
import evolve as ev
import mg6
import run_mg5 as R5
import run_mg6 as R6
from lattice import P17
from runner import seed_for

GATE = 7
ROOT = os.path.join(os.path.dirname(__file__), "..", "results")
OUT = os.environ.get("MG7_OUT") or os.path.join(ROOT, "MG7")
LOCAL = os.path.join(OUT, "local")
P0DEC = os.path.join(OUT, "part0", "decisions.json")
EP_MANIFEST = os.path.join(ROOT, "MG7", "endpoint_manifest.csv")     # committed before any MG7 run
WORKERS = int(os.environ.get("MODELG_WORKERS", "2"))

B_END = (27.5, 28.5)
NU, T_END, DELTA, DT = R6.NU, R6.T_END, 0.1, 0.05                     # MG6 as amended (MG021)
NX, NY = mg6.NX, mg6.NY                                               # MG6 production grid (398 x 402)
THETA = R5.THETA_28                                                   # fixed detector threshold, both endpoints
THETA_B = {27.5: -3.857025988052371, 28.5: -4.840265200775234}       # MG4c b-specific thresholds (descriptive)
N0 = R5.N0
DRUNS = R5.DRUNS
KAPPA = 2.0
ARMS = [("SH2", True), ("CP2", False)]                                # launch order: nulls first
RUN_CAP, CUTOFF = 22, 90.0
MARGIN = 0.14                                                         # 0.7 * phi_lim(2)
COMPLETE, INCOMPLETE, FAILED, INVALID_BOUND = R6.COMPLETE, R6.INCOMPLETE, R6.FAILED, R6.INVALID_BOUND
MG6_REVISION = "c795770"


def bkey(b):
    return f"b{b:.2f}"


def run_key(arm, b, slot, pl):
    return f"{arm}_{bkey(b)}_D{slot}_p{pl}"


def cand_seed(slot, pl):
    """One candidate stream per start, shared across both endpoints and both arms (proposal §3)."""
    return seed_for(GATE, 0, 0, 1, 2 * slot + pl)


def shift_seed(slot, pl):
    return seed_for(GATE, 0, 1, 1, 2 * slot + pl)


# ------------------------------------------------------------------ start fields (review item 1)
def field_relpath(b, slot, pl):
    return f"MG4c/fields/b{b:.6f}_D{slot}_p{pl}.npz"


def read_manifest(path=EP_MANIFEST):
    with open(path, encoding="utf-8") as f:
        return {r["path"]: r for r in csv.DictReader(f)}


def write_endpoint_manifest(path=EP_MANIFEST):
    """Run once before registration-time commit: SHA-256 of the eight MG4c endpoint files."""
    rows = []
    for b in B_END:
        for slot, pl in DRUNS:
            rp = field_relpath(b, slot, pl)
            p = os.path.join(ROOT, rp)
            rows.append([rp, os.path.getsize(p), R5.sha256(p)])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["path", "bytes", "sha256"])
        wr.writerows(rows)


def verify_endpoint_files(manifest=None):
    """(ok, rows): every endpoint file present, listed and matching its SHA-256."""
    man = read_manifest() if manifest is None else manifest
    rows, ok = [], True
    for b in B_END:
        for slot, pl in DRUNS:
            rp = field_relpath(b, slot, pl)
            p = os.path.join(ROOT, rp)
            sha = R5.sha256(p) if os.path.exists(p) else None
            good = rp in man and sha is not None and sha == man[rp]["sha256"]
            ok &= good
            rows.append(dict(path=rp, sha256=sha, listed=rp in man, ok=bool(good)))
    return bool(ok), rows


def start_state7(slot, pl, b, manifest=None):
    """MG4c's own stored float64 final state at the requested b on panel graph D<slot> placement p<pl>, with the panel's
    h, and NO steady-state shift (MG6's run_mg5.start_state always loads b = 28 and shifts it; it is not used here).
    The file must match its manifest SHA-256. Returns (graph, h, F, sha256)."""
    import run_mg4c as R4
    if b not in B_END:
        raise ValueError(f"b = {b} is not an MG7 endpoint")
    man = read_manifest() if manifest is None else manifest
    rp = field_relpath(b, slot, pl)
    p = os.path.join(ROOT, rp)
    sha = R5.sha256(p)
    if rp not in man or man[rp]["sha256"] != sha:
        raise RuntimeError(f"start field {rp}: SHA-256 {sha} does not match the endpoint manifest")
    dec = json.load(open(R5.DEC4))
    g, h, _ = R4.graph_for(R4.job(b, "D", slot, pl), dec)
    with np.load(p) as z:
        raw = z["final"]
    if raw.dtype != np.float64 or raw.shape != (3, g["N"]):
        raise RuntimeError(f"start field {rp}: dtype {raw.dtype}, shape {raw.shape}")
    F = np.array(raw, dtype=np.float64, copy=True)
    with np.load(p) as z:
        if not np.array_equal(F, z["final"]):
            raise RuntimeError(f"start field {rp}: loaded array differs from the file")
    return g, float(h), F, sha


def params_for(b):
    P = {**P17, "b": float(b)}
    assert P["b"] == float(b) and b in B_END, "reaction parameter b is not the requested endpoint"
    return P


# ------------------------------------------------------------------ registered runs
def run_specs():
    """The 16 registered runs in launch order: SH2 (27.5, 28.5), then CP2 (27.5, 28.5); four starts each."""
    return [dict(key=run_key(arm, b, slot, pl), arm=arm, b=b, kappa=KAPPA, shifted=sh, slot=slot, pl=pl)
            for arm, sh in ARMS for b in B_END for slot, pl in DRUNS]


def part0_decision(path=P0DEC):
    if not os.path.exists(path):
        raise SystemExit("MG7 Part 0 has not run (results/MG7/part0/decisions.json missing)")
    d = json.load(open(path))
    bpath = os.path.join(os.path.dirname(path), "budget_decision.json")   # a recorded PI decision, if one is needed
    if not d.get("proceed") and os.path.exists(bpath):
        bd = json.load(open(bpath))
        stop = d.get("stop") or ""
        if (stop.startswith("projected cost") and bd.get("pi_approved")
                and d["timing"]["projected_core_h_16_runs"] <= bd["cutoff_core_h"]):
            d = dict(d, proceed=True, cutoff_core_h=bd["cutoff_core_h"], budget_decision=bd, stop=None)
    if not d.get("proceed"):
        raise SystemExit(f"Part 0 did not authorise the registered runs: {d.get('stop')}")
    if (d.get("delta"), d.get("dt")) != (DELTA, DT):
        raise SystemExit("Part 0 decisions do not carry the registered Delta 0.1, dt 0.05")
    return d


def _set_aside(path):
    """MG021 lesson: never overwrite an output silently. An existing file is moved to superseded/ with its SHA-256 in the
    name and the move is logged."""
    if not os.path.exists(path):
        return None
    d = os.path.join(OUT, "superseded")
    os.makedirs(d, exist_ok=True)
    sha = R5.sha256(path)
    dst = os.path.join(d, f"{os.path.basename(path)}.{sha[:16]}")
    shutil.move(path, dst)
    with open(os.path.join(d, "moves.csv"), "a", encoding="utf-8") as f:
        f.write(f"{os.path.relpath(path, OUT)},{sha},{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
    return dst


def run_one(spec):
    path = os.path.join(OUT, "runs", spec["key"] + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    dec = part0_decision()
    b = spec["b"]
    g, h, F0, sha = start_state7(spec["slot"], spec["pl"], b)
    Lx, Ly = g["Lx"], g["Ly"]
    m0 = ev.Mesh(g["pos"], Lx, Ly)
    assert np.array_equal(m0.g["I"], g["I"]) and np.array_equal(m0.g["J"], g["J"]), "start mesh differs from MG4c graph"
    P = params_for(b)
    grid = mg6.Grid(Lx, Ly, NX, NY)
    rho0 = N0 / grid.A
    seed = cand_seed(spec["slot"], spec["pl"])
    args = (N0, NU, T_END, Lx, Ly, seed)
    cand = mg6.gen_candidates(*args)
    hsh, hsh2 = mg6.candidates_hash(cand), mg6.regen_hash_fresh_process(args)
    shift = mg6.draw_shift(shift_seed(spec["slot"], spec["pl"]), Lx, Ly) if spec["shifted"] else None
    rec = dict(key=spec["key"], spec=spec, b=b, P_b=P["b"], nu=NU, kappa=KAPPA, h=h, dt=DT, delta=DELTA, t_end=T_END,
               theta=THETA, Lx=Lx, Ly=Ly, rho0=rho0, F=mg6.FBOUND, phi_lim=mg6.phi_lim(KAPPA),
               shift=None if shift is None else shift.tolist(), grid=[grid.nx, grid.ny],
               start_field=dict(path=field_relpath(b, spec["slot"], spec["pl"]), sha256=sha, shift_applied=False),
               platform=R5.platform_info(), code_revision=R5.code_revision(),
               registered="GATES.md MG7; spec/MG7_proposal.md draft 2 (5b44cde)",
               replay=dict(seed_for=[GATE, 0, 0, 1, 2 * spec["slot"] + spec["pl"]], seed=seed,
                           bit_generator="PCG64 (numpy.random.default_rng)", numpy=np.__version__,
                           schema="N0 <i8, then life0, tb, x, y, u, ell as <f8 arrays; SHA-256 over that byte stream",
                           n_candidates=int(len(cand["tb"])), sha256=hsh, fresh_process_sha256=hsh2))
    os.makedirs(os.path.join(OUT, "runs"), exist_ok=True)
    if hsh != hsh2:
        rec.update(status=INCOMPLETE, note="candidate-list hash mismatch in a fresh process", wall_s=0.0)
        json.dump(rec, open(path, "w"), indent=1, default=float)
        return rec
    os.makedirs(os.path.join(LOCAL, "candidates"), exist_ok=True)
    cpath = os.path.join(LOCAL, "candidates", spec["key"] + ".npz")
    _set_aside(cpath)
    np.savez(cpath, **{k: cand[k] for k in mg6.CAND_FIELDS})
    fdir = os.path.join(LOCAL, "fields", spec["key"])
    if os.path.isdir(fdir):
        for f in sorted(os.listdir(fdir)):
            _set_aside(os.path.join(fdir, f))
    os.makedirs(fdir, exist_ok=True)

    def snap(t, pos, ids, F):
        np.savez_compressed(os.path.join(fdir, f"t{t:06.1f}.npz"),
                            pos=pos, ids=ids, F=F if R5.is_multiple(t, T_END) else F.astype(np.float32))

    coupled = mg6.Coupled(cand, KAPPA, shift)
    res, checks, ticks = R6.run_ticks6(g["pos"], F0, h, P, THETA, coupled, DELTA, DT, T_END, Lx, Ly, KAPPA, shift,
                                       grid, rho0, seed_audit=seed, snap=snap)
    os.makedirs(os.path.join(OUT, "checks"), exist_ok=True)
    kpath = os.path.join(OUT, "checks", spec["key"] + ".npz")
    _set_aside(kpath)
    R6.save_checks6(kpath, checks, ticks, coupled.accepted)
    rec.update(res)
    rec["n_accepted"] = int(coupled.accepted.sum())
    rec["bound_margin_ok"] = bool(res["max_phiG"] <= MARGIN)
    if res["status"] == COMPLETE:
        rec["summary"] = R6.summarize6(checks, Lx, Ly)
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
    return sum(float(r["wall_s"]) for r in rows) / 3600, [r["key"] for r in rows]


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
    ok, rows = verify_endpoint_files()
    if not ok:
        raise SystemExit(f"endpoint start files do not match the manifest: {[r for r in rows if not r['ok']]}")
    dec = part0_decision()
    cutoff = dec.get("cutoff_core_h", CUTOFF)
    os.makedirs(os.path.join(OUT, "runs"), exist_ok=True)
    todo = [s for s in run_specs() if not os.path.exists(os.path.join(OUT, "runs", s["key"] + ".json"))]
    stopped = None
    with ProcessPoolExecutor(WORKERS) as ex:
        live = {}
        while todo or live:
            while todo and len(live) < WORKERS:
                ch, done = spent()
                n = len(done) + len(live)                     # every started run counts toward the cap
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
    print("MG7 driver finished:", log["status"], f"{ch:.2f} core-h, {len(done)} runs", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "part0":
        import run_mg7_part0
        run_mg7_part0.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "--write-endpoint-manifest":
        write_endpoint_manifest()
    else:
        main()
