"""MG8 production driver (spec/MG8_proposal.md draft 4, registered at 3ab4ec8; §3.2, §6, §7).

Usage (PI's PC):  $env:MODELG_WORKERS=8; python src/run_mg8.py
Resumable at run level: a run whose record exists is skipped. Requires results/MG8/part0/decisions.json with
proceed = true (Part 0, src/part0_mg8.py). Production cutoff 20 core-h (ledger of per-run wall time), cap 64 runs.
Order: B1 k = 0, then the rest of the ladder, then B2. Exit codes: 0 done, 2 budget stop (BUDGET-SHORT), 5 no Part 0.
"""
import os, sys, json, time, hashlib
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
import mg8
from mg8 import J, R, L, P, NAMES, TAU0, H
from runner import seed_for
from measure import smooth_periodic, detect, analyze_run

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results", "MG8")
RUNS, FIELDS, SNAP, FULL = (os.path.join(OUT, d) for d in ("runs", "fields", "snap", "full"))
LEDGER = os.path.join(OUT, "ledger.json")
CUTOFF_CORE_H = 20.0
CAP = 64
KMAX = 6
T_B1, T_CONT, T_LCHK6, T_B2 = 300, 600, 50, 500
OMEGA = 1000.0
N_B2 = 8
CHECK_EVERY = 50


# ------------------------------------------------------------------ job list (proposal §3.2)
def jobs():
    b1 = []
    for k in range(KMAX + 1):
        T = T_CONT if k == KMAX else T_B1
        b1.append(dict(kind="det", sched=J, k=k, rep=0, T=T))
        for j in range(3):
            b1.append(dict(kind="det", sched=R, k=k, rep=j, T=T if j == 0 else T_B1))
        if k in (0, KMAX):
            b1.append(dict(kind="det", sched=L, k=k, rep=0, T=T_B1 if k == 0 else T_LCHK6))
    b1.sort(key=lambda d: (d["k"] != 0, d["k"]))         # k = 0 first
    b2 = [dict(kind="sto", sched=s, k=0, rep=r, T=T_B2) for s in (J, R, P) for r in range(N_B2)]
    b2.append(dict(kind="sto", sched=L, k=0, rep=0, T=T_B2))
    return b1 + b2


def key(job):
    return f"{job['kind']}_{NAMES[job['sched']]}_k{job['k']}_r{job['rep']}"


def seed(job):
    if job["kind"] == "det":
        return seed_for(8, 1, job["k"], 0, job["rep"])     # L uses R rep 0's seed (rep = 0)
    a = {J: 0, R: 1, P: 2, L: 1}[job["sched"]]               # L mirror run uses R's first seed
    return seed_for(8, 0, a, 0, job["rep"])


# ------------------------------------------------------------------ one deterministic run
def det_run(job):
    Pm = mg8.params()
    ref, c = mg8.reference_start()
    S0 = ref["S0"]
    x = ref["x"]
    if job["sched"] == L:
        c = c[:, ::-1]
    tau = TAU0 / 2 ** job["k"]
    spu = mg8.steps_per_unit(tau)
    G, X, Y = (np.ascontiguousarray(c[i], dtype=float).copy() for i in range(3))
    sd = seed(job)
    t_rec, cen, ndet, ckpt = [], [], [], {}
    clips = np.zeros(3, dtype=np.int64)
    t0 = time.time()
    for u in range(1, job["T"] + 1):
        clips += mg8.run_det(G, X, Y, job["sched"], Pm["a"], Pm["b"], Pm["g"], Pm["p"], Pm["q"], Pm["dx"], Pm["dy"],
                             H, tau, spu, sd if u == 1 else -1)
        phiX_s = smooth_periodic(X - S0[1], H)
        t_rec.append(u * spu * tau)
        cen.append(mg8.centre(phiX_s, x, H))
        ndet.append(int(len(detect(phiX_s, x, H, ref["thr"], ref["lam"]))))
        if u % CHECK_EVERY == 0:
            ckpt[u] = np.stack([G - S0[0], X - S0[1], Y - S0[2]])
        if not np.all(np.isfinite(X)):
            break
    wall = time.time() - t0
    np.savez_compressed(os.path.join(FIELDS, key(job) + ".npz"),
                        **{f"t{u}": v for u, v in ckpt.items()})
    return dict(job=job, key=key(job), seed=int(sd), tau=tau, steps_per_unit=spu, t=t_rec, centre=cen, n_dets=ndet,
                clips=clips.tolist(), wall=wall, finite=bool(np.all(np.isfinite(X))), platform=sys.platform,
                numpy=np.__version__)


# ------------------------------------------------------------------ one stochastic run
def track_fail(dets, pos, Ld=100.0):
    """TRACK-FAIL (§3.3): the x0 track ends while a detection lies within 3 units of its last position."""
    last = None
    for s in range(len(pos)):
        if np.isfinite(pos[s]):
            last = pos[s]
        elif last is not None:
            d = dets[s]
            if len(d) and np.min(np.minimum(np.abs(d - last) % Ld, Ld - np.abs(d - last) % Ld)) < 3.0:
                return True, s
            return False, s
    return False, None


def sto_run(job):
    Pm = mg8.params()
    ref, c = mg8.reference_start()
    S0, x = ref["S0"], ref["x"]
    n = [np.rint(OMEGA * c[i]).astype(np.int64) for i in range(3)]
    if job["sched"] == L:
        n = [v[::-1].copy() for v in n]
    tau = TAU0
    spu = mg8.steps_per_unit(tau)
    sd = seed(job)
    t0 = time.time()
    rG, rX, rY, rs, clips = mg8.run_stoch(n[0], n[1], n[2], job["sched"], OMEGA, Pm["a"], Pm["b"], Pm["g"], Pm["p"],
                                          Pm["q"], Pm["dx"], Pm["dy"], H, tau, job["T"] * spu, spu, sd)
    wall = time.time() - t0
    times = rs * tau
    out = dict(job=job, key=key(job), seed=int(sd), tau=tau, steps_per_unit=spu, wall=wall, clips=clips.tolist(),
               clip_frac=(clips / (500 * job["T"] * spu)).tolist(), platform=sys.platform, numpy=np.__version__)
    full = os.path.join(FULL, key(job) + ".npz")
    np.savez_compressed(full, G=rG, X=rX, Y=rY, times=times)
    out["full_sha256"] = hashlib.sha256(open(full, "rb").read()).hexdigest()
    arr = np.stack([rG, rX, rY]).astype(np.int64)
    if job["sched"] == L:
        arr = arr[:, :, ::-1]                                # map back to R's orientation for the mirror check
    out["counts_sha256"] = hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()
    if job["sched"] == L:
        return out                                           # mirror check only
    a = analyze_run(rG, rX, rY, times, OMEGA, H, S0, ref["A_s"], ref["lam"], t_discard=150.0, K=7, T=500.0,
                    min_snap_per_block=1)
    phiX_s = smooth_periodic(rX / OMEGA - S0[1], H)
    dets = [detect(phiX_s[s], x, H, ref["thr"], ref["lam"]) for s in range(len(times))]
    pos = np.array(a["positions"], float)
    tf, tf_at = track_fail(dets, pos)
    phys = [e for e in a["events"] if e[0] in ("loss", "replication", "merge")]
    eligible = (not tf) and len(phys) == 0
    v = None
    if eligible:
        w = (times >= 150 - 1e-9) & np.isfinite(pos)
        v = mg8.slope(times[w], mg8.unwrap(pos[w], 100.0)) if w.sum() >= 2 else None
    nb = np.array(a["n_per_block"], float)
    prof = None
    if eligible and nb.sum() > 0:
        bl = np.nan_to_num(np.array(a["blocks"]))
        prof = (np.tensordot(nb, bl, axes=(0, 0)) / nb.sum()).tolist()
    np.savez_compressed(os.path.join(SNAP, key(job) + ".npz"), phiX_s=phiX_s[4::5].astype(np.float32),
                        times=times[4::5])
    out.update(events=a["events"], track_fail=tf, track_fail_at=tf_at, eligible=eligible, v=v, profile=prof,
               n_per_block=a["n_per_block"], positions=a["positions"], dets=[d.tolist() for d in dets],
               farvar=a["farvar"], n_dets=a["n_dets"])
    return out


def do(job):
    r = det_run(job) if job["kind"] == "det" else sto_run(job)
    with open(os.path.join(RUNS, key(job) + ".json"), "w") as f:
        json.dump(r, f)
    return key(job), r["wall"]


# ------------------------------------------------------------------ main
def ledger():
    return json.load(open(LEDGER)) if os.path.exists(LEDGER) else dict(runs={}, core_h=0.0)


def main():
    for d in (RUNS, FIELDS, SNAP, FULL):
        os.makedirs(d, exist_ok=True)
    dec = os.path.join(OUT, "part0", "decisions.json")
    if not os.path.exists(dec) or not json.load(open(dec)).get("proceed"):
        print("Part 0 decisions missing or proceed = false"); sys.exit(5)
    mg8.reference_start()                                    # cache the reference in the main process
    led = ledger()
    todo = [j for j in jobs() if not os.path.exists(os.path.join(RUNS, key(j) + ".json"))]
    assert len(jobs()) <= CAP
    nw = int(os.environ.get("MODELG_WORKERS", "2"))
    with Pool(nw) as pool:
        for i in range(0, len(todo), nw):
            if led["core_h"] >= CUTOFF_CORE_H:
                print(f"BUDGET-SHORT: ledger {led['core_h']:.2f} core-h >= {CUTOFF_CORE_H}"); sys.exit(2)
            for k_, w in pool.imap_unordered(do, todo[i:i + nw]):
                led["runs"][k_] = w
                led["core_h"] = sum(led["runs"].values()) / 3600.0
                json.dump(led, open(LEDGER, "w"), indent=1)
                print(f"{k_}: {w:.0f} s; ledger {led['core_h']:.2f} core-h", flush=True)
    print("done")


if __name__ == "__main__":
    main()
