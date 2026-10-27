"""MG004: run orchestration for Phase 2 (protocol v1.3, frozen).

- seed_for(...): SeedSequence stream keys (protocol §7).
- reference(b, h): lattice-ODE soliton reference (phi*, A*, lambda*), cached in results/MG003/ref_profiles.npz.
- one_run(...): one tau-leap trajectory (soliton or homogeneous start) + §5 analysis.
- run_cell(...): attempts runs in stream order (2 in parallel) until n_eligible eligible or the attempt cap.
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from lattice import tau_leap_euler_1d, lattice_ode_1d, steady_state, P17
from measure import analyze_run, analyze_homogeneous, smooth_periodic
from mg_continuum import zeros

ENTROPY = 20261027
REF_FILE = os.path.join(os.path.dirname(__file__), "..", "results", "MG003", "ref_profiles.npz")


def seed_for(gate, b_idx, om_idx, stage, replicate, segment=0):
    ss = np.random.SeedSequence(entropy=ENTROPY, spawn_key=(gate, b_idx, om_idx, stage, replicate, segment))
    return int(ss.generate_state(1, dtype=np.uint32)[0])


def reference(b, h, L=100.0):
    key = f"h{h}_b{b}"
    prof = dict(np.load(REF_FILE)) if os.path.exists(REF_FILE) else {}
    if key not in prof:
        N = int(round(L / h))
        x, pG, pX, pY, mdt = lattice_ode_1d({**P17, "b": b}, N, h, 400.0, seed_scale=1.5)
        prof[key] = np.vstack([x, pG, pX, pY])
        np.savez(REF_FILE, **prof)
    x, pG, pX, pY = prof[key]
    N = x.size
    i0 = N // 2
    z = zeros(x[i0:], pX[i0:])
    lam = 2 * float(np.mean(np.diff(z[z < 30])[1:]))
    return dict(x=x, phi=np.vstack([pG, pX, pY]), A=float(pX[i0]), lam=lam,
                S0=np.array(steady_state({**P17, "b": b})))


def one_run(args):
    """args: dict(b, Om, h, tau, T, start ('soliton'|'homog'), seed_scale, rng, kind ('mg2'|'homog'))."""
    b, Om, h, tau, T = args["b"], args["Om"], args["h"], args["tau"], args["T"]
    P = {**P17, "b": b}
    ref = reference(b, h)
    S0 = ref["S0"]
    N = ref["x"].size
    if args.get("start", "soliton") == "soliton":
        c = S0[:, None] + ref["phi"]
        nG, nX, nY = (np.rint(Om * c[k]).astype(np.int64) for k in range(3))
    else:
        rs = np.random.default_rng(args["rng"] ^ 0x5A5A5A5A)
        nG, nX, nY = (rs.poisson(Om * S0[k], N).astype(np.int64) for k in range(3))
    nsteps = int(round(T / tau)); rec = int(round(0.5 / tau))
    t0 = time.time()
    rG, rX, rY, rt, clips = tau_leap_euler_1d(nG, nX, nY, float(Om), P["a"], P["b"], P["g"], P["p"], P["q"],
                                              P["dx"], P["dy"], h, tau, nsteps, 0.0,
                                              float(args.get("seed_scale", 0.0)), ref["x"], rec, args["rng"])
    wall = time.time() - t0
    out = dict(rng=args["rng"], clips=int(clips), clip_frac=clips / (3 * N * nsteps), wall=wall)
    if args.get("kind", "mg2") == "mg2":
        a = analyze_run(rG, rX, rY, rt, Om, h, S0, ref["A"], ref["lam"])
        out.update(a)
    elif args.get("kind") == "homog":
        a = analyze_homogeneous(rG, rX, rY, rt, Om, h, S0, ref["A"], ref["lam"])
        out.update(dict(births=a["births"], first_birth=a["first_birth"], n_dets=a["n_dets"]))
        if args.get("keep_fields"):
            out["phi"] = np.stack([rG / Om - S0[0], rX / Om - S0[1], rY / Om - S0[2]], axis=1).astype(np.float32)
            out["times"] = rt
    elif args.get("kind") == "stats":
        out["fields"] = np.stack([rG, rX, rY], axis=1)
        out["times"] = rt
    return out


def run_cell(pool, base, key, n_eligible, cap, start_rep=0, log=None):
    """Attempt runs (replicate = start_rep, start_rep+1, ...) until n_eligible eligible or cap attempted.
    base: dict of one_run args without rng; key: (gate, b_idx, om_idx, stage). Returns (eligible, all)."""
    elig, allr = [], []
    rep = start_rep
    while len(elig) < n_eligible and (rep - start_rep) < cap:
        nb = min(2, n_eligible - len(elig), cap - (rep - start_rep))
        batch = [dict(base, rng=seed_for(*key, replicate=rep + i)) for i in range(nb)]
        res = pool.map(one_run, batch)
        for i, r in enumerate(res):
            r["replicate"] = rep + i
            allr.append(r)
            if r.get("eligible"):
                elig.append(r)
            if log:
                log(f"  {key} rep {rep + i}: eligible={r.get('eligible')} events={r.get('events')} "
                    f"clip_frac={r['clip_frac']:.1e} wall={r['wall']:.0f}s")
        rep += nb
    return elig, allr
