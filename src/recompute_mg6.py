"""MG6 independent recompute (proposal §8, audit scope). Written separately from mg6 / run_mg6 / evaluate_mg6 / evolve:
from the committed records only (plus the MG4c start graphs), it
  1. regenerates each run's candidate list from its seed with its own code and checks the recorded SHA-256;
  2. rebuilds the node set at every check from the start nodes, the initial lifetimes, the candidates and the accepted
     bitmap (a node is present in M_(k-1) iff born at or before t_(k-1) and dying after it), checks that no ephemeral
     candidate was accepted, and recomputes N(t), N(<R) and the coarse annular counts around the recorded centre;
  3. re-detects every check from the sparse records (brute-force distances, union-find) and classifies the run;
  4. recomputes D0, Dc, the run-level window means, the paired contrasts, the t statistics and every §5 outcome;
and compares all of it with the committed records and results/MG6/evaluation.json. Exits 1 on any difference.
Quadrature quantities (P0, Pc, M-, M+, radial phi_G) are not rebuilt here (they need the PC-only fields); they were
cross-checked in-run.
    python src/recompute_mg6.py                     (on the records' platform: full check, hash included)
    python src/recompute_mg6.py --cross-platform    (elsewhere: hash differences reported, everything else checked)
"""
import os, sys, json, hashlib
import numpy as np

HERE = os.path.dirname(__file__)
OUT = os.environ.get("MG6_OUT") or os.path.join(HERE, "..", "results", "MG6")
RADII = (1.0, 1.5, 2.0, 3.0, 5.0, 8.0)
ANN = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
CAND, MERGE, CENT, STEP = 1.0, 1.5, 0.6, 2.0
TCRIT = 2.353
F_BOUND = float(np.exp(0.4))


def wrap(d, L):
    return d - L * np.round(d / L)


def regen(N0, nu, t_end, Lx, Ly, seed):
    rng = np.random.default_rng(seed)
    life0 = rng.exponential(1.0 / nu, N0)
    n = int(rng.poisson(nu * N0 * F_BOUND * t_end))
    tb = np.sort(rng.random(n) * t_end)
    x, y, u, ell = rng.random(n) * Lx, rng.random(n) * Ly, rng.random(n), rng.exponential(1.0 / nu, n)
    h = hashlib.sha256(np.array([N0], "<i8").tobytes())
    for a in (life0, tb, x, y, u, ell):
        h.update(np.ascontiguousarray(a, "<f8").tobytes())
    return dict(life0=life0, tb=tb, x=x, y=y, u=u, ell=ell), h.hexdigest()


def detect(pos, s, thr, Lx, Ly):
    cands = []
    for i in np.flatnonzero(s < thr):
        dd = np.hypot(wrap(pos[:, 0] - pos[i, 0], Lx), wrap(pos[:, 1] - pos[i, 1], Ly))
        nb = np.flatnonzero(dd <= CAND)
        m = s[nb].min()
        if s[i] < m or (s[i] == m and i == nb[s[nb] == m].min()):
            cands.append(i)
    par = list(range(len(cands)))

    def root(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    for a in range(len(cands)):
        for b in range(a + 1, len(cands)):
            p, q = pos[cands[a]], pos[cands[b]]
            if np.hypot(wrap(p[0] - q[0], Lx), wrap(p[1] - q[1], Ly)) <= MERGE:
                par[root(a)] = root(b)
    groups = {}
    for a, c in enumerate(cands):
        groups.setdefault(root(a), []).append(c)
    out = []
    for mem in groups.values():
        mem = np.array(mem)
        i0 = mem[np.argmin(s[mem])]
        dd = np.hypot(wrap(pos[:, 0] - pos[i0, 0], Lx), wrap(pos[:, 1] - pos[i0, 1], Ly))
        nb = np.flatnonzero(dd <= CENT)
        w = np.maximum(thr - s[nb], 0.0)
        dv = np.column_stack([wrap(pos[nb, 0] - pos[i0, 0], Lx), wrap(pos[nb, 1] - pos[i0, 1], Ly)])
        c = pos[i0] + (w[:, None] * dv).sum(0) / w.sum()
        out.append(dict(centre=np.mod(c, (Lx, Ly)), depth=float(s[i0])))
    return sorted(out, key=lambda p: p["depth"])


def classify(lists, Lx, Ly):
    n = [len(x) for x in lists]
    if n[0] != 1:
        return "unresolved"
    n = n[1:]
    if any(n[i] >= 2 and n[i + 1] >= 2 for i in range(len(n) - 1)):
        return "multiplied"
    if any(n[i] == 0 and n[i + 1] == 0 for i in range(len(n) - 1)):
        return "dissolved"
    if not all(m == 1 for m in n):
        return "unresolved"
    cur = np.array(lists[0][0]["centre"])
    for x in lists[1:]:
        st = np.array([wrap(x[0]["centre"][0] - cur[0], Lx), wrap(x[0]["centre"][1] - cur[1], Ly)])
        if np.hypot(*st) > STEP:
            return "unresolved"
        cur = np.array(x[0]["centre"])
    return "single"


def start_pos(slot, pl):
    sys.path.insert(0, HERE)
    import run_mg4c as r4                                   # the MG4c start graph (positions only)
    dec = json.load(open(r4.DEC))
    g, _, _ = r4.graph_for(r4.job(28.0, "D", slot, pl), dec)
    return np.asarray(g["pos"], float)


def tstat(c):
    c = np.asarray(c, float)
    s = c.std(ddof=1)
    return c.mean() / (s / 2.0) if s > 0 else (-np.inf if c.mean() < 0 else np.inf)


CROSS_PLATFORM = "--cross-platform" in sys.argv
# MG022: on another platform the candidate hash can differ, because numpy's exponential sampler falls back to the C
# library's log for a small fraction of draws (Windows and Linux differ in the last bit). With --cross-platform a hash
# difference is reported and the rest of the audit runs on the regenerated lists; the node-set rebuild then shows whether
# any difference matters. The full check (hash included) is run on the platform that produced the records.


def main():
    bad, notes = [], []
    ev_ = json.load(open(os.path.join(OUT, "evaluation.json")))
    w0, t_end_w = ev_["window"]
    Y, Yc, single_w, status = {}, {}, {}, {}
    for path in sorted(os.listdir(os.path.join(OUT, "runs"))):
        rec = json.load(open(os.path.join(OUT, "runs", path)))
        key = rec["key"]
        status[key] = rec["status"]
        if rec["status"] != "COMPLETE":
            continue
        with np.load(os.path.join(OUT, "checks", key + ".npz")) as z:
            d = {k: z[k] for k in z.files}
        Lx, Ly, nu, delta, t_end = rec["Lx"], rec["Ly"], rec["nu"], rec["delta"], rec["t_end"]
        sp = rec["spec"]
        cand, hsh = regen(18901, nu, t_end, Lx, Ly, rec["replay"]["seed"])
        if hsh != rec["replay"]["sha256"]:
            if not CROSS_PLATFORM:
                bad.append(f"{key}: candidate hash")
                continue
            notes.append(key)
        acc = np.unpackbits(d["accepted_bits"])[:int(d["n_candidates"])].astype(bool)
        k_of = np.ceil(cand["tb"] / delta - 1e-12)                 # birth tick index (tb in ((k-1)Delta, kDelta])
        eph = cand["tb"] + cand["ell"] <= k_of * delta
        if np.any(acc & eph):
            bad.append(f"{key}: ephemeral candidate accepted")
        p0 = start_pos(sp["slot"], sp["pl"])
        lists, ok_counts = [], True
        for i, t in enumerate(d["t"]):
            a, b = d["off"][i], d["off"][i + 1]
            lists.append(detect(d["pos"][a:b], d["s"][a:b], rec["theta"], Lx, Ly))
            tp = t - delta if t > 0 else 0.0                       # the check at t_k sees M_(k-1)
            alive0 = cand["life0"] > tp if t > 0 else np.ones(len(p0), bool)
            born = acc & (cand["tb"] <= tp) & (cand["tb"] + cand["ell"] > tp) if t > 0 else np.zeros(len(acc), bool)
            pos = np.concatenate([p0[alive0], np.column_stack([cand["x"][born], cand["y"][born]])])
            if len(pos) != int(d["m_N"][i]):
                ok_counts = False
            if lists[-1]:
                c = lists[-1][0]["centre"]
                if not np.allclose(c, d["centre"][i], atol=1e-9):
                    bad.append(f"{key}: centre at t = {t}")
                r = np.hypot(wrap(pos[:, 0] - c[0], Lx), wrap(pos[:, 1] - c[1], Ly))
                NR = np.array([(r < R_).sum() for R_ in RADII])
                ann = np.array([((r >= lo) & (r < hi)).sum() for lo, hi in zip(ANN[:-1], ANN[1:])])
                if not (np.array_equal(NR, d["m_N_R"][i].astype(int)) and
                        np.array_equal(ann, d["m_ann_count"][i].astype(int))):
                    ok_counts = False
        if not ok_counts:
            bad.append(f"{key}: rebuilt node counts differ")
        if [len(x) for x in lists] != d["n_parts"].tolist():
            bad.append(f"{key}: particle counts")
        cls = classify(lists, Lx, Ly)
        if cls != ev_["runs"][key]["cls"]["cls"]:
            bad.append(f"{key}: class {cls}")
        A = Lx * Ly
        rho0 = 18901 / A
        t = d["t"]
        win = (t >= w0 - 1e-9) & (t <= t_end_w + 1e-9)
        single_w[key] = bool(np.all(d["n_parts"][win] == 1))
        if single_w[key]:
            NR15 = d["m_N_R"][win, 1]
            Y[key] = float(np.mean(NR15 - rho0 * np.pi * 1.5 ** 2))
            Yc[key] = float(np.mean(NR15 - d["m_N"][win] / A * np.pi * 1.5 ** 2))
            for mine, theirs in ((Y[key], ev_["runs"][key]["Y_D0"]), (Yc[key], ev_["runs"][key]["Y_Dc"])):
                if abs(mine - theirs) > 1e-9 * (1 + abs(theirs)):
                    bad.append(f"{key}: window mean")
        print(f"{key}: {cls}, single in window {single_w[key]}", flush=True)
    starts = ["D0_p0", "D0_p1", "D1_p0", "D1_p1"]
    prim = [f"{a}_{s}" for a in ("CP2", "SH2", "DEC") for s in starts]
    if any(status.get(k) != "COMPLETE" for k in prim):
        gate = "INCOMPLETE"
    elif not all(single_w[k] for k in prim):
        gate = "NOT SCORED (feedback)"
    else:
        ok = True
        for null in ("SH2", "DEC"):
            c = [Y[f"CP2_{s}"] - Y[f"{null}_{s}"] for s in starts]
            ok = ok and all(x < 0 for x in c) and tstat(c) <= -TCRIT
        gate = "DEMONSTRATED" if ok else "NOT DEMONSTRATED"
    print(f"gate: {gate} ({'agrees' if gate == ev_['gate'] else 'DIFFERS from ' + ev_['gate']})")
    if gate != ev_["gate"]:
        bad.append("gate")
    for a in ("CP2", "CP4"):
        ks = [f"{a}_{s}" for s in starts]
        if all(status.get(k) == "COMPLETE" for k in ks):
            n_single = sum(ev_["runs"][k]["cls"]["cls"] == "single" for k in ks)
            fb = "SURVIVES" if n_single >= 3 else "CHANGED"
            if fb != ev_["feedback"][a]["outcome"]:
                bad.append(f"feedback {a}")
            print(f"feedback {a}: {fb}")
    if notes:
        print(f"candidate hash differs on this platform for {len(notes)} run(s) (--cross-platform): {notes}")
    print("Everything agrees." if not bad else f"DIFFER: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
