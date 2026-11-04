"""MG7 independent recompute (proposal §8). Uses MG6's independent audit routines (src/recompute_mg6.py: candidate
regeneration, brute-force detection, classification; written separately from mg6 / run_mg6 / evaluate / evolve) and
nothing from the MG7 driver or evaluator. From the committed records only (plus the MG4c start graphs and files) it
  1. checks every run's start file against the endpoint manifest and that the run used the requested b with no shift;
  2. regenerates each run's candidate list from its seed and checks the recorded SHA-256 (shared across b and arms);
  3. rebuilds the node set at every check from the start nodes, lifetimes, candidates and the accepted bitmap, checks
     that no ephemeral candidate was accepted, and recomputes N(t), N(<R) and the coarse annular counts;
  4. re-detects every check from the sparse records and classifies the run;
  5. recomputes D0 window means, the endpoint contrasts, Delta_i, the t statistics and intervals, the residuals against
     the stored predictions, the feedback outcomes and the gate;
and compares all of it with the committed records and results/MG7/evaluation.json. Exits 1 on any difference.
    python src/recompute_mg7.py                     (records' platform: full check, hash included)
    python src/recompute_mg7.py --cross-platform    (elsewhere: candidate hash differences reported, rest checked)
"""
import os, sys, json, csv, hashlib
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
import recompute_mg6 as A                                     # the MG6 independent audit routines

ROOT = os.path.join(HERE, "..", "results")
OUT = os.environ.get("MG7_OUT") or os.path.join(ROOT, "MG7")
STARTS = ["D0_p0", "D0_p1", "D1_p0", "D1_p1"]
B_END = (27.5, 28.5)
T975 = 3.182
CROSS_PLATFORM = "--cross-platform" in sys.argv


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def stats(c):
    c = np.asarray(c, float)
    m, s = c.mean(), c.std(ddof=1)
    return [*c, m, s, A.tstat(c), m - T975 * s / 2, m + T975 * s / 2]


def main():
    bad, notes = [], []
    ev_ = json.load(open(os.path.join(OUT, "evaluation.json")))
    w0, w1 = ev_["window"]
    man = {r["path"]: r["sha256"] for r in csv.DictReader(open(os.path.join(ROOT, "MG7", "endpoint_manifest.csv")))}
    Y, P_start, Y_P, single_w, status, cls_of = {}, {}, {}, {}, {}, {}
    for path in sorted(os.listdir(os.path.join(OUT, "runs"))):
        rec = json.load(open(os.path.join(OUT, "runs", path)))
        key, sp = rec["key"], rec["spec"]
        status[key] = rec["status"]
        b = sp["b"]
        fp = f"MG4c/fields/b{b:.6f}_D{sp['slot']}_p{sp['pl']}.npz"
        if (rec["start_field"]["path"] != fp or rec["start_field"]["sha256"] != man.get(fp)
                or sha(os.path.join(ROOT, fp)) != man.get(fp) or rec["start_field"]["shift_applied"]
                or rec["b"] != b or rec["P_b"] != b or not key.startswith(f"{sp['arm']}_b{b:.2f}_")):
            bad.append(f"{key}: start field, b or key")
        if rec["status"] != "COMPLETE":
            continue
        with np.load(os.path.join(OUT, "checks", key + ".npz")) as z:
            d = {k: z[k] for k in z.files}
        Lx, Ly, nu, delta, t_end = rec["Lx"], rec["Ly"], rec["nu"], rec["delta"], rec["t_end"]
        cand, hsh = A.regen(18901, nu, t_end, Lx, Ly, rec["replay"]["seed"])
        if hsh != rec["replay"]["sha256"]:
            if not CROSS_PLATFORM:
                bad.append(f"{key}: candidate hash")
                continue
            notes.append(key)
        acc = np.unpackbits(d["accepted_bits"])[:int(d["n_candidates"])].astype(bool)
        k_of = np.ceil(cand["tb"] / delta - 1e-12)
        if np.any(acc & (cand["tb"] + cand["ell"] <= k_of * delta)):
            bad.append(f"{key}: ephemeral candidate accepted")
        p0 = A.start_pos(sp["slot"], sp["pl"])
        lists, ok_counts = [], True
        for i, t in enumerate(d["t"]):
            a, bb = d["off"][i], d["off"][i + 1]
            lists.append(A.detect(d["pos"][a:bb], d["s"][a:bb], rec["theta"], Lx, Ly))
            tp = t - delta if t > 0 else 0.0
            alive0 = cand["life0"] > tp if t > 0 else np.ones(len(p0), bool)
            born = acc & (cand["tb"] <= tp) & (cand["tb"] + cand["ell"] > tp) if t > 0 else np.zeros(len(acc), bool)
            pos = np.concatenate([p0[alive0], np.column_stack([cand["x"][born], cand["y"][born]])])
            if np.isnan(d["m_N"][i]):                       # no particle, so no centred measurement was stored
                if lists[-1] or int(d["n_parts"][i]) != 0:
                    bad.append(f"{key}: NaN measurement at t = {t} although a particle is present")
            elif len(pos) != int(d["m_N"][i]):
                ok_counts = False
            if lists[-1]:
                c = lists[-1][0]["centre"]
                if not np.allclose(c, d["centre"][i], atol=1e-9):
                    bad.append(f"{key}: centre at t = {t}")
                r = np.hypot(A.wrap(pos[:, 0] - c[0], Lx), A.wrap(pos[:, 1] - c[1], Ly))
                NR = np.array([(r < R_).sum() for R_ in A.RADII])
                ann = np.array([((r >= lo) & (r < hi)).sum() for lo, hi in zip(A.ANN[:-1], A.ANN[1:])])
                if not (np.array_equal(NR, d["m_N_R"][i].astype(int))
                        and np.array_equal(ann, d["m_ann_count"][i].astype(int))):
                    ok_counts = False
        if not ok_counts:
            bad.append(f"{key}: rebuilt node counts differ")
        if [len(x) for x in lists] != d["n_parts"].tolist():
            bad.append(f"{key}: particle counts")
        cls_of[key] = A.classify(lists, Lx, Ly)
        if cls_of[key] != ev_["runs"][key]["cls"]["cls"]:
            bad.append(f"{key}: class {cls_of[key]}")
        rho0 = 18901 / (Lx * Ly)
        t = d["t"]
        win = (t >= w0 - 1e-9) & (t <= w1 + 1e-9)
        single_w[key] = bool(np.all(d["n_parts"][win] == 1))
        P_start[key] = float(d["m_P0"][0, 1])
        if single_w[key]:
            Y[key] = float(np.mean(d["m_N_R"][win, 1] - rho0 * np.pi * 1.5 ** 2))
            Y_P[key] = float(np.mean(d["m_P0"][win, 1]))
            if abs(Y[key] - ev_["runs"][key]["Y_D0"]) > 1e-9 * (1 + abs(Y[key])):
                bad.append(f"{key}: window mean")
        print(f"{key}: {cls_of[key]}, single in window {single_w[key]}", flush=True)
    k = lambda arm, b, s: f"{arm}_b{b:.2f}_{s}"
    allk = [k(a, b, s) for a in ("CP2", "SH2") for b in B_END for s in STARTS]
    unaudited = [x for x in allk if status.get(x) == "COMPLETE" and x not in single_w]
    if unaudited:                                             # strict mode: a candidate hash failed, run not audited
        gate = None
        print(f"not audited (candidate hash): {unaudited}; gate not recomputed")
    elif any(status.get(x) != "COMPLETE" for x in allk):
        gate = "INCOMPLETE"
    elif not all(single_w[x] for x in allk):
        gate = "NOT SCORED (feedback)"
    else:
        c = {b: np.array([Y[k("CP2", b, s)] - Y[k("SH2", b, s)] for s in STARTS]) for b in B_END}
        D = c[28.5] - c[27.5]
        gate = "LARGER AT HIGHER b" if (np.all(D < 0) and A.tstat(D) <= -A.TCRIT) else "NOT DEMONSTRATED"
        st = ev_["primary"]
        if not np.allclose(stats(D), list(st["values"]) + [st["mean"], st["sd"], st["t"]] + list(st["ci95"]),
                           rtol=1e-9, atol=1e-9):
            bad.append("primary statistics")
        for b in B_END:
            s2 = ev_["contrasts"][f"b{b:.2f}"]
            if not np.allclose(stats(c[b]), list(s2["c"]) + [s2["mean"], s2["sd"], s2["t"]] + list(s2["ci95"]),
                               rtol=1e-9, atol=1e-9):
                bad.append(f"contrast b {b}")
        for name, P in (("start_field", P_start), ("in_run", Y_P)):
            pc = {b: np.array([P[k("CP2", b, s)] - P[k("SH2", b, s)] for s in STARTS]) for b in B_END}
            e = D - (pc[28.5] - pc[27.5])
            m, s = e.mean(), e.std(ddof=1)
            rs = ev_["residuals"][name]
            if not np.allclose(list(e) + [m, m - T975 * s / 2, m + T975 * s / 2],
                               list(rs["e"]) + [rs["mean"]] + list(rs["ci95"]), rtol=1e-9, atol=1e-9):
                bad.append(f"residuals {name}")
        print(f"Delta_i = {np.round(D, 3).tolist()}, t = {A.tstat(D):.4f} (stored {st['t']:.4f})")
    if gate is not None:
        print(f"gate: {gate} ({'agrees' if gate == ev_['gate'] else 'DIFFERS from ' + ev_['gate']})")
        if gate != ev_["gate"]:
            bad.append("gate")
    for b in B_END:
        ks = [k("CP2", b, s) for s in STARTS]
        if all(status.get(x) == "COMPLETE" and x in cls_of for x in ks):
            fb = "SURVIVES" if sum(cls_of[x] == "single" for x in ks) >= 3 else "CHANGED"
            if fb != ev_["feedback"][f"b{b:.2f}"]["outcome"]:
                bad.append(f"feedback b {b}")
            print(f"feedback b {b}: {fb}")
    if notes:
        print(f"candidate hash differs on this platform for {len(notes)} run(s) (--cross-platform): {notes}")
    print("Everything agrees." if not bad else f"DIFFER: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
