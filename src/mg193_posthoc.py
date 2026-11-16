"""MG193 post-hoc analyses (not registered tests; written after the registered run was seen).

1. Re-derive the time-check verdicts from the saved histories in window.json, applying the registered wording of
   spec/MG193_existence_window.md section 6. That wording says "another deep minimum" for outside and "core depth
   within 1%" for inside. The driver instead required exactly one minimum below -1, but the stationary soliton has
   three: the core plus two tail dips. The start count is measured here from the stationary profiles.
2. Follow the diagonal strip of existence to a = 12, 13, 15, 16. Step along b = b* + k (a - 14), where k is the strip
   slope from the two registered a-windows. Then map the b-window at each a with the same continuation code. This is
   descriptive context only; the registered "could not reach" results for those a values stand.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mg192_mobility as m  # noqa: E402
import mg193_window as w  # noqa: E402

OUT = os.path.join(m.ROOT, "results", "MG193_existence_window")


def minima_count(u, a, b):
    dX = u[1] - m.hom(a=a, b=b)[1]
    return int(np.sum((dX < -1.0) & (dX <= np.roll(dX, 1)) & (dX <= np.roll(dX, -1))))


def main():
    t0 = time.time()
    rec = json.load(open(os.path.join(OUT, "window.json")))
    U14 = np.load(os.path.join(m.ROOT, "results", "MG192_gravity_mobility", "soliton_adjoint_L100_N1024.npz"))["U"]
    ctx = w.Ctx()
    out = dict(note="post hoc; not registered tests", verdicts={}, strip={})

    # 1. verdicts from saved histories with the registered wording
    start_minima = {}
    for name, c in rec["time_checks"].items():
        a, b = c["a"], c["b"]
        if name.startswith("a_"):
            u = w.solution_at(ctx, a / 14.0 - 1.0, None, U14, 14.0, 29.0, "a") if "inside" in name else None
        else:
            u = w.solution_at(ctx, b - 29.0, None, U14, 14.0, 29.0, "b") if "inside" in name else None
        if u is not None:
            start_minima[name] = minima_count(u, a, b)
    n0 = sorted(set(start_minima.values()))
    out["stationary_start_minima"] = start_minima
    base = max(n0) if n0 else 1
    for name, c in rec["time_checks"].items():
        d0, last = c["start_depth"], c["history"][-1]
        if "inside" in name:
            v = "PERSISTS" if abs(last["depth"] - d0) <= 0.01 * abs(d0) else "FAILS"
        elif last["min_depth"] > -1.0:
            v = "DISSOLVED"
        elif last["deep_minima"] > base or abs(last["depth"] - d0) > 0.1 * abs(d0):
            v = "CHANGED"
        else:
            v = "INCONCLUSIVE"
        out["verdicts"][name] = dict(driver=c["verdict"], registered_wording=v, start_depth=d0,
                                     final_depth=last["depth"], final_minima=last["deep_minima"])
        print(name, c["verdict"], "->", v, flush=True)
    rp = rec["recommended_point"]
    ustar = w.solution_at(ctx, rp["b"] - 29.0, None, U14, 14.0, 29.0, "b")
    pr = rp["perturbed_run"]
    out["recommended_perturbed"] = dict(start_minima=minima_count(ustar, 14.0, rp["b"]), driver_ok=pr["ok"],
                                        registered_wording_ok=bool(abs(pr["final_min_depth"] - pr["start_depth"])
                                                                   <= 0.01 * abs(pr["start_depth"])))
    print("recommended perturbed", out["recommended_perturbed"], flush=True)

    # 2. follow the strip
    al = rec["a_lines"]
    keys = sorted(al.keys())
    c29 = 0.5 * (al[keys[0]]["a_low"] + al[keys[0]]["a_high"])
    cst = 0.5 * (al[keys[1]]["a_low"] + al[keys[1]]["a_high"])
    b29, bst = 29.0, rec["a14"]["b_star"]
    k = (bst - b29) / (cst - c29)
    out["strip_slope_db_da"] = k
    print(f"strip slope db/da = {k:.4f}", flush=True)
    for target in (13.0, 12.0, 15.0, 16.0):
        u = U14.copy()
        a, b = 14.0, 29.0
        # first move to b* at a = 14
        u = w.solution_at(ctx, bst - 29.0, None, U14, 14.0, 29.0, "b")
        a, b = 14.0, bst
        n = int(round(abs(target - 14.0) / 0.01))
        ok = True
        for i in range(1, n + 1):
            a_i = 14.0 + (target - 14.0) * i / n
            b_i = bst + k * (a_i - 14.0)
            u2, r, ok, why = w.newton(ctx, u, a_i, b_i)
            if not ok:
                out["strip"][str(target)] = dict(reached=False, failed_at=[a_i, b_i], why=why)
                print(f"a={target}: strip lost at a={a_i:.3f} b={b_i:.4f} ({why})", flush=True)
                break
            u, a, b = u2, a_i, b_i
        if not ok:
            continue
        bc = w.turing_threshold(a)
        up_rows, up = w.continue_line(ctx, u, a, b, "b", +1, 0.02, 5.0, stop_value=bc, log=lambda s: None)
        dn_rows, dn = w.continue_line(ctx, u, a, b, "b", -1, 0.02, 5.0, log=lambda s: None)
        blo = b + 0.5 * sum(dn["bracket"]) if "bracket" in dn else b + dn.get("last_ok_s", dn.get("last_s", 0))
        bhi = b + 0.5 * sum(up["bracket"]) if "bracket" in up else b + up.get("last_ok_s", up.get("last_s", 0))
        out["strip"][str(target)] = dict(reached=True, start_b=b, b_c=bc, b_low=blo, b_high=bhi,
                                         width=bhi - blo, edge_low=dn["edge"], edge_high=up["edge"])
        print(f"a={target}: b-window [{blo:.4f}, {bhi:.4f}] width {bhi - blo:.4f} edges {dn['edge']}/{up['edge']}"
              f" b_c={bc:.4f}", flush=True)
        json.dump(out, open(os.path.join(OUT, "posthoc.json"), "w"), indent=1, default=float)
    out["wall_seconds"] = time.time() - t0
    json.dump(out, open(os.path.join(OUT, "posthoc.json"), "w"), indent=1, default=float)
    print("done", out["wall_seconds"], flush=True)


if __name__ == "__main__":
    main()
