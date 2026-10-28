"""MG005 post-gate diagnostic (descriptive; not a gate criterion).

Question: are the 'replication' events recorded in MG2 cell L_b29.16_300 real particle splits, or detector artifacts?
Method: re-run three replicates that recorded replication, with identical streams (bit-identical trajectory), keeping
the fields. For each, (i) check that the detector events reproduce; (ii) take the field 20 time units after the event
and apply the field-level confirmation rule (protocol v1.3.3 §5.4; seed-free lattice ODE to t = 300 at b = 29.16).
A real split gives >= 2 persistent particle-depth structures. Also confirms a field from one eligible (no-event) run.
Usage: python src/diag_mg2_replication.py
"""
import os, sys, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from runner import one_run, reference
from measure import analyze_run
from confirm import confirm_field

TAU = 6.666666666666667e-4
OUT = "results/MG2/diag_replication.txt"


def main():
    cell = pickle.load(open("results/MG2/cells/L_b29.16_300.pkl", "rb"))
    by_rep = {r["replicate"]: r for r in cell["runs"]}
    reps = [r["replicate"] for r in cell["runs"] if r["events"]][:3]
    reps += [next(r["replicate"] for r in cell["runs"] if r["eligible"])]
    ref = reference(29.16, 0.2)
    lines = ["MG2 post-gate diagnostic: replication events in L_b29.16_300 (b = 29.16, Omega = 300, h = 0.2)."]
    for rep in reps:
        rec = by_rep[rep]
        r = one_run(dict(b=29.16, Om=300.0, h=0.2, tau=TAU, T=100.0, start="soliton", kind="stats", rng=rec["rng"]))
        F, t = r["fields"], r["times"]
        a = analyze_run(F[:, 0], F[:, 1], F[:, 2], t, 300.0, 0.2, ref["S0"], ref["A_s"], ref["lam"])
        same = [(e[0], round(e[1], 1)) for e in a["events"]] == [(e[0], round(e[1], 1)) for e in rec["events"]]
        t_ev = rec["events"][0][1] if rec["events"] else 80.0
        tc = min(t_ev + 20.0, 100.0)
        s = int(np.argmin(np.abs(t - tc)))
        phi = F[s] / 300.0 - ref["S0"][:, None]
        c = confirm_field(phi, 29.16, 0.2, ref["A"], ref["lam"])
        lines.append(f"rep {rep}: recorded events {rec['events']}; reproduced: {same}; field at t = {t[s]:.1f} -> "
                     f"{c['cls']} (persistent structures {c['n_persistent']}, deep minima at t=300 {c['n_deep_end']}, "
                     f"min phi_X at 300 {c['min_phiX_end']:.3f}, max|dphi/dt| at 300 {c['maxdt_end']:.1e})")
        print(lines[-1], flush=True)
    open(OUT, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
