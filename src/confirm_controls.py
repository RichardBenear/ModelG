"""MG004g: control checks of the finalized confirmation rule (src/confirm.py). Ground truth for every control
state = its deterministic fate after a long continuation (t = 2000): 'persistent' if particle-depth minima
remain, 'decays' if max|phi| < 1e-3. The rule's classification (checks at t = 200..300) is compared to it.
Writes results/MG004g/confirm_controls.txt."""
import sys
import numpy as np
sys.path.insert(0, "src")
from runner import reference
from confirm import confirm_field, _deep_minima
from lattice import lattice_ode_1d, P17
from measure import fourier_shift

h, N = 0.2, 500
ref = reference(29.0, h); A, lam = ref["A"], ref["lam"]; phi = ref["phi"]
lines = ["Confirmation-rule controls. Truth = deterministic fate at t = 2000. Rule = classification from checks at t = 200..300."]
def truth(field, b):
    x, pG, pX, pY, mdt = lattice_ode_1d({**P17, "b": b}, N, h, 2000.0, seed_scale=0.0, y0=np.concatenate(field))
    nd = _deep_minima(pX, -0.5 * abs(A)).size
    m = np.max(np.abs(np.concatenate([pG, pX, pY])))
    return ("persistent" if nd > 0 else ("decays" if m < 1e-3 else "other")), nd, mdt
agree = 0; total = 0
def run(name, field, b=29.0, det=None):
    global agree, total
    r = confirm_field(field, b, h, A, lam, det_pos=det)
    tr, nd, mdt = truth(field, b)
    match = (r["cls"] == "persistent" and tr == "persistent") or (r["cls"] == "excursion" and tr == "decays")
    total += 1; agree += match
    lines.append(f"{name:58s} rule: {r['cls']:10s} (n_persist={r['n_persistent']}, deep@300={r['n_deep_end']}, maxdt@300={r['maxdt_end']:.1e}) | "
                 f"truth@2000: {tr} (deep={nd}, maxdt={mdt:.1e}) | {'AGREE' if match else 'DISAGREE'}")
    print(lines[-1], flush=True)
run("(a) stationary single soliton", phi, det=0.0)
run("(b1) two solitons, separation 20", fourier_shift(phi, h, 10.0) + fourier_shift(phi, h, -10.0), det=10.0)
run("(b2) two solitons, separation 5", fourier_shift(phi, h, 2.5) + fourier_shift(phi, h, -2.5), det=2.5)
run("(b3) two solitons, separation 6.6", fourier_shift(phi, h, 3.3) + fourier_shift(phi, h, -3.3), det=3.3)
run("(b4) two solitons, separation 9", fourier_shift(phi, h, 4.5) + fourier_shift(phi, h, -4.5), det=4.5)
run("(b5) four solitons, irregular spacing", sum(fourier_shift(phi, h, s) for s in (-33.0, -12.0, 7.0, 29.0)), det=7.0)
for frac in (0.5, 0.9, 0.99, 0.999):
    x, pG, pX, pY, _ = lattice_ode_1d({**P17, "b": 29.0}, N, h, 20.0, seed_scale=0.949 * frac)
    run(f"(c) subthreshold seed {frac} x threshold (state t=20)", np.vstack([pG, pX, pY]), det=0.0)
for frac in (1.01, 1.1):
    x, pG, pX, pY, _ = lattice_ode_1d({**P17, "b": 29.0}, N, h, 20.0, seed_scale=0.949 * frac)
    run(f"(c') suprathreshold seed {frac} x threshold (state t=20)", np.vstack([pG, pX, pY]), det=0.0)
for b in (28.928, 28.92, 28.90):
    run(f"(d) phi*(29) continued at b={b} (near/below b_s=28.929)", phi, b=b, det=0.0)
lines.append(f"Agreement with long-run truth: {agree}/{total}")
print(lines[-1])
open("results/MG004g/confirm_controls.txt", "w").write("\n".join(lines) + "\n")
