"""MG020c: UNREGISTERED centre-estimator diagnostic on the committed MG6 Part 0 sparse records (no dynamics are run).

Question (PI, 2026-11-02): is the Part 0 tick-check failure (part a, centre differences) caused by the graph's
irregular node placement adding noise to the registered centre estimate?

The registered centre (measure_graph.particles) is a (theta - s)-weighted centroid of the below-threshold nodes within
CENT_R = 0.6 of the deepest node. Here a second estimator is computed from the same stored records: a least-squares
quadratic surface fitted to the smoothed phi_X at the stored nodes within radius R of the registered centre, whose
stationary point is taken as the centre. Distances between the Delta runs are then compared under both estimators.
This is exploratory only and does not re-score Part 0 (protocol rule 2).

    python review/MG020_centre_diagnostic.py
"""
import platform
import numpy as np

P = "results/MG6/part0/check_%s.npz"


def mi(d, L):
    return d - L * np.round(d / L)


def est(z, i, R):
    a, b = z["off"][i], z["off"][i + 1]
    pos, s = z["pos"][a:b], z["s"][a:b]
    Lx, Ly = float(z["Lx"]), float(z["Ly"])
    c0 = z["centre"][i]
    d = np.c_[mi(pos[:, 0] - c0[0], Lx), mi(pos[:, 1] - c0[1], Ly)]
    k = (d ** 2).sum(1) < R * R
    x, y, v = d[k, 0], d[k, 1], s[k]
    A = np.c_[np.ones_like(x), x, y, x * x, x * y, y * y]
    co = np.linalg.lstsq(A, v, rcond=None)[0]
    H = np.array([[2 * co[3], co[4]], [co[4], 2 * co[5]]])
    off = np.linalg.solve(H, -co[1:3])
    return c0 + off, int(k.sum())


runs = {"D0.2": "d0.2_dt0.05", "D0.1": "d0.1_dt0.05", "D0.05": "d0.05_dt0.05"}
Z = {k: np.load(P % v) for k, v in runs.items()}
Lx, Ly = float(Z["D0.2"]["Lx"]), float(Z["D0.2"]["Ly"])


def dist(p, q):
    return float(np.hypot(mi(p[0] - q[0], Lx), mi(p[1] - q[1], Ly)))


print(f"platform: {platform.node()} {platform.platform()} python {platform.python_version()} numpy {np.__version__}")
print("files:", ", ".join(P % v for v in runs.values()))
print("registered tick tolerance: centre 0.15")
for R in (1.0, 1.5):
    print(f"--- quadratic fit radius {R}")
    print("chk  t   | registered d(0.2,0.1) d(0.1,0.05) | fit d(0.2,0.1) d(0.1,0.05) | nodes | fit-vs-registered shift (0.2 run)")
    for i in range(len(Z["D0.2"]["t"])):
        f = {k: est(z, i, R) for k, z in Z.items()}
        r = {k: z["centre"][i] for k, z in Z.items()}
        print(f"{i:3d} {Z['D0.2']['t'][i]:5.0f} | {dist(r['D0.2'], r['D0.1']):10.3f} {dist(r['D0.1'], r['D0.05']):11.3f} | "
              f"{dist(f['D0.2'][0], f['D0.1'][0]):14.3f} {dist(f['D0.1'][0], f['D0.05'][0]):11.3f} | {f['D0.2'][1]:5d} | "
              f"{dist(f['D0.2'][0], r['D0.2']):.3f}")
print("--- registered centre tracks")
for k, z in Z.items():
    print(k, np.round(z["centre"], 3).tolist())
