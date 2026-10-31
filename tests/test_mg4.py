"""Tests for MG4 Part 0 code (src/graphs.py). Run: python tests/test_mg4.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import graphs as gr


def test_torus_geometry():
    t = gr.torus(0.2, 122, 142)
    assert abs(t["Lx"] - 19.923) < 1e-3 and abs(t["Ly"] - 20.082) < 1e-3
    t3 = gr.torus(0.3, 82, 94)
    assert abs(t3["Lx"] - 20.086) < 1e-3 and abs(t3["Ly"] - 19.940) < 1e-3


def test_triangular_structure_and_corrector():
    t = gr.torus(0.2, 20, 22)
    g = gr.triangular(t)
    assert g["N"] == 440 and len(g["I"]) == 3 * 440
    assert np.all(gr.degrees(g) == 6)
    assert np.allclose(np.linalg.norm(g["dx"], axis=1), t["a"], atol=1e-12)
    L = gr.laplacian(g)
    assert abs(L - L.T).max() < 1e-15 and np.abs(L.sum(1)).max() < 1e-12
    c = gr.corrector(g)
    a2 = t["a"] ** 2
    assert abs(c["h2"] - 1.5 * a2) < 1e-10 and abs(c["eig"][0] - c["eig"][1]) < 1e-10
    assert abs(c["h2_naive"] - 1.5 * a2) < 1e-10
    assert abs(c["h2"] - 0.04) < 1e-12        # h_t = 0.2 exactly


def test_triangular_spectrum_matches_analytic():
    t = gr.torus(0.2, 24, 28)
    g = gr.triangular(t)
    h = 0.2
    nu = gr.spectrum(g, h, cutoff=30.0)
    # analytic: nu(q) = (1/h^2) sum_{3 dirs} 2 (1 - cos(q.d)) over allowed q of the lattice's own periodicity
    a = t["a"]; d = a * np.array([[1, 0], [0.5, gr.SQ3 / 2], [-0.5, gr.SQ3 / 2]])
    # reciprocal lattice of the supercell: q = 2 pi (m / Lx, n / Ly), with the lattice's Brillouin zone counted once:
    # build the operator explicitly on Fourier modes by evaluating the Rayleigh quotient of each node-sampled plane wave
    vals = []
    for m in range(-t["nx"] // 2, t["nx"] // 2):
        for n in range(-t["ny"] // 2, t["ny"] // 2):
            q = 2 * np.pi * np.array([m / t["Lx"], n / t["Ly"]])
            vals.append((2 * (1 - np.cos(d @ q))).sum() / h**2)
    vals = np.sort(np.array(vals))
    ref = vals[vals <= 30.0]
    # the node-sampled modes double-count nothing only if the set above is the lattice's full Fourier basis;
    # compare the low part, where aliasing is impossible
    n = min(len(ref), len(nu), 40)
    assert np.allclose(nu[:n], ref[:n], atol=1e-8), (nu[:8], ref[:8])


def test_delaunay_is_torus_triangulation():
    rng = np.random.default_rng(5)
    Lx, Ly = 10.0, 9.0
    pos = gr.poisson_points(600, Lx, Ly, rng)
    g = gr.delaunay(pos, Lx, Ly)
    assert len(g["I"]) == 3 * g["N"]        # Euler: V - E + F = 0 and 3F = 2E on the torus
    assert abs(gr.degrees(g).mean() - 6.0) < 1e-12
    assert np.linalg.norm(g["dx"], axis=1).max() < min(Lx, Ly) / 2
    L = gr.laplacian(g)
    assert np.abs(L.sum(1)).max() < 1e-12
    c = gr.corrector(g)
    assert c["eig"][1] > 0 and c["h2"] <= c["h2_naive"] + 1e-12   # homogenized <= Voigt bound


def test_rgg_component_and_regular():
    rng = np.random.default_rng(6)
    pos = gr.poisson_points(800, 8.0, 8.0, rng)
    g = gr.rgg(pos, 8.0, 8.0, radius=np.sqrt(12 / (np.pi * 800 / 64)))
    assert 0.9 < g["kept_fraction"] <= 1.0
    assert np.linalg.norm(g["dx"], axis=1).max() <= g["radius"] + 1e-12
    from scipy.sparse.csgraph import connected_components
    assert connected_components(gr.laplacian(g) != 0, directed=False)[0] == 1
    r = gr.random_regular(200, 6, seed=3)
    assert np.all(gr.degrees(r) == 6)


def test_torus_count():
    kstar = 2 * np.pi / 3.08
    assert gr.torus_count(kstar, 20.0, 20.0) == 137
    t = gr.torus(0.2, 122, 142)
    assert gr.torus_count(kstar, t["Lx"], t["Ly"]) == 137
    assert gr.torus_count(0.0, 20.0, 20.0) == 1


def test_balances_exact_on_random_state():
    import graph_ode as go
    rng = np.random.default_rng(7)
    pos = gr.poisson_points(500, 6.0, 6.0, rng)
    g = gr.delaunay(pos, 6.0, 6.0)
    y = rng.normal(0, 1, 3 * g["N"])
    rhs, _, _ = go.make_system(g, 0.2)
    b = go.balances(y, rhs(0.0, y), g["N"])
    assert b["ok"], b
    z = np.zeros(3 * g["N"])                       # exact homogeneous state: M = 0, residual 0, still defined
    bz = go.balances(z, rhs(0.0, z), g["N"])
    assert bz["ok"] and np.abs(rhs(0.0, z)).max() < 1e-12


def test_jacobian_matches_finite_difference():
    import graph_ode as go
    rng = np.random.default_rng(8)
    g = gr.delaunay(gr.poisson_points(60, 3.0, 3.0, rng), 3.0, 3.0)
    rhs, jac, _ = go.make_system(g, 0.2)
    y = rng.normal(0, 0.5, 3 * g["N"])
    J = jac(0.0, y).toarray()
    eps = 1e-6
    Jn = np.column_stack([(rhs(0, y + eps * e) - rhs(0, y - eps * e)) / (2 * eps) for e in np.eye(len(y))])
    assert np.abs(J - Jn).max() < 1e-5 * max(1.0, np.abs(J).max())


def test_nd_order_matches_default_solver():
    """MG010 speed-up: the nested-dissection path gives the same steps and the same solution to rounding."""
    import graph_ode as go
    rng = np.random.default_rng(9)
    g = gr.delaunay(gr.poisson_points(900, 6.0, 6.0, rng), 6.0, 6.0)
    o = gr.nd_order(g)
    assert sorted(o) == list(range(g["N"]))
    y0 = go.soliton_ic(g, np.array([3.0, 3.0]))
    a, _ = go.integrate(g, 0.2, y0, 2.0, order=None)
    b, _ = go.integrate(g, 0.2, y0, 2.0, order="nd")
    assert a.nfev == b.nfev and a.nlu == b.nlu
    assert np.abs(a.y[:, -1] - b.y[:, -1]).max() < 1e-10


def test_cotan_operator_construction():
    """MG4c Part 3 construction checks (proposal §4) on a small periodic Delaunay graph and on T."""
    import graph_ode as go
    rng = np.random.default_rng(11)
    Lx, Ly = 6.0, 5.0
    g = gr.delaunay(gr.poisson_points(700, Lx, Ly, rng), Lx, Ly)
    c = gr.cotan(g)
    assert len(c["faces"]) == 2 * g["N"]                              # Euler: F = 2V on the torus
    assert abs(c["face_area"].sum() / (Lx * Ly) - 1) < 1e-12
    assert abs(c["A"].sum() / (Lx * Ly) - 1) < 1e-12 and c["A"].min() > 0
    ge = set(map(tuple, np.column_stack([g["I"], g["J"]])))
    fe, cnt = np.unique(c["edges"], axis=0, return_counts=True)
    assert set(map(tuple, fe)) == ge and np.all(cnt == 2)          # every graph edge in exactly two kept faces
    assert c["W"].data.min() >= -1e-12                               # Delaunay: nonnegative cotangent weights
    K = c["K"]
    assert abs(K - K.T).max() < 1e-12
    assert np.abs(c["L"] @ np.ones(g["N"])).max() < 1e-10
    f = rng.normal(size=g["N"])
    Lf = c["L"] @ f
    assert abs((c["A"] * Lf).sum()) < 1e-12 * (c["A"] * np.abs(Lf)).sum()
    # on T the cotangent operator equals MG4's L_G = (A - D)/h^2 with h^2 = 1.5 a^2
    tor = gr.torus(0.2, 30, 34)
    gT = gr.triangular(tor)
    cT = gr.cotan(gT)
    assert len(cT["faces"]) == 2 * gT["N"]
    LT = gr.laplacian(gT) / (1.5 * tor["a"] ** 2)
    assert abs(cT["L"] - LT).max() < 1e-12 * abs(LT).max()
    # Jacobian with the operator against finite differences, and area-weighted balances
    gs = gr.delaunay(gr.poisson_points(60, 3.0, 3.0, rng), 3.0, 3.0)
    cs = gr.cotan(gs)
    rhs, jac, _ = go.make_system(gs, None, op=cs["L"])
    y = rng.normal(0, 0.5, 3 * gs["N"])
    J = jac(0.0, y).toarray()
    eps = 1e-6
    Jn = np.column_stack([(rhs(0, y + eps * e) - rhs(0, y - eps * e)) / (2 * eps) for e in np.eye(len(y))])
    assert np.abs(J - Jn).max() < 1e-5 * max(1.0, np.abs(J).max())
    assert go.balances(y, rhs(0, y), gs["N"], w=cs["A"])["ok"]


def test_measurement_on_synthetic_fields():
    import measure_graph as mg
    g = gr.triangular(gr.torus(0.2, 62, 70))
    S = mg.smoother(g)
    assert abs(S.sum(1).max() - 1) < 1e-12
    well = lambda c: -6.0 * np.exp(-np.hypot(*gr.min_image(g["pos"] - c, g["Lx"], g["Ly"]).T) ** 2 / 0.5)
    c1 = np.array([3.0, 3.0]); c2 = np.array([7.5, 6.0]); c3 = np.array([3.6, 3.0])
    one = mg.particles(g, S @ well(c1), -2.0)
    assert len(one) == 1 and np.hypot(*gr.min_image(one[0]["centre"] - c1, g["Lx"], g["Ly"])) < 0.05
    assert len(mg.particles(g, S @ (well(c1) + well(c2)), -2.0)) == 2
    assert len(mg.particles(g, S @ (well(c1) + well(c3)), -2.0)) == 1      # 0.6 apart: merged
    assert mg.particles(g, S @ np.zeros(g["N"]), -2.0) == []
    P = lambda c: [dict(centre=np.asarray(c, float))]
    assert mg.classify(g, [P(c1)] * 11)[0] == "single"
    assert mg.classify(g, [P(c1) + P(c2)] + [[]] * 10)[0] == "multiplied"   # precedence over decayed
    assert mg.classify(g, [P(c1)] * 10 + [[]])[0] == "decayed"
    assert mg.classify(g, [P(c1)] * 5 + [[]] + [P(c1)] * 5)[0] == "unresolved"
    assert mg.classify(g, [P(c1)] * 5 + [P(c2)] * 6)[0] == "unresolved"      # broken link
    rc, prof, cnt = mg.radial_profile(g, np.cos(2 * np.pi * np.hypot(*gr.min_image(g["pos"] - c1, g["Lx"], g["Ly"]).T) / 3.0), c1)
    assert abs(mg.wavelength(rc, prof) - 3.0) < 0.05
    assert mg.wavelength(rc, np.ones_like(prof)) is None
    assert mg.profile_l2(rc, prof, prof) == 0.0


if __name__ == "__main__":
    import time
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            t0 = time.time(); f(); print(f"{name}: ok ({time.time() - t0:.1f}s)")
