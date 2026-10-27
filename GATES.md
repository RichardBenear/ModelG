# Gates

Each gate's prediction, pass criterion and null control are filled in and
committed **before** the run that tests it. Status: PENDING / PASS / FAIL / REVISED.

| Gate | Phase | Question | Status |
|---|---|---|---|
| MG1 | 1 | Continuum PDE reproduces a published stationary Model G soliton | PROPOSED (awaiting PI sign-off) |
| MG2 | 2 | Stochastic token rules converge to the MG1 PDE solution as Ω → ∞ | PENDING |
| MG3 | 2 | Particles nucleate spontaneously from intrinsic fluctuations at finite Ω | PENDING |
| MG4 | 3 | Particles survive irregular connectivity; graph Laplacian ≈ ∇² at soliton scale | PENDING |
| MG5 | 4 | Behavior of solitons on an evolving hypergraph (persist / drift / dissolve) | PENDING |
| MG6 | 5 | G-coupled growth produces a volume-growth deficit; falloff and mass scaling | PENDING |
| MG7 | 6 | Multiway branch divergence and foliation dependence | PENDING |

## Template

```
### MGn — <title>
Pre-registered: <date>, commit <hash>
Prediction:
Pass criterion (quantitative):
Fail criterion:
Null control(s):
Configs:
Result:
Replication (ChatGPT):
PI sign-off:
```

---

### MG1 — Continuum reproduction of the published Model G soliton
Status: PROPOSED 2026-10-27 (MG001). Becomes pre-registered on PI sign-off, before any MG1 run.

Reference: `spec/model_g_published.md` (Pulver & LaViolette 2013, eq. 17 parameters, Table 1).

Setup: dimensionless eqs. 7a–7c; eq. 17 parameters; seed eq. 18 added to the X equation;
φ = 0 initially; t ∈ [0, 100]. 1D Cartesian x ∈ [−50, 50]; 2D and 3D radial r ∈ [0, 50];
Dirichlet φ = 0 at the outer edge, Neumann at the center. Method-of-lines finite differences
with a stiff implicit integrator (our choice; the paper does not report its method).

Prediction: a stationary soliton forms in 1D, 2D and 3D, matching Table 1.

Pass criteria (all must hold):
1. Turing wavelength λ₀ = 3.08 ± 2% in 1D, 2D and 3D (measured as in the paper: twice the
   mean zero spacing).
2. Core amplitudes (Table 1 row a) for φ_G, φ_X, φ_Y within ±3% in 1D, 2D, 3D.
3. Core radii (row b) within ±3%.
4. Stationarity: max |∂φ/∂t| < 10⁻⁴ at t = 100.
5. Grid convergence: criteria 1–3 change by < 0.5% when the grid spacing is halved.

Fail criterion: any pass criterion fails at the converged resolution. Partial reproduction is
reported row by row, not rounded to a pass.

Null controls:
- N1: no seed → max |φ| < 10⁻⁶ at t = 100 (homogeneous state stable; consistent with the
  MG001 linear check, max Re σ = −0.195).
- N2 (exploratory, not pass/fail): seed amplitude scaled by 0.25, 0.5, 0.75 → report whether a
  soliton nucleates and estimate the nucleation threshold amplitude. This feeds MG3.

Out of scope for MG1 (later stages): the 2Dδ variant, multi-particle bonding (3.303 / 3.314),
gradient drift (v = −0.0599 m).

Configs: `configs/MG1_*.json` (to be written)
Result:
Replication (ChatGPT):
PI sign-off:
