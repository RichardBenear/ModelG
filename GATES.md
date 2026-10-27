# Gates

Each gate's prediction, pass criterion and null control are filled in and
committed **before** the run that tests it. Status: PENDING / PASS / FAIL / REVISED.

| Gate | Phase | Question | Status |
|---|---|---|---|
| MG1 | 1 | Continuum PDE reproduces a published stationary Model G soliton | PENDING |
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
