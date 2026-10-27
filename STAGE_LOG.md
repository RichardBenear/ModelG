# Stage Log

| Stage | Date | Phase | Summary | Gate | Commit |
|---|---|---|---|---|---|
| MG000 | 2026-10-26 | 0 | Repo scaffold; charter (README), rules spec v1.0 frozen, gate table created | — | 1bc0015 |
| MG001 | 2026-10-27 | 0 | Published Model G pinned (Pulver & LaViolette 2013); linear check confirms subcritical, 3.5% below Turing threshold; rules v1.1 adds G⇌X back-reaction; MG1 criteria proposed | MG1 proposed | (this commit) |

## MG000 — Project scaffold
- Created charter, `spec/rules_v1.0.md`, `GATES.md`, folder layout.
- Next (MG001): pin LaViolette's published Model G equations, reverse-step variants and parameter sets from primary sources into `spec/model_g_published.md`.


## MG001 — Published Model G pinned
- Source: Pulver & LaViolette, Int. J. Gen. Syst. 42 (2013) 519–541. Transcribed in two independent passes → `spec/model_g_published.md`.
- Reference set (eq. 17): d_x=1, d_y=12, a=14, b=29, g=1/10, p=q=1, s=u=w=0. Only reverse step used: G ⇌ X.
- Independent linear check (`src/linear_stability.py`, `results/MG001/linear_stability.txt`): steady state G₀=X₀=15.556, Y₀=1.864; homogeneous state linearly stable (max Re σ = −0.195); Turing threshold b_c = 30.044, so b = 29 is 3.5% subcritical; least-stable λ = 3.149 vs published nonlinear λ₀ = 3.08 (2.3% apart).
- Findings affecting the plan:
  1. Rules v1.0 lacked the G ⇌ X back-reaction → `spec/rules_v1.1.md` (v1.0 superseded).
  2. Full-space φ_G integral is slightly positive in standard Model G; the G-well is a core feature. MG6's "G-deficit" must be defined on the core or use the 2Dδ variant (PI decision at Phase 5).
  3. X core is ~94% depleted in 3D → largest relative noise at the core; relevant to Ω choice in MG2/MG3.
  4. Paper gives no grid/step/method → MG1 includes an explicit grid-convergence criterion.
- MG1 pass/fail criteria proposed in `GATES.md`, awaiting PI sign-off.
- Next (MG002): implement the continuum solver and run MG1 after sign-off.
