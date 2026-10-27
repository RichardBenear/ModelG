# Rewrite-Rule Specification v1.0 (FROZEN)

Status: frozen 2026-10-26 (MG000). Changes require `rules_v1.1.md` or later.

## 1. Source system: Model G (forward reactions)

| Reaction | Rate constant |
|---|---|
| A → G | k₁ |
| G → X | k₂ |
| B + X → Y + Z | k₃ |
| 2X + Y → 3X | k₄ |
| X → Ω | k₅ |

Mean-field reaction–diffusion form:

- ∂G/∂t = k₁A − k₂G + D_G ∇²G
- ∂X/∂t = k₂G + k₄X²Y − k₃BX − k₅X + D_X ∇²X
- ∂Y/∂t = k₃BX − k₄X²Y + D_Y ∇²Y

Reversible variants (e.g. G ⇌ X) and exact published parameter values are
pinned in Phase 0 (`spec/model_g_published.md`) and are **not** part of v1.0.
If a reverse step is adopted, it becomes rules v1.1.

## 2. Substrate

- **Space:** a hypergraph of binary relations {u, v}. Phase-dependent:
  fixed lattice (Phases 1–2), frozen irregular graph (Phase 3),
  evolving under a spatial rewrite rule (Phases 4–5).
- **Etherons:** unary label hyperedges {v, S}, S ∈ {G, X, Y}, attached to node v.
  n_S(v) denotes the count of species S at node v.
- **Reservoirs:** A and B are uniform and enter only through effective rates
  (k₁A, k₃B). Z and Ω are not represented (deleted on production).
- **Volume parameter:** Ω (tokens per node at unit concentration). Concentrations
  are n_S/Ω.

## 3. Rules

| ID | Rule | Propensity at node v |
|---|---|---|
| R1 | {v} → {v}, {v,G} | Ω · k₁A |
| R2 | {v,G} → {v,X} | k₂ · n_G |
| R3 | {v,X} → {v,Y} | k₃B · n_X |
| R4 | {v,X},{v,X},{v,Y} → {v,X},{v,X},{v,X} | (k₄/Ω²) · n_X(n_X−1) · n_Y |
| R5 | {v,X} → ∅ | k₅ · n_X |
| D_S | {v,S},{v,u} → {u,S},{v,u} | (D_S / h²) · n_S per incident edge |

h is the effective lattice spacing. On irregular graphs, D_S hopping uses the
graph Laplacian; its normalization is fixed per phase in that phase's config.

## 4. Updating scheme

- Pure Wolfram models have no rates. v1.0 selects one weighted path through the
  multiway system by continuous-time stochastic simulation (exact Gillespie for
  small systems, tau-leaping for large ones).
- This choice is explicitly a choice of foliation and is examined in Phase 6.

## 5. Expected mean-field limit

As Ω → ∞, the master equation of R1–R5 + D_S converges (Kurtz) to the Model G
equations above with ∇² replaced by the graph Laplacian. Gate MG2 tests this.

## 6. Reserved for Phase 5 (not part of v1.0)

Coupling of local n_G to the firing weight of the spatial rewrite rule.
Candidate coupling functions will be proposed in a later spec version and
chosen by the PI.
