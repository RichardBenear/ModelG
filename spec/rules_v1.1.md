# Rewrite-Rule Specification v1.1 (FROZEN)

Status: frozen 2026-10-27 (MG001). Supersedes v1.0 for all runs.
Change from v1.0: **adds R2r (X → G back-reaction)**, required by the published
parameter set (g = 1/10). Rates are now given in the dimensionless units of
`spec/model_g_published.md`, so eq. 17 parameters plug in directly.

## 1. Substrate (unchanged from v1.0)

- Space: hypergraph of binary relations {u, v}, phase-dependent.
- Etherons: unary hyperedges {v, S}, S ∈ {G, X, Y}; n_S(v) = count at node v.
- Reservoirs A, B uniform and implicit; Z, Ω not represented.
- Ω: tokens per node at unit dimensionless concentration (concentration = n_S/Ω).
- Time unit T = 1/(k₋₂ + k₅); length unit L = √(D_G T).

## 2. Rules and propensities (dimensionless)

| ID | Rule | Reaction | Propensity at node v |
|---|---|---|---|
| R1 | {v} → {v}, {v,G} | A → G | Ω · a |
| R2 | {v,G} → {v,X} | G → X | p · n_G |
| **R2r** | **{v,X} → {v,G}** | **X → G** | **g · n_X** |
| R3 | {v,X} → {v,Y} | B + X → Y + Z | b · n_X |
| R4 | {v,X},{v,X},{v,Y} → {v,X},{v,X},{v,X} | 2X + Y → 3X | n_X (n_X − 1) n_Y / Ω² |
| R5 | {v,X} → ∅ | X → Ω | (1 − g) · n_X |
| R1r | {v,G} → ∅ | G → A | (q − p) · n_G |
| D_S | {v,S},{v,u} → {u,S},{v,u} | diffusion | (d_S / h²) · n_S per incident edge |

d_G = 1, d_X = d_x, d_Y = d_y. h = effective lattice spacing (in units of L).

R1r is listed for completeness. For the reference set (p = q) its rate is zero.
Reverse steps 3c, 3d, 3e (u, s, w) are **not** implemented in v1.1. They are zero
in the reference set.

R5's coefficient follows from the scaling: X's total first-order loss is
(k₋₂ + k₃B + k₅)/K = g + b + (1 − g) = 1 + b.

## 3. Mean-field limit (check)

Summing propensity × stoichiometry and dividing by Ω reproduces eqs. 7a–7c
(with s = u = w = 0) exactly, with ∇² → graph Laplacian:

- dG/dt = a − qG + gX
- dX/dt = pG − (1 + b)X + X²Y   (R4's n_X(n_X−1) → X² as Ω → ∞)
- dY/dt = bX − X²Y

## 4. Updating scheme

Continuous-time stochastic simulation selects one weighted path through the
multiway system (exact Gillespie for small systems, tau-leaping for large ones).
This is a choice of foliation and is examined in Phase 6.

## 5. Reserved for Phase 5

Coupling of local n_G to the spatial rewrite rate.

Erratum (MG001a, non-normative note only; the rules above are unchanged): the earlier note
here is superseded. With conservative hopping (D_S as specified), a stationary localized
soliton has ∫φ_G = ∫φ_X = 0 in infinite space, given integrable perturbations,
vanishing integrated time derivatives and vanishing limiting fluxes, whenever q − gp ≠ 0.
See `model_g_published.md` §9.2.
How "mass" is defined for Phase 5 is an open PI decision.
