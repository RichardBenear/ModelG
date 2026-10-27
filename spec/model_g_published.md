# Model G — Published Reference Specification (MG001)

Primary source for all numbers below:
**M. Pulver & P. A. LaViolette, "Stationary dissipative solitons of Model G,"
*International Journal of General Systems* 42(5), 519–541 (2013).**
doi:10.1080/03081079.2013.776313. Author-hosted PDF: https://etheric.com/Downloads/model_g.pdf

Transcribed 2026-10-27 in two independent read passes; values agreed in both.
Equation numbers refer to that paper. Anything that is our inference rather than
the paper's statement is marked **[ours]**.

## 1. Reaction scheme (eqs. 3a–3e)

| | Reaction | Forward | Reverse |
|---|---|---|---|
| 3a | A ⇌ G | k₁ | k₋₁ |
| 3b | G ⇌ X | k₂ | k₋₂ |
| 3c | B + X ⇌ Y + Z | k₃ | k₋₃ |
| 3d | 2X + Y ⇌ 3X | k₄ | k₋₄ |
| 3e | X ⇌ Ω | k₅ | k₋₅ |

In the published soliton runs, **only k₋₂ is nonzero among the reverse rates**
(g = 1/10; s = u = w = 0, and p = q implies k₋₁ = 0). The G ⇌ X back-reaction is
the feedback that couples X back into G and creates the G-well.

## 2. Dimensional PDEs (eqs. 4a–4c)

- ∂G/∂t = D_G∇²G − (k₋₁ + k₂)G + k₋₂X + k₁A
- ∂X/∂t = D_X∇²X + k₂G − (k₋₂ + k₃B + k₅)X + k₋₃ZY − k₋₄X³ + k₄X²Y + k₋₅Ω
- ∂Y/∂t = D_Y∇²Y + k₃BX − k₋₃ZY + k₋₄X³ − k₄X²Y

## 3. Nondimensionalization

With K = k₋₂ + k₅:

- Scales: T = 1/K, L = √(D_G T), C = 1/√(k₄ T)
- d_x = D_X/D_G, d_y = D_Y/D_G
- a = k₁A√k₄ / K^{3/2}, b = k₃B/K, g = k₋₂/K, p = k₂/K, q = (k₋₁ + k₂)/K,
  s = k₋₄/k₄, u = k₋₃Z/K, w = √(k₄k₋₅) Ω / K^{3/2}

## 4. Dimensionless PDEs (eqs. 7a–7c)

- ∂G/∂t = ∇²G − qG + gX + a
- ∂X/∂t = d_x∇²X + pG − (1 + b)X + uY − sX³ + X²Y + w
- ∂Y/∂t = d_y∇²Y + bX − uY + sX³ − X²Y

Potentials: φ_G = G − G₀, φ_X = X − X₀, φ_Y = Y − Y₀ (eqs. 13a–c).

**Homogeneous steady state:** G₀ = (a + gw)/(q − gp), X₀ = (pa + qw)/(q − gp),
Y₀ = X₀(sX₀² + b)/(X₀² + u).

## 5. Reference parameter set (eq. 17)

| d_x | d_y | a | b | g | p | q | s | u | w |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 12 | 14 | 29 | 1/10 | 1 | 1 | 0 | 0 | 0 |

Steady state: G₀ = X₀ = 140/9 ≈ 15.5556, Y₀ = b/X₀ ≈ 1.86429.

## 6. Numerical setup (as published)

- 1D: x ∈ [−50, 50]; 2D/3D: radial r ∈ [0, 50] with circular/spherical symmetry
  imposed, ∇²φ = φ_rr + ((n−1)/r)φ_r.
- t ∈ [0, 100]. Initial condition φ = 0 everywhere.
- Boundaries: Dirichlet φ = 0 at the outer edge; Neumann ∂φ/∂r = 0 at r = ε, ε = 10⁻⁹.
- Positivity: φ_S ≥ −S₀.
- **Not reported:** software, grid spacing, time step, integration scheme.
  MG1 must therefore establish its own convergence.

**Seed (eq. 18), added to the RHS of the X equation:**
χ(r, t) = −exp(−r²/2) · exp(−(t − 10)²/18)

Multi-particle seeds: χ₂ = χ(x + d₂/2) + χ(x − d₂/2); χ₃ = χ(x + d₃) + χ(x) + χ(x − d₃).

Core-broadening variant "2Dδ": ∇² → δ(r)∇² with δ(r) = 1 + exp(−4r²/9).

## 7. Published results (Table 1 and text)

| Row | Field | 1D | 2D | 3D | 2Dδ |
|---|---|---|---|---|---|
| (a) core amplitude | φ_Y | 0.930 | 1.50 | 1.70 | 1.53 |
| | φ_G | −0.161 | −0.308 | −0.411 | −0.320 |
| | φ_X | −8.36 | −13.6 | −14.6 | −13.8 |
| (b) core radius (1st zero) | φ_Y | 0.724 | 1.21 | 1.68 | 1.67 |
| | φ_G | 0.867 | 1.37 | 1.85 | 1.87 |
| | φ_X | 0.642 | 1.13 | 1.61 | 1.58 |
| (c) core RMS radius | φ_Y | 0.302 | 0.665 | 1.07 | 0.929 |
| | φ_G | 0.363 | 0.739 | 1.14 | 1.02 |
| | φ_X | 0.272 | 0.634 | 1.04 | 0.892 |
| (d) core integral | φ_Y | 0.806 | 3.12 | 13.2 | 6.31 |
| | φ_G | −0.167 | −0.710 | −3.17 | −1.44 |
| | φ_X | −6.59 | −26.3 | −111 | −53.6 |
| (e) full integral | φ_Y | 0.281 | 0.893 | 3.30 | 2.96 |
| | φ_G | 8.39e−5 | 2.88e−3 | 2.87e−2 | −0.480 |
| | φ_X | 6.01e−4 | 0.128 | 1.93 | −13.7 |
| (f) Turing wavelength λ₀ | all | 3.08 | 3.08 | 3.08 | 3.08 |

λ₀ was measured as twice the mean spacing of successive zeros (centrally aligned).

Other reported behavior:
- No seed → the system stays at the homogeneous steady state.
- Two seeds settle at separation 3.303 (also from an initial d₂ = 4); three seeds settle at 3.314;
  a seed pair at d₂ = 6.628 spontaneously creates a third particle.
- In an imposed G-gradient (switched on at t = 100), particle velocity tends to a constant
  v = −0.0599 m (m = gradient slope), linear over 11 slopes. Acceleration was not shown.
- 3D tail amplitude decays roughly as r^−3.7 (≈2λ₀), r^−7 (≈4λ₀), r^−10 (≈6λ₀).

## 8. Independent checks [ours] — `src/linear_stability.py`, `results/MG001/`

- Steady-state residuals vanish at the eq. 17 values (transcription self-consistent).
- The homogeneous state is **linearly stable**: max Re σ(k) = −0.195 at k* = 1.995
  (λ* = 3.149). The system is **subcritical**, matching the paper's claim that it needs a seed.
- Turing threshold: b_c = 30.044 (λ_c = 2.993). The published b = 29 is **3.5% below threshold**.
- The least-stable linear wavelength (3.15) is within 2.3% of the reported nonlinear λ₀ = 3.08.

## 9. Implications for this project [ours]

1. **Rule spec needs the G ⇌ X back-reaction.** Rules v1.0 omitted k₋₂, but the
   reference set uses g = 0.1. This is fixed in `spec/rules_v1.1.md`.
2. **The "G-well" is a core feature, not a net deficit.** φ_G is negative in the core,
   but its full-space integral is slightly **positive** in standard 1D/2D/3D. It is negative
   (the authors' "positive gravitational mass") only in the core-broadened 2Dδ variant.
   MG6's "integrated G-deficit" must therefore be defined on the core (row d), or use the δ variant.
   This decision is deferred to the Phase 5 spec.
3. **The X core is nearly depleted.** In 3D, X ≈ X₀ + φ_X ≈ 15.56 − 14.6 ≈ 0.96 at the center,
   about 6% of background. Relative stochastic fluctuations will be largest exactly there,
   which matters for choosing Ω in MG2/MG3.
4. **Resolution requirement.** λ₀ ≈ 3.08 L. Lattices must resolve this, e.g. spacing h ≲ 0.3
   (≥ 10 sites per wavelength), which is what makes Phase 2–5 graphs large.
5. **Proximity to threshold.** At 3.5% below b_c, finite-Ω noise may nucleate particles
   (MG3) or drive spontaneous patterning. Distance to b_c is the natural control variable for MG3.

## 10. Secondary sources (not used for numbers)

- P. A. LaViolette, "An introduction to subquantum kinetics," *Int. J. General Systems* 11 (1985):
  the original Model G papers.
- P. A. LaViolette, "The electric charge and magnetization distribution of the nucleon,"
  *Physics Essays* 21 (2008): https://starburstfound.org/downloads/physics/nucleon.pdf.
  Gives only the forward scheme and no numerical parameters.
- Model G Vortical Motion Group, "Explanation of Model G in subquantum kinetics" (2020):
  https://starburstfound.org/downloads/physics/explanation-Model-G.pdf. Same reversible scheme;
  states that only the G ⇌ X reverse rate is set large.
