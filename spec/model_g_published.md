# Model G — Published Reference Specification (MG001)

Primary source for all numbers below:
**M. Pulver & P. A. LaViolette, "Stationary dissipative solitons of Model G,"
*International Journal of General Systems* 42(5), 519–541 (2013).**
doi:10.1080/03081079.2013.776313. Author-hosted PDF: https://etheric.com/Downloads/model_g.pdf

Transcribed 2026-10-27 in two independent read passes; values agreed in both.
Revised MG001a (2026-10-27) after ChatGPT review: `review/MG001_chatgpt.md`, `review/MG001_response.md`.
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

In the published reference PDEs, **only the X → G reverse contribution remains active**
(g = 1/10). The other reverse contributions vanish as *combinations*: p = q gives k₋₁ = 0;
s = 0 gives k₋₄ = 0; u = 0 gives k₋₃Z = 0; w = 0 gives k₋₅Ω = 0. The last two do not fix
k₋₃ or k₋₅ individually without assuming nonzero reservoir concentrations Z and Ω.
(Corrected in MG001a per review item 2.) The G ⇌ X back-reaction couples X back into G.

## 2. Dimensional PDEs (eqs. 4a–4c)

- ∂G/∂t = D_G∇²G − (k₋₁ + k₂)G + k₋₂X + k₁A
- ∂X/∂t = D_X∇²X + k₂G − (k₋₂ + k₃B + k₅)X + k₋₃ZY − k₋₄X³ + k₄X²Y + k₋₅Ω
- ∂Y/∂t = D_Y∇²Y + k₃BX − k₋₃ZY + k₋₄X³ − k₄X²Y

## 3. Nondimensionalization

With K = k₋₂ + k₅:

- Scales: T = 1/K, L = √(D_G T), C = 1/√(k₄ T)
- d_x = D_X/D_G, d_y = D_Y/D_G
- a = k₁A√k₄ / K^{3/2}, b = k₃B/K, g = k₋₂/K, p = k₂/K, q = (k₋₁ + k₂)/K,
  s = k₋₄/k₄, u = k₋₃Z/K, **w = √k₄ · k₋₅Ω / K^{3/2}**

Check: a source term σ in a dimensional equation scales to (T/C)σ, and T/C = √k₄ / K^{3/2}.
This gives a = √k₄ k₁A / K^{3/2} and w = √k₄ k₋₅Ω / K^{3/2}. (MG001a: the earlier
√(k₄k₋₅) was a transcription error and is dimensionally inconsistent. It has no effect on the
reference set, where w = 0.)

## 4. Dimensionless PDEs (eqs. 7a–7c)

- ∂G/∂t = ∇²G − qG + gX + a
- ∂X/∂t = d_x∇²X + pG − (1 + b)X + uY − sX³ + X²Y + w
- ∂Y/∂t = d_y∇²Y + bX − uY + sX³ − X²Y

Potentials (eq. 12): φ_G = G − G₀, φ_X = X − X₀, φ_Y = Y − Y₀. Potential-form PDEs: eqs. 13a–c.

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

Core-broadening variant "2Dδ" (eqs. 24–25): ∇² → δ(r)∇² with δ(r) = 1 + exp(−4r²/9),
applied to all three Laplacian terms (the paper states one blanket substitution into eqs. 19).
**Note:** δ∇²φ is *not* conservative diffusion. It differs from ∇·(δ∇φ) by ∇δ·∇φ, so it does
not conserve ∫φ under transport. This matters for graph translation (§9) and for the integral
identity in §9.2. A graph hopping rule with a node-dependent rate δ(source) gives yet another
form, ∇²(δφ). The three forms must not be conflated.

### 6a. Benchmark measurement definitions (§3.2, eqs. 20–21)

- **Core radius r₀:** each field's own innermost zero crossing. φ_G, φ_X, φ_Y each have their own r₀.
- **Core integral (eq. 21):** 1D ∫₋ᵣ₀^{r₀} φ dr; 2D 2π∫₀^{r₀} φ r dr; 3D 4π∫₀^{r₀} φ r² dr.
- **Core RMS radius (eq. 20):** √(∫φ r² dS / ∫φ dS) over the core, with the same geometric measure.
- **Full integral:** the same as eq. 21 with r₀ → 50, the domain edge. It is a **finite-domain** quantity.
- **Turing wavelength:** twice the mean spacing of successive zeros, using centrally aligned points.

### 6b. Gradient benchmark (§4.3, eqs. 26–27): separate setup

This benchmark does **not** use the zero-initial-condition setup above.
- 1D, x ∈ [−50, 50], t ∈ [0, 200].
- Gradient γ_m(x,t) = 0 for t < 100, and m·x + 10⁻⁴/6 for t ≥ 100 (eq. 26).
  It is added as a source to the φ_G equation (13a) **and** sets the G boundary:
  φ_G(±50,t) = γ_m(±50,t) (eq. 27e). φ_X and φ_Y keep zero boundaries (27f–g).
- Gaussian initial conditions (27b–d), α = 3875/4096:
  φ_G = −0.161α·e^{−½(x/0.363)²}, φ_X = −8.37α·e^{−½(x/0.272)²}, φ_Y = 0.930α·e^{−½(x/0.302)²}.
  (The paper uses −8.37 here and −8.36 in Table 1.)
- Slopes m = −10⁻⁵k/3, k = 0…10.
- Position: midpoint of the two central zeros of φ_X flanking its minimum. Velocity: averaged
  over t = 150–200 (the differencing method is not stated).
- Result: v = −0.0599 m.

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
  (λ* = 3.149). The operating point is therefore below the linear Turing threshold.
  This alone does **not** establish a subcritical bifurcation, which needs nonlinear branch
  analysis. If MG1 confirms that a finite-amplitude stationary soliton coexists with this stable
  homogeneous state, that shows bistability. Bistability is consistent with, but not proof of, a
  subcritical branch. (Wording tightened in MG001a per review item 6. The numbers were
  independently reproduced by ChatGPT.)
- Turing threshold: b_c = 30.044 (λ_c = 2.993). The published b = 29 is **3.5% below threshold**.
- The least-stable linear wavelength (3.15) is within 2.3% of the reported nonlinear λ₀ = 3.08.

## 9. Implications for this project [ours]

1. **Rule spec needs the G ⇌ X back-reaction.** Rules v1.0 omitted k₋₂, but the
   reference set uses g = 0.1. This is fixed in `spec/rules_v1.1.md`.
2. **Integral identity (MG001a, wording per review round 2).** Integrate eqs. 13a–c over a
   domain V with constant-coefficient (conservative) diffusion. The X²Y, b, s and u terms cancel in
   the sum of all three equations, leaving two linear balances (F_S = d_S∮∇φ_S·n dS):
   - G equation: dI_G/dt = F_G − q I_G + g I_X
   - Sum: d(I_G + I_X + I_Y)/dt = F_G + F_X + F_Y − (q − p) I_G − (1 − g) I_X

   **Assumptions:** (i) the perturbations φ_S are integrable over ℝⁿ; (ii) the integrated time
   derivatives vanish (stationarity); (iii) the boundary fluxes F_S → 0 in the limit V → ℝⁿ.
   Under these, the balances become

   ⎡ q      −g  ⎤ ⎡I_G⎤
   ⎣ q − p  1 − g⎦ ⎣I_X⎦ = 0,  with determinant q − gp.

   So **I_G = I_X = 0 whenever q − gp ≠ 0** (reference set: q − gp = 0.9). The balances impose
   **no condition on I_Y**. (Determinant form and assumption list from ChatGPT review round 2.)

   Consequences, stated conservatively:
   - **Table 1 row (e).** If the published R = 50 states were stationary, the balances would
     *require* boundary fluxes F_G = −0.1643 and F_total = +1.737 in 3D
     (`results/MG001/integral_identity.txt`). These are inferred required fluxes, not measured
     ones. The nonzero full integrals are compatible with finite-boundary flux, but nonstationarity,
     numerical error or inconsistent tabulation remain alternative explanations. MG1 will
     distinguish them by measuring the fluxes.
   - **2Dδ.** The negative I_G = −0.480 (the authors' "positive gravitational mass") is not
     excluded by the unweighted zero-integral identity, because δ∇² changes that balance. Finite
     boundaries and incomplete stationarity still need checking before it is attributed to δ.
   - **Phase 5 (conditional constraint, not a conclusion).** Under the assumptions above, the
     integrated G perturbation of a standard Model G soliton is zero. Whether that rules out a 1/r
     geometric response depends on the coupling and field equation still to be specified. For
     example, a linear Poisson-type equation sourced proportionally to φ_G would have no monopole
     term, but a nonlinear growth coupling need not inherit the zero integral. Any Phase 5
     coupling model must state how it interacts with this identity. Defining "mass" from the
     core integral would be a new project convention, not the paper's, and is a PI decision.
   - MG1 tests the identity numerically, including the time-derivative terms (see `GATES.md`).

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
- P. A. LaViolette, "The electric charge and magnetization distribution of the nucleon:
  Evidence of a subatomic Turing wave pattern," *Int. J. General Systems* 37(6) (2008),
  649–676: https://starburstfound.org/downloads/physics/nucleon.pdf. (Citation corrected in MG001a.
  Verify volume and pages against the journal record before the paper draft.)
  Gives only the forward scheme and no numerical parameters.
- Model G Vortical Motion Group, "Explanation of Model G in subquantum kinetics" (2020):
  https://starburstfound.org/downloads/physics/explanation-Model-G.pdf. Same reversible scheme;
  states that only the G ⇌ X reverse rate is set large.
