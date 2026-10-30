# MG009 — ChatGPT design review of the MG4 proposal, draft 2 (verbatim, relayed by the PI 2026-10-30)

I support the 2D, noise-free scope, but would revise the criteria before registration. There are several concrete design issues—not merely residual risks.

1. C1’s confidence calculation is wrong, and its sampling unit needs clarification.
A two-sided 95% Clopper–Pearson lower bound ≥ 0.75 requires 30/32 successes, not 28/32. At 28/32 the lower bound is 0.710.
Also, eight placements on each of four graphs are clustered observations. They cannot generally be treated as 32 independent trials for a graph-family probability. An economical alternative is to retain the 4×8 design as a fixed benchmark panel, report results per graph, and use an explicit count threshold without a family-population confidence claim.
2. C3 measures less than its title claims.
Matching the eigenvalue count at one cutoff does not establish isotropy, plane-wave dispersion, localization, or the operator’s action on a soliton. Either narrow C3 to a spectral-count benchmark, or promote a finite-wavelength operator test to a criterion.
Finite-box effects matter: the exact square-torus count at λ = 3.08 is 137 modes, versus Weyl’s 132.47—already a 1.7% difference in inferred k. Use the exact finite-domain spectrum or T as the reference.
The plane-wave diagnostic must use allowed periodic wavevectors, not arbitrary directions at exactly k*; otherwise seam discontinuities can create artificial spectral width even on T.
3. The comparator geometry and calibration need precise definitions.
An unstrained equilateral triangular lattice cannot have an exactly square periodic supercell. Use a commensurate, nearly square rectangular torus for every family, or explicitly specify the strain and calibrate that geometry.
I support homogenized h, but report the full diffusivity tensor. The proposed `|Dxx−Dyy|/trace` misses off-diagonal anisotropy; use its eigenvalues. Also specify how random graphs attain the target h within a declared tolerance.
4. Detection and outcome classes need an implementation specification.
Define how scattered-node minima become distinct particles, including ties, merging/separation, empty radial bins and missing zero crossings. The current classes overlap: two particles at t = 200 followed by none at t = 300 satisfies both “multiplied” and “decayed.” Give an explicit precedence rule.
5. Survival is not stationarity.
The zero node-sum identities cannot be unconditional implementation checks at t = 300. First check the time-dependent balances:
\[ \dot I_G=-I_G+0.1I_X,\qquad \dot I_{\rm total}=-0.9I_X. \]
Require zero sums only after establishing stationarity, or retain the derivative terms. This is the same distinction that mattered in MG1.

For the proposal’s explicit questions: homogenized calibration and the continuum-profile survival initialization are appropriate. Seeded formation can remain descriptive. C2’s tolerances are reasonable benchmark choices once measurements are unambiguous. I would leave W optional/deferred; a true finite-volume version needs node volumes and weighted conservation, not just unspecified edge weights.
Finally, soften the failure interpretations: failed survival does not uniquely identify an existence-window shift, and failed profile matching does not mean a solution of graph Model G is “not Model G.”
The full review and arithmetic checks are in [MG009 design review](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG009_replication/README.md).
