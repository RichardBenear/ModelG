# ChatGPT review of MG001 (round 1) — verbatim
Received 2026-10-27 via PI. Reviewed file: spec/model_g_published.md @ 110bf8b.

The specification is mostly faithful, but I would request changes before using it as the implementation reference. I checked [the current file](https://github.com/RichardBenear/ModelG/blob/main/spec/model_g_published.md) against the [cited PDF](https://etheric.com/Downloads/model_g.pdf), including the rendered equations and Table 1.

1. Incorrect nondimensionalization of \(w\) — section 3, line 39.
The square root should cover only \(k_4\):
\[ \boxed{w=\frac{\sqrt{k_4}\,k_{-5}\Omega}{K^{3/2}}} \]
This follows independently by scaling the dimensional source: \(w=(T/C)k_{-5}\Omega\). The current \(\sqrt{k_4k_{-5}}\) is dimensionally inconsistent. It has no effect on the reference run with \(w=0\), but would break general implementations. See paper eqs. (6), (8), (9).
2. Reverse-rate claim is stronger than the parameters establish — section 1.
From \(u=0\), we know \(k_{-3}Z=0\); from \(w=0\), we know \(k_{-5}\Omega=0\). Neither establishes the individual rate constant without assuming nonzero reservoir concentrations. Replace the claim with “only the \(X\to G\) reverse contribution remains active in the reference PDEs.”
3. Gradient experiment needs its own numerical setup — sections 6–7.
Paper eqs. (26)–(27) use Gaussian initial conditions, altered G-boundaries, and forcing through \(t=200\). Applying the general zero-initial-condition setup will not reproduce this benchmark. The velocity coefficient is transcribed correctly. [Paper §4.3](https://etheric.com/Downloads/model_g.pdf#page=18).
4. The mass interpretation needs qualification — section 9.2.
Table 1 reports finite-domain integrals, not demonstrated infinite-space integrals. Defining mass from the core alone would be a new project convention. [Paper eq. (21), Table 1 and §4.2](https://etheric.com/Downloads/model_g.pdf#page=13).
My additional concern: integrating the G-equation gives
\[ \frac{dI_G}{dt}=\oint_{\partial V}\nabla\phi_G\cdot n\,dS-qI_G+gI_X. \]
Zero Dirichlet values do not imply zero boundary flux. Before interpreting small residual integrals physically, check domain size, boundary flux, quadrature error, and convergence time.
5. Define the benchmark measurements explicitly.
Add eqs. (20)–(21): each field uses its own first-zero radius, with geometric integration weights. [Paper §3.2](https://etheric.com/Downloads/model_g.pdf#page=13). Also specify that \(\delta(r)\nabla^2\phi\) differs from conservative diffusion \(\nabla\cdot(\delta\nabla\phi)\); they differ by \(\nabla\delta\cdot\nabla\phi\). This matters when translating diffusion into graph rules.
6. Tighten the stability conclusion — section 8.
I independently reproduced \(\max\Re\sigma=-0.195054\), \(k_*=1.9953\), and \(b_c=30.04408\). Those numbers check out. However, linear stability establishes that this operating point is below the linear instability threshold; it does not alone establish a subcritical bifurcation. That requires nonlinear branch analysis.

The reaction scheme, PDEs, steady state, reference parameters, seed, and every Table 1 entry match. Minor corrections: potential definitions are eq. (12); the secondary 2008 citation has the wrong journal metadata. [Paper](https://etheric.com/Downloads/model_g.pdf).
I have not changed the repository.
