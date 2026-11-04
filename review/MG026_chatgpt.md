# MG026 review of MG9 Draft 1

**Recommend Draft 2 before registration.** I support a deterministic study of the Turing normal form and front motion. The proposed methods can test local subcriticality and operational front arrest. They cannot yet identify the single-particle window with a homoclinic-snaking region or establish a physical Higgs potential. Those distinctions must be reflected in the outcomes before running the study.

Reviewed revision 99229bc0beb849f877fae17054aba50cddcb770e on 2026-11-04. My new calculations were static Jacobian/dispersion checks and run-budget arithmetic only. I did not compute L₃, evolve nonlinear modes, continue branches, or run the primary front experiment. The static script and results are in review/MG026_replication.

1. **Narrow the interpretation of P1 and P2 (§§1, 5, 11).**

   L₃ > 0 establishes local subcriticality of the selected primary Turing branch under the stated normalization. It does not establish a stable saturated pattern, a fold, a pinning interval, or snaking. The cubic equation as written has no saturating term when L₃ > 0. A stabilizing higher-order term or the full nonlinear system is needed to obtain the other state.

   Conversely, a validated L₃ < 0 would establish a locally supercritical primary branch. Coexisting finite-amplitude patterns below onset could lie on a disconnected branch or a branch that turns farther from onset. Their existence does not contradict that local result. Replace the proposed automatic contradiction/error reading with this possibility.

   Even if front arrest approximately matches the old single-dip window, report consistency with a pinning interpretation. Demonstrating homoclinic snaking requires the organization of localized branches and successive folds, which the draft does not compute. The cited Swift–Hohenberg study explicitly separates small-amplitude normal-form behavior from the large-amplitude localized states. A negative front result would likewise leave the local amplitude reduction intact. [Burke and Knobloch 2006](https://doi.org/10.1103/PhysRevE.73.056211)

   Replace “the system is driven and dissipative, so it has no energy functional” with “no global variational structure has been established for Model G.” Dissipation alone does not exclude a Lyapunov functional. A potential for a reduced real amplitude does not supply physical energy units or prove a global potential for three coupled fields. Non-variational extensions require separate analysis of local normal forms, large-amplitude branches and stability. [Burke and Dawes 2011](https://arxiv.org/abs/1108.1984)

2. **Distinguish continuous wavenumbers from the finite periodic box (§§4–5).**

   I independently obtain the following thresholds from the reaction Jacobian and the discrete Laplacian.

   | Dispersion problem | b_c | critical physical k |
   |---|---:|---:|
   | Continuum, continuous k | 30.04407748 | 2.09927549 |
   | h = 0.2, continuous k | 30.04407748 | 2.11500857 |
   | h = 0.1, continuous k | 30.04407748 | 2.10314949 |
   | h = 0.2, L = 100 discrete modes | 30.04576994 | 2π·34/100 |
   | h = 0.1, L = 100 discrete modes | 30.04757318 | 2π·33/100 |

   For these spacings the continuous-k lattice dispersion can realize the same optimum Laplacian eigenvalue as the continuum. The quoted 30.0458 benchmark belongs specifically to the finite h = 0.2, L = 100 box. It should not be used to fail a continuous-k calculation. The apparent rise on mesh refinement comes from changing the available mode nearest the optimum.

   Freeze N, h, L, mode number and the corresponding threshold for every validation domain. “An integer number of critical wavelengths” does not determine a realizable mesh at fixed h. For mesh comparisons of fronts, hold the physical pattern wavelength fixed unless wavelength dependence is explicitly part of the test.

3. **Make L₃ and its validation numerically identifiable (§5).**

   Specify the right critical eigenvector normalization, left/right pairing, Fourier convention, and whether A is a complex Fourier coefficient or a real cosine amplitude. A factor of two in amplitude changes the numerical cubic coefficient by four. Also specify original time t versus slow time T and the scaling of b − b_c. Alpha from the physical dispersion derivative must be converted consistently.

   Initialize the full critical eigenmode, with either its slaved mean and second harmonic or a defined transient-discard procedure. Measure A by projection onto the left eigenvector. A bare cosine has transient quadratic responses before the cubic reduction applies. Fit dA/dt = sigma A + L₃|A|²A using the actual residual sigma of the implemented mode. Establish that residual linear growth and solver error are small compared with the cubic signal.

   Freeze amplitudes, fitting windows, amplitude-convergence checks, solver refinement and the NOT EVALUATED rule before the sign is seen. A fitted slope at very small A can be dominated by a tiny threshold error. Keep the textbook coefficient check. In all autonomous integrations using lattice_ode_1d, explicitly set seed_scale = 0. Its default adds a localized forcing pulse, even when y0 is supplied.

4. **A time-stepping sweep does not locate a branch fold (§5).**

   Descending b in steps of 0.02 and waiting a finite time detects a protocol-dependent loss of an attracting pattern. Critical slowing near a fold, sideband instability, a wavelength change or inadequate settling can change that location. Specify stationarity and solver-success criteria, maximum settling time, how unresolved steps are classified, the initial mode, and treatment of phase slips.

   For the proposed budget, report a finite-time persistence bracket. To claim an actual fold, use steady-state continuation with a phase condition and a fold/stability check. When carrying a field to another b, preserve the intended concentrations by accounting for the changed homogeneous Y reference. Record the convention. This also makes a discontinuous field replacement distinguishable from continuation.

5. **Strengthen the front observable and classifications (§6).**

   The current categories leave gaps. A count change of one with a front displacement above one wavelength is neither pinned, invading nor retreating. Endpoints alone also miss motion followed by reversal. Internal replication can change cell count while neither front translates. Fronts may meet or fill the box before the fit window. Fix those cases in advance.

   Track both interfaces with an envelope or a phase-aware pattern-amplitude measure, unwrap their positions on the periodic domain, and retain cell count as a separate observable. Specify the pattern phase relative to the mesh and the tanh splice, and check both interfaces rather than assuming they remain mirror images. Define A*_s(b), prominence/separation rules for cell minima, and behavior when the reference periodic branch does not exist.

   Reserve an UNRESOLVED classification for slow motion, oscillations, missing fronts, inadequate settling and borderline thresholds. The proposed 500-unit window can hide a displacement approaching a wavelength. At minimum, extend apparent pinned/borderline cases under a predeclared rule, check final residuals and late-time motion, and repeat representative cases with a longer domain or different initial patch length. If no such validation is affordable, call the result finite-time arrest under the specified resolution.

   The velocity conversion dn/dt times lambda/2 is appropriate only when two fronts advance symmetrically and count changes represent addition of pattern periods. Direct interface velocities should be primary. The cubic subcritical amplitude equation alone does not supply the stable invaded state needed for the proposed front-speed prediction.

6. **Treat P2 as a comparison of measured intervals and finite-time criteria (§§2, 6).**

   The old b_s value came from a particular initial field, a centre-depth criterion and a finite t = 300 run. The old b_r came from seeded runs and the number of deep minima at t = 400. These are useful benchmarks, but their narrow bisection brackets do not establish equally precise stationary bifurcation boundaries. State that provenance and check time/detector sensitivity near the comparison edges before asserting agreement with mathematical folds.

   Retain brackets for front transitions rather than selecting an arbitrary point as an edge. Define COINCIDES only when both edge uncertainties are compatible with the tolerance in the specified manner. Add UNRESOLVED if an edge is outside the scan, a classification is ambiguous, or several disconnected arrest intervals occur. A failure to observe arrest on a coarse finite grid does not establish an empty pinning interval.

   The h = 0.1 sampling plan must be deterministic and broad enough to bracket both transitions. Define exactly which twelve points are selected. Compare each mesh with its own particle-window benchmarks. A two-spacing change of more than half is a mesh-sensitivity flag. It cannot uniquely identify grid locking, since pattern wavelength, phase, front speed resolution and ordinary truncation error also change. Targeted phase shifts and further refinement are needed for that attribution.

7. **Fix the run budget and disclosure (§9).**

   The listed 80 continuation integrations, 33 front integrations and approximately 70 seed-bisection integrations already total 183, before V2, refinement or settling runs. The cap is 120. D1 additionally asks for both mesh spacings and continuum results beyond the tasks explicitly scheduled in A2. Align the deliverables with the count and reserve refinement slots.

   The “no computation of any kind” sentence also needs qualification because §9 reports a newly timed reference run. List the timing-only work and its inputs. No scientific outcome need be hidden or relabeled.

I recommend the PI approve the following direction for Draft 2.

| Decision | Recommendation |
|---|---|
| Scope | Keep deterministic Part A. Define P1 as the local onset type and P2 as the relation between measured front arrest and the operational single-particle window. Defer a claim of snaking to branch continuation. |
| P2 tolerance | Keep 0.10 provisionally, but use explicit edge brackets, targeted refinement and UNRESOLVED. The accuracy of the observation must support the comparison. |
| Detector | Use a direct envelope/interface measurement for motion and unsmoothed minima as secondary cell counts. Freeze definitions, thresholds and phase conventions. |
| Budget | Give P1 and P2 priority, explicitly count every integration and reserve validation runs. Defer D2 if needed because its trend is already partly known. |
| Stochastic follow-up | Keep MG9b separate. A decreasing deterministic seed threshold is only a protocol-specific basin-access proxy. A stochastic action depends on the noise law and escape path as well as deterministic dynamics. |

This stage can make the earlier analogy much more precise by determining how b changes a critical mode and whether fronts arrest near the known particle window. It should allow either local bifurcation sign and either front outcome to stand as a valid finding.
