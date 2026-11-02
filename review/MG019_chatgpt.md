MG019 — MG6 proposal Draft 1 review

Reviewed [Draft 1](https://github.com/RichardBenear/ModelG/blob/53a2adc89e330d2ca2686af9bde1665b0f57003e/spec/MG6_proposal.md) at commit `53a2adc89e330d2ca2686af9bde1665b0f57003e`, with the design calculation, roadmap, glossary and MG5 findings.

**I support the investigation, but would revise this draft before registration.** Its strongest feature is the explicit acknowledgement that suppressing births in a G-well largely builds in a density deficit. The interesting questions are whether the particle survives the feedback, how large the response is, and whether the controls and quantitative predictions distinguish the intended mechanism.

**Changes needed before registration**

1. **Specify a valid coupled birth–death process.** “Deaths are the MG5 stream” is insufficient once births are accepted or rejected. An MG5 death can name a node that MG6 never created; F1 can also create nodes absent from the original MG5 process.

   A workable construction is a dominating candidate-birth process, with an independent exponential lifetime of rate ν assigned to every initial node and every candidate newborn. Accept births using the boundary field; ignore deaths of rejected candidates. Equivalently, simulate an upper process containing rejected nodes as ghosts, while distinguishing its population from the actual population. Define persistent IDs, ephemeral nodes, candidate marks and accepted-set records precisely.

   This can provide common random numbers across arms, but hash-identical reproduction of the old MG5 generator is a separate compatibility requirement. The current text does not establish it. Preserve MG5's tick-boundary transfer convention explicitly.

2. **Separate a node-density proxy from curvature, and define the integration measure.** Fixed h in the graph operator does not by itself prove that node counts inside an embedding-coordinate disc measure intrinsic curvature. The disc radius uses stored coordinates, whereas the proposed physical area uses node count. Those are different constructions.

   I recommend calling the primary quantity an **embedded node-density deficit**, adopted as a volume proxy. A curvature interpretation needs an additional metric/volume argument or graph-distance measurement. This does not require adding that investigation to the first campaign.

   Also distinguish the birth-rate reference density, ρ₀=N₀/A, from the snapshot mean, ρ̄(t)=N(t)/A. Record both:
   
   D₀(R)=N(<R)−ρ₀πR², and D_c(R)=N(<R)−ρ̄(t)πR².

   They answer different questions: absolute change versus local redistribution relative to the current background. F1 generally changes the total expected population. In my static starting-state calculation, subtracting the predicted global mean changes the κ=4, R=8 excess from roughly 29–35 nodes to 14–17. This is material for the outer profile.

3. **Correct the coupling arithmetic and the prediction.** F2's statement that clipping is inactive at the proposed strengths is false. At κ=4, the minimum of 1+κφ_G is approximately −0.105, −0.107, −0.053 and −0.020 in the four starting states. Respectively 13, 14, 6 and 1 nodes already require clipping. The exact zero-net-change argument therefore does not apply to the clipped F2 panel.

   F1 at κ=2 and 4 is also not a linear-response test. Using the registered particle detector to centre the starting states, the predicted core-deficit ratio is about **1.70–1.73**, rather than 2. That is consistent with the exponential response; it is not an error in F1.

   For later coupled states, do not approximate a uniform-coordinate area integral by assigning A/N to every node: the nodes are deliberately sampled nonuniformly. Use fixed-coordinate quadrature with barycentric field evaluation, or justified periodic cell-area weights. Otherwise the density variation biases the integral itself. The stationary graph identity Σφ_G=0 concerns the node measure; it does not automatically establish a zero coordinate-area integral on a nonuniform mesh.

   Finally, ρ₀f is an equilibrium target for a prescribed, stationary field. For a prescribed time-dependent field, the mean density obeys ∂ₜn=ν[ρ₀f−n], so it depends on field history. Treat the instantaneous target as a quasistatic approximation, not an exact prediction for a moving, coupled particle. Agreement with a prediction based on the particle's own reshaped field would not rule out feedback.

4. **Replace the current statistical decision rule.** Four run averages do not give a calibrated “3σ” test merely by dividing their mean by an estimated standard error. Even under an ideal independent Gaussian four-run model, a statistic of 3 has a two-sided Student-t p-value of about **0.058**, not 0.0027. The approximately ten independent snapshots per run are also a power assumption, not something already established by MG5.

   Average within each run, then analyse the four run-level results. If arms share candidate streams, use paired treatment-minus-control differences and their covariance rather than an independent-arm combined SE. Report the individual contrasts and explicitly conditional uncertainty intervals. If formal significance is required, register its statistical assumptions and calibration.

   A nonzero decoupled density profile is not automatically an INVALID-NULL: the detected centre depends on fields that themselves respond to the graph, so particle-centred sampling need not be spatially unbiased. The treatment-minus-null contrast is the relevant test. Reserve invalidity for specified failures of the control construction or implementation.

   Failure to establish a deficit should be **INCONCLUSIVE / NOT DEMONSTRATED**, unless a separate equivalence criterion establishes that any deficit is small. It does not establish NO DEFICIT.

5. **Fit falloff to a defined local quantity.** ΔN(R) is cumulative. A localized density perturbation can produce a plateau, a zero crossing or cancellation in this quantity, rather than the same decay as φ_G. Taking its absolute value can hide the compensating halo. In two dimensions, dΔN/dR=2πR·δn(R); a power exponent fitted to the cumulative count is not interchangeable with a density exponent.

   Keep the cumulative radii, but measure annular density and radial φ_G as well. Prespecify the fit range, treatment of signs/zero crossings and minimum signal requirement. If the outer signal is unresolved, report an unresolved exponent rather than selecting a fit. Neither an exponential nor a power law should be a foregone conclusion.

6. **Complete the numerical and audit rules.** MG5's full-event replay is no longer automatic: changing dt or Δ can change acceptance decisions through the field. Specify what is held fixed in each refinement check. A fixed accepted-edit history can isolate ODE time-step error; a fully coupled refinement checks the complete algorithm and may need a different comparison criterion.

   Before production, define quantitative Part 0 tolerances and stop/fallback rules, the F1 bound and breach policy, and the numerical/incomplete status precedence. Multiplication or dissolution under coupling is a physical feedback result, not “no deficit”; do not silently discard those runs or average only survivors.

   Commit enough new measurement data to recompute the claims: centres/classes, annular counts and areas, both density baselines, G-deficit integrals, quadrature definitions and predictions, and accepted-candidate identifiers or bitmaps. Candidate hashes alone cannot reconstruct field-dependent acceptances. Full fields can remain local under a checksum manifest.

**My recommendations for the Part 9 decisions**

| Decision | Recommendation | Why |
|---|---|---|
| 1. Coupling function | **F1**, with a registered dominating bound and explicit bound-failure policy. | It is smooth, positive and retains both the well and halo. F2 clips at κ=4; F3 deliberately removes the halo and changes the question to a one-sided coupling. |
| 2. κ panel | **Keep {2,4}: κ=2 as the primary response arm, κ=4 as a stronger feedback/stress arm.** | The expected core birth factors are about 0.58–0.60 and 0.33–0.36 respectively. Both are informative, but neither the panel nor its ratio should be advertised as a small-κ linearity test. |
| 3. Randomized null | **Fixed spatial shift, matched at both κ values; defer shuffled.** | A shift retains spatial smoothness while breaking alignment with the particle. Shuffling changes spatial correlations and interpolation statistics as well, making it a different null. Specify the periodic offset distribution, separation and streams in advance. |
| 4. Integrated deficit / mass and b-scan | **Define an operational well integral M₋=∫max(−φ_G,0)dA; also report the positive halo and signed integral. Defer the b-scan to a separately registered scaling study.** | This measures a nonzero well without hiding compensation. For the proposed embedded-density experiment, use coordinate-area weights consistently. Call it a G-well integral, not established gravitational mass. Varying κ tests coupling strength; it cannot establish scaling with particle mass. Deferring the b-scan must be recorded as deferring that roadmap objective. |
| 5. Decoupled arm | **Rerun under the MG6 pipeline, with declared MG6 streams and matched starts.** | This tests the complete new code path and allows deliberate pairing with the new arms. Reusing MG5 remains useful as an external compatibility check. Repeating the same seed does not create independent evidence. |
| 6. Rates | **ν=0.2 only in the first campaign.** | It has MG5 survival support at the highest tested turnover, a five-unit density relaxation time, and an established cost baseline. Adding ν=0.02 introduces slower equilibration and another interpretation before the primary coupling is validated. |
| 7. Thresholds, radii and time window | **Keep R=1.5 primary, retain the other radii descriptively, and provisionally keep t_end=200/window 100–200. Replace the uncalibrated 3σ/NO DEFICIT rule as above.** | The core has the best predicted signal. The time window is a reasonable starting design at ν=0.2, conditional on Part 0 checking transient effects and temporal dependence. More snapshots are not more independent runs. |
| 8. Budget | **20 main runs, cap 24, provisional 40-hour launch cutoff, restated after full-pipeline Part 0 timing.** | Coupled κ=2/4, shifted κ=2/4 and decoupled give five arms × four starts. The MG5-rate estimate is about 30 ledger hours, leaving four run slots as headroom. A launch cutoff is not a hard spending ceiling; include event generation/replay and any interrupted work in the accounting. |

This gives a focused first experiment: **does the chosen G-dependent birth rule produce the predicted density response while preserving or changing the particle?** It is a worthwhile question. It should not yet be presented as a completed demonstration of curvature, a gravitational falloff law, or particle-mass scaling.

I performed static checks only, using the four stored MG4c states; no MG6 dynamics or event streams were generated. The checks and pinned inputs are documented in `review/MG019_replication/`. No upstream proposal or gate was changed.
