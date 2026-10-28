# ChatGPT review of Amendment 3 and the glossary — verbatim
Received 2026-10-28 via PI. Reviewed @ 8b2f742.

**I support the revised Ω scan and the exploratory spectrum study. The confirmation change needs a clearer operational definition and control checks before approval.** The glossary is useful, but a few entries currently state interpretations more strongly than the evidence supports.

For [Amendment 3](https://github.com/RichardBenear/ModelG/blob/main/spec/phase2_amendment_3.md):

1. **Correct the large-Ω statement.**

   The physics note says particles at large Ω do not survive long without decaying. That contradicts the high-Ω survival diagnostics and the independently verified deterministic fixed point.

   The supported expectation is: **increasing Ω suppresses intrinsic fluctuations, making spontaneous nucleation rarer and noise-induced loss less likely for a stable seeded soliton.** Finite runs cannot establish that spontaneous nucleation never occurs.

2. **V7 establishes detector crowding, not yet stochastic Turing patterns or a nucleation window.**

   “Consistent with stochastic Turing patterning” is appropriate. “Nucleation is everywhere” and “isolated rare nucleation exists only around Ω=100–300” are premature: V7 confirmed **zero events under the registered rule**, and Ω=100’s persistent detections remain mostly unresolved.

   The cited Brusselator paper establishes the phenomenon in its model, providing a useful precedent rather than verification of Model G’s mechanism. [Biancalani, Fanelli & Di Patti](https://arxiv.org/abs/0910.4984).

   Keep E1. Until its results are available, describe these as **crowded, transitional, and quiet detector regimes over the observed exposure**.

3. **The replacement confirmation rule measures finite-time persistence, not proven basin membership.**

   Removing a whole-domain stationarity requirement is reasonable when several structures interact slowly. However, a deep minimum present at eleven times could still be a long transient, and the minimum need not belong to the same structure each time.

   Before freezing the rule, check it against:

   - a known stationary single soliton;
   - a known slowly evolving multi-soliton state;
   - subthreshold perturbations that eventually decay.

   Call the outcome **“deterministically persistent particle-depth structure over continuation times 200–300.”** Treat basin membership as its interpretation, with the finite-horizon limitation stated.

   Also distinguish **confirmation of the whole field** from confirmation of an individual detection. V7’s JSON contains repeated identical continuation outcomes for multiple detections from the same snapshot. Those are not independent confirmations.

4. **The new screen is a defensible sampling rule, but needs precise application.**

   More than two detections per snapshot is a reasonable **crowding exclusion** for an isolated-event study. It does not itself diagnose noise-sustained patterning.

   Apply the screen to every added Stage-B candidate, including inserted midpoints, and prohibit excluded cells from re-entering the rate fit through the neighbour-selection rule. Complete the required clip pilots for newly introduced Ω values.

   Stopping descent is an acceptable predeclared budget choice; it does not prove that every lower Ω is more patterned. If fewer than three usable rate cells remain, INCONCLUSIVE is appropriate.

For the [glossary](https://github.com/RichardBenear/ModelG/blob/main/spec/GLOSSARY_AND_ASSUMPTIONS.md), I recommend three corrections:

- **Separate the limits.** At fixed h and τ, Ω→∞ gives the Euler lattice map. Recovering the continuous-time lattice ODE additionally requires τ→0; recovering the continuum PDE requires spatial refinement. The stationary fixed-point equivalence remains valid under the stated no-clipping conditions.
- **Qualify the noise formula.** √n fluctuations are a Poisson benchmark, not a general property of this interacting system. The homogeneous X LNA variance is substantially larger than its Poisson value. The weak-noise \(1/\sqrt{\Omega}\) scaling is the relevant general statement here.
- **Keep hypotheses labelled consistently.** Change “stochastic Turing pattern—seen in V7” to “candidate interpretation of V7, pending E1.” Keep the physical Ω estimate explicitly conditional and its density citation unverified.

One carry-over inconsistency also needs removal: protocol v1.3.2’s V2 section states the amended bound-based rule, then repeats the old “FAIL if any point estimate lies outside its margin” sentence.

**Recommendation:** approve the scan refinement and E1; finalize the confirmation wording, controls, and factual corrections before using Amendment 3 to classify gate events. These are targeted changes, not a reason to reopen the accepted integrator or statistical framework.