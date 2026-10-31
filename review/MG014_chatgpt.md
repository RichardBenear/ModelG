# ChatGPT review of MG5 draft 1 (MG014)

*Saved verbatim from the PI's message in the project thread, 2026-10-31 21:28 UTC. The reviewer's full review and
reproducible checks are in `review/MG014_replication/review.md` on the PI's machine, to be committed by the PI.*

---

MG5 is a worthwhile next step, but I would revise Draft 1 before registration. The PI’s choices—interpolating transfer, Part B, and unweighted diffusion—can remain unchanged.
The main issues are:

1. The event rule and cost estimate do not yet match.
At the fastest rate, one Part A run contains about 3.78 million birth/death events, roughly 1,890 per epoch. Births must interpolate on the triangulation after preceding events. The budget appears to allow one global rebuild per epoch.
If every event requires the quoted one-second rebuild, that is about 1,050 hours per run. Efficient local updates could avoid this, but their implementation and timing must be established. Interpolating all births on the old epoch mesh would instead define a different transfer algorithm.
2. MG4c is not an identical frozen control after changing h.
I evaluated the stored initial states without integrating them. Their maximum derivatives change from roughly \(10^{-12}\) under their original operators to 0.018–0.049 with MG5’s proposed \(h=0.2\).
Include matched frozen controls using MG5’s actual spacing, solver and horizon, or retain the original spacings and validate preservation of those equilibria. The existing MG4c result alone cannot establish this baseline.
Also define Part B’s parameter change: preserving absolute concentrations when raising b requires shifting the stored Y potential by −0.0642857. Reusing the potential arrays unchanged gives a different initial condition.
3. ROS2 is acceptable in principle; the complete validation needs clarification.
The proposed method is an established second-order, L-stable scheme. [ROS2 documentation](https://kpp.readthedocs.io/en/2.6.0/tech_info/07_numerical_methods.html#ros-2)
The paired epoch check must replay exactly the same events and define comparisons when particles dissolve or multiply. There is also a concrete alignment issue: halving the fastest epoch gives Δ=0.125, which is 2.5 steps of dt=0.05. Specify the partial-step or aligned-step policy and its factorization cost.
4. A1 would demonstrate recovery under rewriting plus transfer—not uniquely annealing.
Barycentric interpolation preserves constant fields but systematically changes curved profiles. Repeated remapping can therefore affect particle survival independently of changing connectivity. Keep the operational A1 result, but describe “averaged regular medium” as a possible explanation rather than an established mechanism. No weighted campaign is needed merely to narrow that claim.
5. The saved data must support independent classification.
Particle lists every 10 units, with fields only every 50, leave most decisive checks unauditable. Save at least the smoothed X vector at every check, with reconstructible geometry and node ordering. Add birth/death jump accounting; between-event balance identities cannot detect incorrect field transfer. Define numerical-failure and incomplete-run statuses before applying the scientific gate outcomes.

For your explicit reviewer questions:

* The 1.5 mobility radius is reasonable as a declared descriptive threshold measured from the settled starting centre.
* Four runs per rate are adequate for this fixed-panel benchmark, without claiming a general 75% survival probability.
* The solver and epoch approach are reasonable conditionally, once the implementation, controls and validation above are fixed.

The [full review and reproducible checks](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG014_replication/review.md) are saved locally. No MG5 trajectories were run, and nothing has been committed.
