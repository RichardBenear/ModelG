# MG5 proposal draft 1 review (MG014)

Reviewed 2026-10-31 at commit `3179342a9b559af06beabdfea98890c627a97e12`.
Proposal: spec/MG5_proposal.md, blob `9d6dd2b6800ffcbcef61923a232180c0f341c944`.
Read with InitialPlan.md, GATES.md and the MG4c closure/response. MG4c's W1/O1 and recorded numerical deviation are accepted and are not reopened here.

## Recommendation

MG5 is a useful next study: it asks whether the localized structures already found on frozen irregular graphs survive spatial turnover, and whether they become mobile. **I would not register draft 1 yet.** The main problems are the relationship between the specified event/transfer rule and its implementation cost, the changed frozen baseline, and the numerical/reproducibility rules needed for a new evolving substrate.

The PI has already selected interpolating transfer, inclusion of Part B and unweighted diffusion only. This review accepts those choices; it does not ask to reopen them or add a weighted campaign, token noise, or a larger graph panel. The appropriate scope is a fixed-panel benchmark of the specified coupled graph/field process.

## 1. Sequential transfer and the per-epoch timing model do not yet agree

The birth/death process and field-transfer rule are individually understandable. The computational plan is not yet consistent with their composition. Every birth must interpolate on the triangulation **after all preceding events**, including deaths and other births in the same epoch. The timing estimate appears to allow one global Delaunay rebuild per epoch.

At N0=18,901, the expected birth count and death count over a horizon T are each N0*nu*T, since the population mean remains N0. The expected event workload in Part A is:

| nu | Epoch length | Events per epoch | Events per 500-unit run |
|---|---:|---:|---:|
| 0.002 | 1 | 75.604 | 37,802 |
| 0.02 | 1 | 756.04 | 378,020 |
| 0.2 | 0.25 | 1,890.1 | 3,780,200 |

If reusing graphs.delaunay means a full one-second rebuild after each event, the fastest run requires roughly **1,050 hours of rebuild work alone**, rather than 1.6 hours. This is a conditional cost calculation, not a claim that efficient local updates are impossible.

If instead there is just one rebuild per epoch, explain how each birth obtains the correct containing triangle on the sequential intermediate point set. Interpolating every birth on the old epoch mesh is a different transfer algorithm. Newborn nodes can be parents of later births, or can die before the epoch ends; those histories cannot generally be discarded just because the final point set agrees.

**Before registration:** select and specify either an actual dynamic/local triangulation implementation for sequential events, or an explicitly defined batched transfer approximation. Time the implemented epoch with its real event count on an unregistered graph. Re-state the budget from that work. A solver-only timing pilot cannot establish the cost of the prescribed transfer process. Maintain stable node IDs independently of array indices when deletions reorder storage.

## 2. MG4c is not an identical frozen control after changing h and solver

The proposal changes the two graphs' h values from 0.2000738278 and 0.2000276216 to 0.2. That is a small coefficient change, but the stored states are equilibria of the old operators, not the new one. Direct static RHS evaluations give:

| Start | Max RHS at original h | Max RHS at proposed h=0.2 |
|---|---:|---:|
| D0 p0 | 6.25e-13 | 0.0490905 |
| D0 p1 | 7.96e-13 | 0.0491524 |
| D1 p0 | 8.53e-13 | 0.0184669 |
| D1 p1 | 6.82e-13 | 0.0184458 |

No integration was performed to obtain these numbers. They show that the claim of a frozen outcome 'already known exactly' is not valid for the proposed setup. They do not show that the small h change will destroy a particle.

Part A should have matched nu=0 controls using MG5's actual h, ROS2 implementation and 500-unit horizon, with the same starts and detector. Four frozen controls are a direct fix; count them explicitly (24 primary runs rather than 20) if they are added. An alternative is to retain each graph's original h throughout its lineage and explicitly establish that the new solver preserves the stored equilibrium. Do not substitute evidence from an unrelated continuum-start pilot for this baseline.

Part B's four frozen controls are correctly included, because its initial state differs from MG4's. Also define the b=28 to b=29 initialization precisely. The stored files contain potentials relative to the b=28 background. If b is physically quenched while absolute concentrations stay continuous, then

    phiY_new = phiY_old + Y0(28) - Y0(29)
             = phiY_old - 0.0642857142857143.

G0 and X0 do not change here. Reusing the old potential arrays without this shift instead resets the absolute Y background. Either can be a deliberately defined test, but they are different starts; both frozen and evolving Part B runs must use the same stated convention.

## 3. ROS2 is acceptable in principle; validate the complete algorithm

The solver switch is reasonable given restart costs. The cited ROS2 method is an established second-order L-stable method with gamma=1+1/sqrt(2): [Verwer et al.](https://ir.cwi.nl/pub/1252/) and the [maintainers' coefficients](https://kpp.readthedocs.io/en/2.6.0/tech_info/07_numerical_methods.html#ros-2).

For the usual autonomous ROS2 normalization with B=I-gamma*h*A,

    k1 = B^-1 f(y)
    k2 = B^-1 [f(y+h*k1) - 2*k1]
    y_new = y + h*(3*k1+k2)/2,

the O(h^2) contribution involving the approximate Jacobian A cancels, leaving h^2*J*f/2. Thus a correctly implemented approximate-Jacobian version can retain second-order consistency. This does not make arbitrary reuse unconditionally accurate or stable for the evolving nonlinear problem. Freeze the actual stage formulas, Jacobian refresh policy and treatment of graph changes, then test that implementation.

The frozen BDF comparison and paired epoch-halving check are sensible starting validations, with these necessary clarifications:

- Save one event stream with stable IDs and replay **exactly the same events and transfers** at both epoch resolutions. Changing resolution must not change the random process through RNG-consumption order.
- Require compatible particle counts/classes. Centre/depth differences are defined only where identities can actually be matched. Predefine how the check is evaluated if one or both resolutions dissolve or multiply; otherwise the validation is undefined for the very outcomes MG5 is designed to study. Include a field discrepancy measure on a common representation where centroid comparison is unavailable.
- On halving the fastest epoch, Delta becomes 0.125, which is **2.5 steps of dt=0.05**. Specify an aligned reduced dt or the final partial-step rule and its extra factorization. Also fix whether epoch-length changes affect the Jacobian reuse interval. The current 'one LU per epoch' cost assumes a constant stage matrix.
- State the boundary convention and t=0 centre consistently. Comparing resolutions at pre-event check times versus post-event times must not introduce a different state convention.
- Carry the already defined finite retry/stop policy into the statuses. Validation failures, unresolved balance flags, crashes and incomplete budget runs must not automatically become physical FAIL or PARTIAL outcomes.

I would accept the residual risk of a modest predeclared validation set once it exercises the actual algorithm. I am not requesting a broad solver-comparison campaign.

## 4. A1 tests rewriting plus transfer, not uniquely an averaged regular medium

Preserving a homogeneous field is a useful property of the chosen transfer rule. It does **not** imply that transfer is neutral for a particle. Repeated barycentric interpolation changes curved profiles systematically. For example, f=x^2+y^2 is 1 at all vertices of the triangle (-1,0),(1,0),(0,1), while at its centroid the true value is 1/9 and its linear interpolant is 1. This is the familiar curvature-dependent interpolation error; repeating it is not an innocuous zero-mean operation.

In the high-rate Part A runs, the expected birth and death totals each equal 100 times the initial node count. A1 could therefore reflect a combination of changing connectivity, transfer/remapping effects, and an annealed diffusivity that differs from the frozen homogenized diffusivity. The average unweighted random-graph operator is not automatically T's operator merely because one-time snapshots remain D-type.

Keep A1 as the useful operational result **'this evolving-graph and transfer model restores single-particle survival relative to its matched frozen controls.'** It is consistent with annealing, but does not identify that mechanism or show equivalence to a regular medium. A2 likewise means no restoration under this tested model/rate, not a proof that disorder cannot average out. No extra weighted copy is required merely to narrow the interpretation.

The birth/death population formula is correct for an initial iid uniform fixed-count ensemble: mean N0, variance N0*(1-exp(-2*nu*t)), and standard deviation bounded by sqrt(N0)=137.48. Conditional on the realized starting graph, the process initially retains its particular surviving points; the unconditional D-family statement is not exact resampling of each frozen realization. Its count distribution also evolves from fixed N to a Poisson limit. Describe snapshots as the same D family conditional on N rather than claiming identical full graph statistics at all times.

## 5. Store enough data for an independent review of the main statuses

Saving particle lists and a scalar minimum every 10 units, but fields only every 50 units, does not permit independent re-detection at most status-defining checks. A recompute script that rereads those same particle lists would verify bookkeeping rather than the detector's output.

At every 10-unit check save at least the smoothed X vector plus the associated node IDs/order, and enough event-stream/geometry information to reconstruct that snapshot (initial graph, full seed tuple, deterministic replay convention and snapshot hashes, or explicit positions). This does not require saving all three full fields at every check. Preserve the initial centre and float64 final state. Freeze unique stream tuples for rate, start and purpose so they do not accidentally reuse MG4 graph/placement streams.

The between-event RHS balance identity is also insufficient to audit birth/death transfer. Add event or epoch jump accounting:

    Delta I_S = +phi_S(new) on birth;
    Delta I_S = -phi_S(deleted) on death,

for potential sums at fixed parameters, with the corresponding background term when reporting absolute concentrations. Sum these jumps across each batch and verify the observed before/after change. This catches stale-parent interpolation, duplicate/deleted-node mistakes and species misalignment that a correctly coded ODE RHS cannot detect.

Include explicit INCOMPLETE / NOT SCORED handling and control-failure precedence before the PASS/PARTIAL/FAIL and A1/A2/A3 rules. Distinguish tracking-resolution failures from actual disappearance in the interpretation of unresolved runs. 'Split' and 'remote' are proximity-defined descriptions, not demonstrated mechanisms of fission or nucleation; freeze the original-identity matching rule used to assign them.

A3 can include four dissolved fast-rewriting runs, not only a mixture of different outcomes. Label it 'other outcome' or keep the full class vector visible rather than interpreting 'mixed' as heterogeneous evidence. A2 is specifically persistence of multiplication, rather than every possible failure to restore a single particle.

## Direct answers to reviewer questions 4–6

**4. Radius 1.5:** acceptable as a predeclared, descriptive movement threshold of roughly half a wavelength. Measure from each stored settled centre at t=0, not the earlier continuum-placement centre. The MG4c placement-to-pinning shifts justify a scale but are not themselves a null distribution of tracking noise on changing meshes. Interpret MOBILE as movement under this model, not proof of a response to a G gradient. That response remains MG6's question.

**5. Epochs and solver:** acceptable in principle, after resolving sequential event semantics/cost, step alignment, complete paired validation and the matched frozen baseline. The 5% turnover rule alone is not an error bound. The proposed fixed finite validation-and-stop approach is appropriate once its comparisons are defined for all physical outcomes.

**6. Four runs and 3/4 rule:** adequate for a small fixed-panel benchmark. They are not enough to claim a general survival probability of at least 75% or a universal mobility law. Shared graphs/initial states and independent rewriting streams should both be stated. Report every start's result, and treat the rate-wise single-run spread as descriptive.

The MSD fit is useful descriptively, but call it a lag-10-to-100 effective fit unless the lag dependence supports a diffusive regime. Plot/report the MSD and drift; a straight-line fit alone does not establish Brownian motion, and overlapping lag windows are not independent replicates. Mobility estimates conditioned on survival describe the surviving subset. Preserve that conditioning explicitly.

## Smaller background clarification

The Brownian-map discussion should be limited to the relevant uniform random-map ensemble and scaling limit. It does not establish that every combinatorial flip rule on an embedded starting triangulation must produce Hausdorff dimension four. The [primary scaling-limit result](https://arxiv.org/abs/1104.1606) concerns uniform planar maps with specified laws. This qualification does not prevent the proposed geometric study.

## Next step and evidence

Revise the event/transfer implementation and measured workload, establish the actual frozen baseline, and freeze the numerical and audit branches above before registration. The PI's three already accepted modelling choices can remain unchanged. I am not asking for more random realizations merely to strengthen statistical confidence.

The companion design_checks.py/design_checks.json reproduce the workload arithmetic, interpolation example, b-background conversion and static RHS mismatch. All 59 local reference inputs are hash-verified against the reviewed public commit. No MG5 trajectory or new particle outcome was generated, and no upstream file was edited. This review is saved locally under review/MG014_replication/ and has not been committed or pushed.
