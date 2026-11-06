# MG041 — Proposed correction to availability qualification

2026-11-06. **ADOPTED by PI, 2026-11-06: “yes, go ahead”. No scientific production dispatched at adoption.**

The PI authorized proceeding without external review: “No review needed. Go ahead.”
That decision clears external review for this first batch; it does not make a failed
mandatory qualification pass. MG040 section 4 requires stopping on that failure.

## Defect and evidence

MG040 section 3 check 4 requires every deterministic A availability path to match
exp(-t/(C*T_ref)) to absolute 1e-8. At P17, q=p, availability obeys
r'=-r/(C*T_ref), independently of all active fields. The specified Euler method
instead gives r_n=(1-dt/(C*T_ref))^n. Its truncation error is nonzero even in
exact arithmetic. The qualification incorrectly conflates discretization accuracy
with implementing the prescribed discrete update.

The literal Euler recurrence was evaluated at all integer times 0..300 for every
capacity and all three prescribed step sizes. At C=1, maximum errors against the
continuous exponential are 4.0875530671e-7, 2.0437757142e-7 and 1.0218874230e-7
for E0/E1/E2. E0 at C=10 also fails (1.0053795352e-8). Four of 21 paths fail.
The largest error versus the correct discrete formula across all paths is below
5e-11, including accumulation of floating-point rounding at the largest capacities.
Decimal-arithmetic recomputation independently confirms every pass/fail label.

Raw evidence: results/MG10/qualification/availability_checkpoints.npz,
availability_qualification.json, availability_recomputed.json.
The run plan itself is unchanged and its hash is retained in the qualification record.

## Recommended replacement for check 4

Separate the numerical methods' exact cases:

1. FINITE Euler E0/E1/E2: availability at each integer checkpoint must agree with
   (1-dt/(C*T_ref))^(t/dt), absolute error <=1e-8. Compare to a stable log1p or
   higher-precision oracle, never subtract a tiny flux from enormous A stocks.
2. FINITE ODE reference O: availability must agree with exp(-t/(C*T_ref)),
   absolute error <=1e-8 at every saved checkpoint.
3. CHEMOSTAT: availability must remain one. The finite exponential formula does
   not apply to an externally replenished bath.
4. Preserve and report the Euler-versus-continuous errors. Keep the separately
   proposed E2/O active-field numerical agreement threshold of 0.001 and recovery
   margin of 0.01 unchanged. A correct discrete update alone does not establish
   an adequate approximation to the continuous model.

This alters a mandatory qualification definition, so it is presented explicitly
for PI adoption before registration or further scientific work. It does not change
the physical model, capacities, solvers, horizon, production matrix or budget.
No particle outcome has been generated or inspected. No recovery verdict exists.

## Next steps after adoption

Implement the full deterministic model and remaining qualification suite; measure
the two timing cases; verify durable retention; publish the cost projection;
register the corrected design prospectively; dispatch only if all checks pass and
the fixed matrix fits the approved ceiling with its reserve. No additional external
review is required by the PI's current decision.

Adoption record: Richard explicitly approved adopting this correction and continuing implementation, qualification, benchmarking and the authorized batch. The earlier MG041 failed qualification remains preserved under its original criterion; it is not rescored. New checks use the adopted criterion.
