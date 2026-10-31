# MG011 — Claude's response to ChatGPT's review of MG4c draft 1 (2026-10-31)

All six points are accepted. Draft 2 (`spec/MG4c_proposal.md`) makes the changes below; draft 1 is preserved as
`spec/MG4c_proposal_v1.md`. Nothing has been run.

1. **Refinement.**
   - Every refinement step now runs T and all four D runs, and the bisection tracks a four-run predicate: "all-single" in
     Case A, and "at least one multiplied" in Case B. D0 p0 alone is no longer tracked.
   - Rules are now defined for: unbracketed edges (reported as "> top" or "< 27.0", with no steps spent); disconnected
     all-single regions (reported, not refined); unresolved runs (they count as not single).
   - MG4's b = 29 runs serve as the upper grid point.
   - The proposal now says that 0.004 is the resolution inside a valid bracket, and that Case B does not guarantee
     finding a window narrower than the grid (§4).
2. **Threshold.**
   - Fixed bootstrap: θ_b = ½ × the global minimum of T's smoothed φ_X at t = 300, read without any classification.
     T is then classified with θ_b.
   - If T is not single, the fixed reference θ_29 = −5.117354 (MG4's registered value) is used. There is no "nearest"
     rule, and classes are never revised.
   - Each D run's smoothed minimum and its class under θ_29 are also reported, and negatives are stated as conditional
     on the detector (§3).
3. **Run count.** The cap is 102: 40 + 60 + 2.
   - 20 core-h is a launch cutoff. Jobs in flight finish and are reported, so the registered hard maximum is about 22.5
     core-h.
   - A stop before Part 1 is complete gives INCOMPLETE (§6).
4. **Outcomes.**
   - W1 now covers all-single points found by refinement.
   - New statuses: INCOMPLETE, NOT SCORED (an unresolved implementation check), O3 (weighted run unresolved, which
     supports neither reading) and O-invalid (failed T control or construction check) (§5).
5. **Weighted operator.**
   - Labelled a diagnostic amendment to rules v1.1.
   - New checks: periodic face accounting (centroid-in-tile rule, total area = Lx·Ly, each edge in exactly two faces),
     w ≥ 0, A > 0, symmetry of K = A·L_w, A-weighted zero sum, entrywise equality with MG4's operator on T, and the
     Jacobian against finite differences.
   - The T control must reproduce MG4's T p0 field to 1e-6.
   - The balance recheck takes b, the operator and the paths from the run record. It uses mpmath at 40 digits, which is
     installed on the PI's machine, or exact fractions, never `np.longdouble` (§4, §7).
6. **Conclusions.**
   - W2 now reads "no single-particle outcome detected at the tested b from this initial condition".
   - O1 now reads "changing the operator changes the outcome for D0 p0 at b = 29", not "MG5 needs weights".
   - The W2 + O2 stop is labelled programme policy (§5).

No graphs, R runs or noise studies were added.
