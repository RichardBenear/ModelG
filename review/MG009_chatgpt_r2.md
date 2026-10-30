# MG009 r2 — ChatGPT review of the MG4 proposal, draft 3 (verbatim, relayed by the PI 2026-10-30)

Draft 3 resolves my main objections. I support its scientific design and benchmark approach, with a few definitions to finish before registration.
The answers to its two open questions are yes:

* C1’s panel counts are reasonable: ≥28/32 overall with ≥6/8 on every graph is a clear benchmark, without an unsupported population-confidence claim. C3(b)’s 10% PASS / >20% FAIL bounds are defensible predefined tolerances.
* The band-averaged spectral count is sound for its stated purpose. Save all 21 ratios so cancellations across the band remain visible. It measures mode counts; C3(b) separately tests operator action.

Four small clarifications remain:

1. Define C2’s “median.” There is no uniquely ordered “median run” across several metrics. Specify componentwise medians of per-run errors, require both depth measures to pass, and define missing wavelength as a failed wavelength comparison.
2. Fix the eigenvalue sign convention. C3(a) should count nonnegative eigenvalues of
\[ -L_G=(D_{\rm deg}-A)/h^2. \]
The dynamics and C3(b) use the negative-semidefinite generator \(L_G\).
3. Make the balance tolerance safe near zero. A purely relative residual is undefined for the exact homogeneous state and unstable near stationarity. Use an explicit absolute-plus-relative tolerance. Also, a pointwise RHS below \(10^{-8}\) does not automatically guarantee the global zero-sum target; retain the summed derivative terms when interpreting failures.
4. State gate precedence. Specify what happens when one criterion FAILs while another is INCONCLUSIVE, or a positive control fails. Resolve that ordering before results exist.

Two interpretation limits should remain explicit: C3(b) tests coarse-grained action on the chosen sampled profiles, rather than universal operator convergence; and finite-horizon survival alone does not establish stationarity.
The full response is in [Draft 3 review](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG009_replication/draft3_review.md).
These are specification edits, not a request for another statistical study or additional graph placements. Once incorporated, I’m comfortable with registration and the planned Part 0 controls.
