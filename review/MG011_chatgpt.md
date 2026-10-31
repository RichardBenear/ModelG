# MG011 — ChatGPT review of the MG4c proposal, draft 1 (verbatim, as pasted by the PI in the project thread, 2026-10-31 16:07 UTC)

**I support MG4c as the next investigation, but draft 1 needs revision before registration.** The problems are in the decision rules and interpretation; they do not require a larger experiment.

I reviewed [draft 1](https://github.com/RichardBenear/ModelG/blob/7d2b174800adb5f12809a3e1c5867ba35709b4f8/spec/MG4c_proposal.md) and checked several decision examples independently.

1. **The refinement does not locate the claimed panel boundaries.**  
   Brackets are selected using “all four D runs are single,” but bisection follows D0 p0 alone. D0 could remain single at both bracket endpoints while another placement changes class. Either bisect the four-run predicate or explicitly report boundaries for D0 p0 only. Also define handling for unbracketed edges, disconnected single-particle regions and unresolved runs. The stated 0.004 resolution applies to a valid bracket; it does not guarantee discovery of a narrow window.

2. **The threshold rule is circular.**  
   T must be classified as single before its depth supplies the threshold, but classification itself requires a threshold. Specify a fixed bootstrap rule. “Nearest valid threshold” also needs a fixed reference set and tie rule, so later refinement cannot silently change earlier classifications. Below T’s window, a borrowed threshold could miss a shallower D particle; negative findings must remain conditional on the detector.

3. **The maximum run count is 102, not 90.**  
   Section 3 requires T at every newly tested b, including refinement:

   | Part | Maximum runs |
   |---|---:|
   | Coarse scan | 40 |
   | Refinement: 12 × (T + four D) | 60 |
   | Weighted diagnostic | 2 |
   | **Total** | **102** |

   Any additional window-discovery steps must also fit the registered allowance. Clarify whether 20 core-hours is a launch cutoff or a hard total including running jobs.

4. **The outcome categories have a gap.**  
   If no coarse point is single, but refinement finds four single outcomes at one b, none of W1–W3 applies as written. W1 should include eligible refined points. Incomplete studies, failed implementation checks and failed weighted controls also need explicit statuses. A weighted *unresolved* result should not support “neither works.”

5. **The weighted operator is mathematically sound, but its implementation checks need strengthening.**  
   Its triangular normalization is correct, and area-weighted conservation is appropriate for this lumped finite-element operator. This agrees with the [standard operator/mass-matrix construction](https://libigl.github.io/tutorial/#mass-matrix). However, Part 3 changes the diffusion model relative to rules v1.1; it should be labelled a diagnostic amendment. Add periodic-face accounting, total-area and weighted-symmetry checks.

   The existing balance rechecker also hardcodes MG4’s paths, **b=29 and the unweighted operator**. It must check the actual MG4c run and operator. On the PI’s Windows installation, `longdouble` supplies no extra precision; use a genuine higher-precision fallback.

6. **Narrow the conclusions.**  
   W2 means “no single-particle outcome detected at the tested parameters from this initial condition,” not “no window exists.” W2+O1 would show that changing the operator changes the outcome for D0 p0 at b=29—not that every MG5 mapping requires weights. A stop/rethink decision is reasonable programme policy, but these limited tests cannot establish a general obstruction.

**I would support registration once these points are fixed, without adding graphs, R runs or noise studies.**

The [full review and reproducible decision checks](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG011_replication/review.md) are saved locally under `review/MG011_replication/`; they have not been committed.
