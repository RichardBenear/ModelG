# MG029 round 2 — MG9c Draft 2 review

**Recommendation: Draft 2 is ready for PI registration. I found no remaining substantive design blocker.** The outstanding limitations are now appropriately stated as limits of the experiment, rather than presented as conclusions it could establish.

Reviewed `spec/MG9c_proposal.md` and `review/MG029_response.md` at revision `4ae861362fcc20f647ce87b7ba49e951d18d4d5a`, 2026-11-04, against the MG029 review and independent checks. This is a design review. No new Model G integration was performed, and the proposed MG9c implementation and tests are not yet present at this revision.

## Disposition of the previous findings

- **Prior evidence:** the cell-count correction is accurate, including the change to 19 at b = 29.35 from t = 45. The draft no longer assumes that all observed motion was edge-cell relaxation. The independent fits from the first review remain applicable; no refit is needed.
- **Q2:** interval arithmetic correctly distinguishes proximity from physical equality. An edge slightly above the replication benchmark is no longer excluded by assertion. The local ordering assumption is explicit, and Q2 requires a confirmed SETTLED-ROBUST result at 29.35.
- **Classification:** opposite-sign motion is covered; terminal events start at t = 0; structural checks retain their separate window. The main/domain/mesh extensions and fixed perturbation horizon are explicit. Inconclusive preparations and robustness trials have defined statuses. Replay mismatch blocks confirmation, and 29.22 is correctly exempt from a nonexistent baseline comparison.
- **Interpretation and controls:** BOUNDED is an operational finite-time outcome at three specified b values. The domain control holds patch width fixed. The two amplitude perturbations are explicitly limited tests, and the stronger stability, depinning and snaking claims have been removed. The wider patch remains descriptive.
- **Auditability:** the added fields at scoring boundaries and final horizons, together with tracks, hashes and completion records, address the principal audit limitation of MG028.

The point-status ordering in section 5 is reasonable: one resolved FAIL suffices for NOT-ROBUST even if another trial is unresolved. The unresolved trial must still be reported individually; it is not evidence of a physical failure. If there is no resolved failure but a required check is inconclusive, ROBUSTNESS-UNRESOLVED correctly prevents BOUNDED.

## One nonblocking wording correction

Section 10 describes MIXED as a settling region ending inside 29.25–29.35 or a robustness failure. The section 7 catch-all also covers other resolved outcomes, including RECEDES or RETREATS, and does not by itself establish a monotone transition within that interval.

Suggested replacement: **“MIXED: the three points do not jointly meet BOUNDED or SLOW INVASION, with no unresolved condition taking precedence. Report each dynamical outcome and robustness status; do not infer a unique transition from the aggregate label alone.”**

This is an interpretation correction, not a request for a new experiment or a scoring change.

## Recommendation to the PI

Approve the proposed 12000 horizon and 30000 extension, the settling tolerances, the 0.01 operational proximity tolerance, and the descriptive 29.22 run. The design now tests the intended question at modest cost while preserving the distinction between finite-time settling and asymptotic stability.

Registration can proceed with the small MIXED wording correction above; another design-review round is not needed for that edit alone. Before production, complete the promised synthetic tests and freeze the pre-run implementation specification. In particular, preserve the original MG9 replay procedure and document any solver segmentation used for reproducing [0,1800]; replay should be validated rather than weakening its tolerance after a mismatch.

No new replication script accompanies this round: the original `review/MG029_replication/` contains the relevant independent fits and counterexamples, and this review checks how Draft 2 resolves them. Acceptance of the design is not verification of code that has yet to be written.
