# MG030 centre-diagnostic review — Claude's response

*Claude, 2026-11-05. Review: `review/MG030_centre_diagnostic_review.md` (revision 2417833).*

| # | Item | Disposition |
|---|---|---|
| 1 | One matched subset for the mean and the spread | **Accepted; my error.** `src/diag_mg9c_centre.py` now stores the matched identities (`matched_min_ref_x`, `matched_min_x`). All late statistics (mean, n, sd, range, maximum deviation) come from the interior minima matched to t = 1800 that are present at both times. `all_min_x` keeps every detected minimum. |
| 2 | sd is not a bound | **Accepted.** The README now gives the sd (about 1e-4) and the maximum deviation (2.1–3.5e-4), and drops "each … to within 1.5e-4". |
| 3 | Scope of the synthetic precision | **Accepted.** The README states that it applies to the synthetic profiles only, and cites the 1.5e-4 Fourier-location sensitivity. |
| — | Rebased plot | **Done.** `centre.png` has a second row rebased at t = 15,000. |
| — | "Predominantly translation"; unexplained-fraction definition | **Adopted** in the README. |

The conclusion is unchanged, and no MG9c classification changes.
