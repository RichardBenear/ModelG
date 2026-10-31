# MG012 — ChatGPT review of the MG4c result (verbatim, pasted by the PI 2026-10-31 20:21 UTC)

Full review and replication files: `review/MG012_replication/` (committed separately by the PI).

---

**I independently reproduce W1 and the weighted single-particle outcome reported as O1. The scientific conclusions are supported; no new simulations are needed.**

Across all **102 saved runs and 1,122 observations**:

- All four D starts are single at **13 tested b values**, spanning 27.265625–28.7265625.
- Both boundary brackets reproduce exactly:

  | Boundary | Bracket |
  |---|---|
  | Lower | [27.26171875, 27.265625] |
  | Upper | [28.7265625, 28.73046875] |

- Every new D classification agrees under the adaptive and fixed MG4 thresholds.
- Weighted D0 p0 remains single at b=29.
- All independently recomputed balances pass. The ledger confirms **11.327 core-hours**.

**One numerical-tolerance discrepancy should be documented before formal O1 sign-off.** The proposal states entrywise agreement within \(10^{-12}\), while the code uses a scaled tolerance. The measured T-operator difference is \(3.50\times10^{-12}\): it passes the implemented check but misses a literal absolute reading. The final T fields agree within \(4.55\times10^{-13}\), strongly supporting numerical equivalence. An explicit numerical resolution is appropriate; rerunning the experiment is unnecessary.

Keep the interpretation limited:

- **63 of 72 unweighted singles are stationary** under the registered cutoff.
- The weighted D single is **not yet stationary**.
- Both boundary-limiting single runs remain unsettled, especially the lower one. These are finite-time detector boundaries, not precisely established bifurcation points.
