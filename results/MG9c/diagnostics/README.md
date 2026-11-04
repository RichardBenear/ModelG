# MG9c — exploratory centre diagnostic (post-result; no integration; no classification changes)

Code: `src/diag_mg9c_centre.py`. Request: `review/MG030_amplitude_center_prompt.md` (committed by the PI).
Output: `centre.json` (all rows and validation) and `centre.png` (centre changes against t). It uses the committed
snapshots only, at **t = 1800, 3600, 6000, 12,000, 15,000 and 30,000**. These sparse times are a real limit: nothing
between them is seen. Conventions are in the code docstring. The periodic cut at x = ±50 lies in the background, where
the E² weight there is below 2e-8 of the total, so no unwrapping is needed.

## Synthetic validation (before reading the runs)

| Case (true shift) | E²-weighted centre | φ_X²-weighted centre | Interior minima | Profile alignment |
|---|---|---|---|---|
| Translation +0.0083 | +0.008300 | +0.008300 | +0.0083004 | +0.0083000 |
| Translation −0.0083 | −0.008300 | −0.008300 | −0.0083004 | −0.0083000 |
| Translation +0.05 | +0.050000 | +0.050000 | +0.0500005 | +0.0500000 |
| Symmetric amplitude relaxation (0) | 0 | 0 | 0 | 4e-10 |
| Asymmetric amplitude change, peaks fixed (0) | **+0.0358** | **+0.0345** | 3e-5 | 2e-5 |

- Translation is recovered to better than 1e-6 by every measure, about 1e-4 of the 0.008 late drift.
- The last row shows that a weighted centroid can move through amplitude redistribution alone, here by 4× the late
  drift, while the minima and the profile alignment stay put. Centroids alone cannot tell the two apart; minima and
  alignment can.

## Late interval [15,000, 30,000]

| b | Edge midpoint | E²-weighted | φ_X²-weighted | Interior minima, mean (spread of individual shifts) | Profile alignment (unexplained part of the change) | Width change |
|---|---|---|---|---|---|---|
| 29.25 | −0.0081 | −0.0085 | −0.0081 | −0.0084 (−0.0085 to −0.0081, sd 1.2e-4) | −0.0084 (1.5%) | +0.0004 |
| 29.30 | −0.0088 | −0.0093 | −0.0087 | −0.0091 (−0.0092 to −0.0089, sd 0.9e-4) | −0.0091 (1.3%) | −0.0003 |
| 29.35 | −0.0087 | −0.0092 | −0.0085 | −0.0089 (−0.0090 to −0.0086, sd 1.3e-4) | −0.0089 (1.7%) | −0.0002 |
| 29.22 (descriptive) | −0.0065 | −0.0068 | −0.0066 | −0.0067 (−0.0068 to −0.0065, sd 1.1e-4) | −0.0067 (2.0%) | +0.0005 |

## Reading (exploratory)

- **The diagnostics agree.** Over the late interval, all five measures give the same leftward shift to within about
  6%. The 16–17 interior minima each move by the same amount to within about 1.5e-4. A pure translation explains
  98–99% of the change between the 15,000 and 30,000 snapshots. The width changes by about 2–5% of the shift.
- This supports **coherent translation of the whole patch**, not amplitude redistribution, and not just movement of the
  threshold crossings. It does not identify the cause. Nothing here distinguishes slow relaxation of the patch's sub-cell
  position (for example toward a lattice-symmetric placement), a mesh effect, or a genuinely travelling localized state.
  That needs new simulations under separately registered criteria.
- **Earlier intervals behave differently.** Between 1800 and 6000 the change is mostly not a translation (alignment
  explains almost none of it, unexplained part 98–100%), because the edges are still finishing their outward adjustment.
  In the 18-cell runs, the interior-minima mean also jumps about −0.013 between 1800 and 3600 (individual spread sd
  0.15) during that adjustment, then runs parallel to the other measures. The early part is shape change, not drift.
