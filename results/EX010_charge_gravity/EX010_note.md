# EX010 (EXPLORATORY): reciprocal charge chemistry, and gravity plus charge in one particle

**Answer.**
1. **Charge can be made reciprocal by tuning.** Adjusting the "+" source's two rate constants makes it an exact mirror
   of the "−" sink: Q₊ = −Q₋ = 2.2438 and μ₊ = −μ₋ = 0.0086666, matched to 1e-13. With that, the implied pair speeds
   are equal and opposite. Like charges repel and unlike attract, with the same strength for ++, −− and +−. This is a
   match at one operating point, not a built-in symmetry.
2. **Gravity and charge coexist, and gravity is very nearly universal.**
   - With the G sink added, +, − and neutral particles all consume G (mass M > 0) and all drift toward lower G.
   - Their G-mobilities agree within 2% (kg = 0.005) and 1.2% (kg = 0.05 and 0.5), so every type attracts every
     other gravitationally.
   - Gravity slightly breaks the charge mirror: Q₊ and |Q₋| differ by 0.2%, 1.9% and 10% at kg = 0.005, 0.05 and 0.5.
   - The gravity/charge force ratio is about 1e-2 at kg = 0.005 and falls as kg².

Exploratory, our construction, radial 3D with single-particle linear response. No direct pair runs here. Not a claim.

## Set-up
- I wrote a generic mass-action solver (`src/ex010_network.py`). It reproduces EX008 and EX009 to 5e-4 or better.
- **Species:** X, Y, W (core-only), C (conserved charge carrier) and G (conserved gravity carrier). C and G are fixed
  at 1 at the outer boundary, R = 60, and have no ether reactions.
- **Reactions:**
  - Base: the Brusselator (A = 140/9, b = 29) plus W's switch (W → ∅, Y + 2W → Y + 3W, 3W → 2W).
  - "+": Y + W ⇌ C + W, forward rate kf, backward rate kf·ε.
  - "−": C + X + W → W, rate 0.005.
  - Gravity, added to every type: G + W → W + X, rate kg.

## 1. Tuning the + particle (kg = 0)
| Particle | Rates | Q | μ_C |
| --- | --- | ---: | ---: |
| − | km = 0.005 | −2.24377 | +0.0086666 |
| + (before, EX009) | kf = 0.05, ε = 0.1 | +5.883 | −0.000622 |
| + (tuned) | kf = 0.053588, ε = 1.37951 | +2.24377 | −0.0086666 |

- **Why ε matters.** A source's C-mobility comes only from its backward step, C + W → Y + W. Raising ε strengthens
  that step while cutting net production, so |μ|/Q can be set anywhere. ε passes 2.4 before Q changes sign.
- **Pair speeds** from the C fields alone: |v| = μQ/(4πd²) = 1.55e-3/d² for each particle in every pairing.

## 2. Gravity plus charge
| kg | Type | M (G consumed) | Q | μ_G | μ_C |
| ---: | --- | ---: | ---: | ---: | ---: |
| 0.005 | + | 0.2300 | +2.2291 | −7.40e-4 | −8.654e-3 |
| 0.005 | − | 0.2284 | −2.2336 | −7.31e-4 | +8.661e-3 |
| 0.005 | neutral | 0.2171 | 0 | −7.25e-4 | 0 |
| 0.05 | + | 1.992 | +2.113 | −7.14e-3 | −8.548e-3 |
| 0.05 | − | 1.981 | −2.153 | −7.05e-3 | +8.611e-3 |
| 0.05 | neutral | 1.886 | 0 | −7.01e-3 | 0 |
| 0.5 | + | 9.008 | +1.581 | −5.59e-2 | −7.93e-3 |
| 0.5 | − | 8.998 | −1.753 | −5.53e-2 | +8.26e-3 |
| 0.5 | neutral | 8.597 | 0 | −5.54e-2 | 0 |

- Charged particles consume about 6% more G than neutral ones: the extra core reactions slightly enlarge the core.
- **Gravity vs charge.** At kg = 0.005, |μ_G·M| ≈ 1.7e-4 against |μ_C·Q| ≈ 1.9e-2, a ratio of about 9e-3. In EX009
  both M and μ_G scaled as kg, so the ratio should fall as kg². Gravity at 1e-36 of charge would need kg of about
  5e-20 in these units. That is an extrapolation, not a run.

## Limits and open points
- The reciprocity is tuned, and gravity spoils it at the percent level. A built-in mirror symmetry is still missing.
- Single-particle linear response only. A direct three-way check (+, −, neutral pairs in the axisymmetric solver from
  EX009) is the obvious next run.
- Particle "mass" here is the consumption rate M. There is no inertia, since the dynamics are overdamped: velocity is
  proportional to force.
- About 0.1 core-hours.

## Files
- Code: `src/ex010_network.py`.
- Records: `part1_tuned.json`, `part2_combined.json`. The tuning and combined scans ran from inline drivers; their
  commands are logged in `part_runs.md`.
