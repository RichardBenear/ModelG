# MG055 — D1 source-law discrepancy and prospective correction

D1 completed80/80 paths with exact accounting; its frozen verdict is not rescored. Posthoc autonomous A-source check found HIGH fine-step deficits0.300-0.304%. Original report/raw archive and caveat are published separately. PI acknowledged status and authorized next steps, not a claim that the bias was physically acceptable.

## Isolated numerical mechanism

The installed legacy binomial initializes P(0)=(1-p)**n. At n=4200000000000, p=6.9444444444444436e-14, its compiled initial probability is0.74719166734931586, versus the stable log1p-based value0.74701750031042502. This difference is present before drawing any trajectories. The geometric-spacing source avoids both forming1-p and normalizing an inaccurately initialized binomial recurrence. It preserves a finite binomial support; it is not a Poisson approximation.

The independent80-digit finite-binomial CDF agrees with the stable analytic law. Both the original Numba and NumPy Generator fine-case means show deficits relative to the exact mean in this diagnostic; NumPy is therefore not a reliable sole oracle at these enormous stocks/tiny probabilities. The correction's mean agrees within the observed scalar sampling uncertainty. The registered CDF margin.0008 was met by all6 primary comparisons, including legacy: the test is too loose to itself reject this small legacy discrepancy. Do not relabel legacy CDF tests as FAIL or claim a statistically resolved DKW lower bound; no primary legacy lower bound exceeded0. The direct erroneous P(0) computation plus D1's long accumulated source-count discrepancy supplies the numerical mechanism/evidence. All outcomes and thresholds are preserved.

## Qualification and forward scope

Revised2 primary source laws and4 supplemental fixtures meet their frozen bounds; all histograms, raw moments and overflow values retained. Twelve low-stock/non-HIGH same-seed regressions are bit-identical to frozen D1; four corrected CHEMO/HIGH chunk-vs-stream regressions and explicit integer ledgers pass. Bounds/domain failures reject. New32/32 and31/32 synthetic reporting fixtures and192 independent summaries pass. These are qualification data, not new physical outcomes or trajectory-law equivalence.

Only new D2 files use the corrected helper. Frozen D1/S1/S2 sources and installed numerical libraries are untouched. New32-path controlled comparison uses fresh CHEMO/HIGH streams; no replacing D1 observations or pooling denominators. Capacity1 D1 behavior remains descriptive context; no blanket particle/persistence claim. Existing 60-core-hour authority is sufficient:33.333427 used,900seconds qualification ceiling and24-core-hour corrected ceiling total57.583427 worst case. No increase requested or assumed.
