# MG146 — Contrast estimator and compact archive validation

The new estimator passed its bounded implementation checks without changing the token dynamics. For the same seeds, it reproduced every event, state and final random-generator state from the frozen sparse engine while additionally recording a predictable-drift estimate of mean pattern contrast. This provides a promising route to measuring the weak Gpattern more precisely. It is not yet evidence of a particular precision gain on the actual ring or a qualification of particle preservation.

## Method and raw inputs

For a fixed linear contrast, the expected terminal value equals its initial value plus the expected time-integrated accepted-event drift. The implementation groups channels with identical available-count monomials and records each monomial's time exposure. Multiplying those exposures by fixed contrast coefficients reconstructs the estimator. Only affected groups are updated; the final no-event interval is included. The sampler, reservations, rate tree, products and random calls are unchanged.

Exposure arrays are retained as raw decision data, not merely the final estimator. They add payload beyond earlier count-only compact prototypes. Every path still retains count/histogram/checkpoint provenance; predeclared sentinels retain full order/times. A drift-integral estimate is not equal to a direct terminal contrast on every path, and unbiasedness does not guarantee lower variance in every regime.

## Completed development checks

-128same-seed mixed test paths,384ticks: all event IDs/marks, states and RNGstates exactly matched the frozen sampler.
-6scripted checks: known one-channel integral, boundary exclusion/truncation, empty case, zero weights, exposure corruption and histogram corruption.
-4096one-channel paths: exact expected output1-exp(-0.05)=0.0487706. Direct mean0.0493164, integrated mean0.0487192; both satisfy the frozen6MCSE+1e-6diagnostic.
-4mixed fixtures x512paths: direct-minus-integrated means were consistent with zero under the same diagnostic, with fresh full-rate drift and ordered-trace integration checks.
-Independent saved-array/einsum recomputation passed for the stored means and exposures. Same-author work and shared coefficients are disclosed; no external replication is claimed.

The one-channel example had an observed variance ratio of **1138.93** (direct/integrated). This is specific to that solvable short-horizon case and must not be transferred to actual G/X/Yprecision. Empirical diagnostics are implementation evidence, not mathematical proof of universal unbiasedness or formal confidence coverage; the conditional-intensity identity supplies the mathematical basis.

The batch took **3.03126seconds** within60seconds on one worker/librarythread1. Source/seeds were published as`a5caaca` before draws. No failures or retries occurred. Frozen sampler/source hashes remain unchanged. Development draws are excluded from any future pilot/qualification denominator.

Next: actual-ring compact writer/estimator preparation and measured24pathprecision-pilot budget. Five-percent initialpatterncontrast target and approvedcompactstorage choice remain unchanged. No actual-ring precision-pilot production has started.
