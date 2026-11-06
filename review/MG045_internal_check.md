# MG045 — MG10-C1 result check

2026-11-06. ChatGPT / Codex, implementation lead. No external review claimed.

All ten jobs ran after prospective registration commit 1451074 was pushed. The
source/spec hashes remained unchanged during production. No extra mode, amplitude,
trajectory, comparison window or detector threshold was selected after outcomes.
The planned primary linear scan selected mode32; all 250 responses are retained.

Registered outcomes: HOMOGENEOUS REFERENCE RESOLVED; CONTROLLED AMPLIFICATION
SUPPORTED; RESOLVED-MINIMUM QUALIFICATION PASS. Separate raw checking agrees on
all outcomes, coefficients, windows, detector counts, hashes and cost totals.
The original MG10 UNRESOLVED result is not changed. Closure remains PI decision.

Scientific challenge: the enormous formal linear gain is not a physical population;
only the prespecified small-amplitude intervals support nonlinear comparison.
The worst eligible relative error is .03295%, below5%. Roundoff-seeded O+ departs
from exact homogeneity, while the uniform Euler path follows the reduced ODE;
this supports the amplification explanation but supplies no physical noise law.
The diagnostic minima are not self-sustaining particles or a nucleation-rate
measurement. All 301 uniform Euler states have zero new detections, and all 23
prespecified measurement fixtures pass; that does not qualify every possible
depleted-background particle definition.

The checker shares the original smoothing/reference and prior independent ledger
helper, but reimplements periodic prominence by contour walks and does not import
the new model/evaluator. It is same-author recomputation, not independent external
replication. All nonlinear ledgers/oracles pass. No timeout or production cap stop.
The report preserves failed development attempts and physical/modeling limits.
