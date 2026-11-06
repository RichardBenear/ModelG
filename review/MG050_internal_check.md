# MG050 — S2 prospective registration check

2026-11-06. ChatGPT / Codex, lead author; no external review claimed.

PI signed off S1 and authorized the next stage. S2 tests the actual h=.2 and
Omega100/1000 with full500-node H/S snapshot laws and short16-node exact-SSA
comparisons. It is not a long-horizon particle or physical-noise claim.

Five lattice/oracle tests and three apparatus tests pass. They cover unchanged-S1
kernel replay, rate scaling, Fourier normalization/phase/Nyquist, per-token covariance
enumeration, stock/ledger safety, lossless event packing and population-hash recovery,
retention after interruption, signed CDF/decision boundaries, synthetic families,
and actual owned-job cap/timeout/interruption/resume. The separate checker verifies
all32 timing classes, raw hashes, exact integer reconstruction, independently
derived conditional means and noise covariances, and canonical spectral measurements.

Before scientific registration, a static precision review found that the original
candidate R8192 spent more than half the .08 margin on the simultaneous-band
allowance. Draft1, its timings and calibration are retained. The new allocation
uses n>=2*log(4*M/.005)/.04^2, hence R16384 (next power of two) for both families.
This is a preproduction design correction, not a post-outcome sample extension.
No physical S2 ensemble has run. Margins, alpha, physical inputs, methods and
observable definitions are unchanged. Revised halfwidths are .03818084/.03705623.
Revised Bernoulli checks show coverage1.0 and qualification frequency1.0 for both
null and distance.01 families, zero qualifications at.15. These finite examples
check implementation; the confidence statement still relies on DKW/PRNG assumptions.

The 32 unchanged class timing measurements are rescaled arithmetically to the
revised allocation:3062.691866 worker-seconds (0.8507477 core-hours), below5400
allowed after reserve from7200; maximum projected job83.8231 seconds versus600.
Total qualification/timing/check/revision cost102.7158449 seconds, below900.
The fixed denominator is88 jobs,32 logical groups,16384 trajectories/group,
524288 trajectories and1280 unique production chunk streams.

Actual monitoring gets a prospective .25-second scheduled target with recorded
maximum observed gap. The generic controller change is administrative only;
historical samplers/results remain fixed in their commits. Source/spec hashes
are LF-normalized text hashes to avoid Git EOL artifacts; binaries use raw hashes.
All population states are reconstructible from retained event shards, checked
against hashes of the actual prepacking populations. No lossy packing or omitted
decision data is authorized. Durable space and seed uniqueness are verified
before registration; each shard remains below100MiB.

Adversarial limits: an exact snapshot-law match tests implementation, not the
accuracy of that law against continuous time. The short SSA arm supplies a separate
finite-interval check. Its primary step is the finest, not an automatic endorsement
of the coarser production step. Mode powers/cross-products test selected joint
spatial observables, not the entire joint field law. The snapshot centering uses
the known conditional mean and does not rescore MG10's uncentered recovery test.
Short-time success does not establish long-lived particles, rare-event rates,
geometry-derived supply or quantum entanglement. Those remain later studies.

Both analysis paths deliberately share the canonical FFT/scaling/rounding helper
and original seed helper. The separate checker independently derives initial
states, means, integer replay, population hashes, covariance enumeration, CDF
distances and labels without importing the new sampler/evaluator. Same-author
recomputation is not external replication. PI closure follows results.
