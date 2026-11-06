# MG10-S2 — Actual-lattice spatial-noise qualification

MG049, 2026-11-06. **Prospective design; no scientific S2 ensemble has run.**
PI signed off S1 and instructed proceeding. S1 is CLOSED with its limited small-
system qualification. Original MG10 remains UNRESOLVED; C1 is CLOSED. This study
tests h=.2 and Omega100/1000 before any long-horizon stochastic particle claims.
No external review was requested with the instruction to proceed.

## Model and fixed inputs

Use the unchanged MG10-S1 finite-reactant leap and independent SSA kernels. Their
diffusion argument is a per-neighbor hop rate; pass (1,1,12)/h^2=(25,25,300).
P17, finite reservoirs C_A=C_B=1,T_ref300, no external input. At each site
A0=round(Omega*14*300),B0=round(Omega*29*(140/9)*300). Counts G/X/Y are rounded
Omega*(S0+phi), with phi=0 for H and the committed MG003 reference for S. All
stocks start full; Z=W=0. Keep actual denominators and all initial arrays.

Two complementary arms:

1. **Full snapshot law.** N500, H/S, Omega100/1000, one leap at each maximum
   delta in {1/6000,1/12000,1/24000}. Compare with an independently drawn exact
   snapshot-law oracle, not a continuous-time SSA. At these states B>=X, so the
   joint R3/remaining-X construction has exactly the ordinary X-multinomial law.
   The oracle uses NumPy Generator.multinomial for G/X/Y and bounded binomial A,
   without invoking S1 drawing helpers. This tests implementation of the finite-
   step law and spatial noise on the actual500-node states, not its CTMC accuracy.
2. **Short continuous-time comparison.** N16 periodic ring, H start only, the
   same h/densities/stocks, checkpoints {1/8000,1/4000,1/2000}. Exact SSA versus
   the same three leap maxima. Finest versus SSA is primary; coarse/middle
   comparisons are descriptive. The short ring resolves wavelength3.2 near the
   relevant full-ring wavelength3.125, plus its high-frequency end. It does not
   establish long-horizon, full-ring nonlinear accuracy.

Each of the24 snapshot groups and8 short groups has R8192 trajectories. Snapshot
groups run as24 jobs, chunks512; short groups have8 fixed blocks of1024 each,
chunks256, giving64 short jobs. Total88 jobs,32 logical groups,262144 trajectories.
No replacement, extension or data-selected arm. Partitioning avoids a long exact-
SSA job while preserving fixed streams and the full sample denominator.

## Spectral observables and reference moments

For each trajectory define z_s(m)=DFT_m(n_s-mu_s)/sqrt(N*Omega*time).
For snapshot mu is the exact conditional Euler mean from its fixed count state
(including falling-factorial R4 and lattice drift), not a fitted sample mean.
For short H trajectories subtract the uniform initial count; every nonzero mode
has zero expected amplitude by translation symmetry, irrespective of mean drift.
Use NumPy real FFT with Nyquist imaginary component exactly zero. A dedicated
canonical measurement helper, shared by the two analysis paths, freezes scaling
and rounding of observables to8 decimal places to suppress irrelevant FFT ties.
This numerical quantization is far below the physical noise scale; it is not a
noise cutoff or a particle detector. Preserve the raw event/state data to recompute it.

Probe modes {32,250} on N500 and {1,8} on N16. At each mode retain8 observables:
|zG|^2,|zX|^2,|zY|^2,Re(zG*conj(zX)),Re(zG*conj(zY)),Re(zX*conj(zY)),Re(zX),Im(zX).
Power and cross-product distributions probe joint spatial fluctuations; they do
not determine the entire field's joint law or every mode/cross-mode covariance.
Also derive the exact snapshot-law same-mode covariance analytically from the
per-token multinomial covariances and finite-A binomial variance. Its comparison
to the CTMC instantaneous noise matrix is descriptive; finite-interval SSA
covariance is not identical to that instantaneous matrix.

## Frozen statistical decisions

Use independent seeds and empirical-CDF sup distances, including signed continuous
observables. Two simultaneous families, with alpha=.005 each (global union <=.01):
- Snapshot: 2 densities*2 starts*3 deltas*2 modes*8 observables=192 comparisons,
  leap versus snapshot oracle for every delta.
- Short: 2 densities*3 checkpoints*2 modes*8 observables=96 comparisons, finest
  leap versus SSA. Other steps are descriptive, not additional equivalence claims.

For a family of M comparisons set w=2*sqrt(log(4*M/.005)/(2*8192)). Each distance
has simultaneous interval [max(0,d_hat-w),min(1,d_hat+w)]. Require every upper
bound<=.08 in that family. The8% CDF margin and R are frozen here; they are not
concentration, variance-percentage or particle-rate margins. This larger margin
than S1 reflects the cost of the high-density exact oracle and is explicit.
Pseudorandom trajectories are assumed independent; dependence between observables
and shared reference samples is handled conservatively by the union bound.

Separate family labels: QUALIFIED ON TESTED SNAPSHOT LAWS / QUALIFIED ON TESTED
SHORT-TIME SPATIAL LAWS if complete and every corresponding upper<=.08;
NOT QUALIFIED ON TESTED CASES if complete with a lower>.08; otherwise UNRESOLVED.
Invalid transitions/probabilities/ledgers, oracle failure or missing mandatory
qualification gives QUALIFICATION FAILED. Missing/terminated jobs retain their
denominator and cannot support a qualification claim. Overall is QUALIFIED WITHIN
REGISTERED S2 SCOPE only if both families qualify and all exact checks pass.
No physical particle or nucleation-rate verdict belongs to S2.

## Streams, qualification, run limits

Use seed_for(10,case,method,stage,global_chunk). Snapshot cases0..3 are
H100,S100,H1000,S1000; method=2*delta_index+sampler (leap0,oracle1). Short cases4/5
are H100/H1000; methodsSSA0 and three leap levels1/2/3. Stages4901 development,
4902 timing,4903 production; block boundaries cannot change global chunk indices.
Freeze/check all640 production chunk seeds for uniqueness before registration.
The oracle uses PCG64; S1 kernels retain their seeded Numba RNG. Record all streams.

Before production: unchanged-S1 regression; h scaling; no safety adaptation in the
registered snapshot states; independent snapshot oracle/ledger tests; actual-vs-
reconstructed population hashes; FFT/phase/Nyquist/rounding tests; analytic noise
matrix versus explicit per-token Fourier-event enumeration; statistical label and
synthetic-family checks; and actual interruption/timeout/budget adapter checks.
Keep source/spec/reference hashes and package versions. A separate checker reads
raw data without importing sampler/evaluator; shared seed/measurement dependencies
are declared, and mathematical conditional means are independently verified.

Lead-allocated local ceilings under PI's proceed instruction: qualification/timing
900 worker-seconds; production7200 (2 core-hours), per-job600, <=8 workers with
library threads1. Use the qualified owned-process pool; poll<=.5 seconds, include
in-flight elapsed time, kill only owned jobs, preserve partial evidence. No automatic
extension. Time every32 class with a charged one-trajectory warm-up then64 snapshot
or8 short trajectories. Forecast twice measured compute scaling plus conservative
startup/checkpoint I/O overhead for every fixed job; require total<=5400 and each
job<=600 before prospective registration/production. Otherwise stop/report.

## Raw retention and closure

Use persistent local results/MG10S2/. Snapshot shards retain all13 reaction/directed-
hop event arrays as uint32, initial state, stocks, exact conditional mean and seed.
These events reconstruct every final population; record a SHA-256 of the actual
int64 population array before packing and independently verify reconstruction.
Every event count is bounded by its initial source count (<2^32 here); assert the
range before conversion, never truncate. Snapshot shards512*13*500*4=13312000
bytes before compression. Short shards retain full seven-population/seven-reaction/
six-hop arrays (256*3*20*16*8=1966080 bytes). No blob>=100MiB. Keep every shard,
metadata, costs, initial arrays, statistics, environment and hashes; commit/push.
Checkpoint per completed shard; do not silently rerun an interrupted block.

PI closure follows results. No historical result is rescored. Long-horizon particle
tracking/persistence and physical fluctuation or reservoir derivation still require
separate definitions and evidence; a short-time noise pass does not establish them.

Source/spec integrity in S2 uses SHA-256 of UTF-8 text normalized to LF, so Git line-ending conversions do not masquerade as scientific code changes. Binary references and raw files retain byte-exact SHA-256. Older study registrations are preserved unchanged.

Apparatus checks before production:256 Bernoulli families at true CDF distances0,.01,.15 for each family size192/96, with the fixed R/alpha/margin. Require observed coverage>=.98, null qualification frequency>=.70, and no qualifications at.15; inside-margin power is descriptive. This finite check validates implementation, not universal coverage or a basis to retune the design. The conservative simultaneous bands can leave even identical distributions unresolved; that outcome is retained.

Production monitoring uses a .25-second target period, allowing administration within the .5-second requirement. The shared controller now schedules sleep after subtracting administration time and records the maximum observed polling gap. This prospective monitor-only change leaves the historical physical samplers and results untouched; old controller versions remain in their registration commits. Timing used the earlier, slower monitor and remains conservative.
