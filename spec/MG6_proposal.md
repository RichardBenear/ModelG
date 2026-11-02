# MG6 — G-coupled spatial growth (Phase 5): design proposal, DRAFT 1

*Stage MG019. Claude, 2026-11-02. Status: **DRAFT, not registered, for PI decisions and ChatGPT review.** No MG6 run of any
kind has been made. Every threshold below marked **[proposed]** is mine and needs PI sign-off before it is frozen. Items
marked **[PI decision]** are open choices, not defaults.*

*Computation behind this draft (disclosed): one descriptive calculation on MG4c's stored b = 28 start states, with no
integration and no MG6 dynamics: `spec/MG6_design_calc.py` (§5). It reads φ_G and node positions only.*

*A correction to my own earlier wording: in MG5 §1 and in my last message I described MG6 as testing whether a particle
"moves down a G gradient". The registered roadmap (`InitialPlan.md`, Phase 5) is different: **measure the volume-growth
(curvature) deficit around a particle, its falloff exponent, and whether it scales with the integrated G-deficit**, with
a decoupled run and a randomized-coupling run as null controls. This proposal follows the roadmap. Particle motion
under coupling is recorded as a secondary descriptive quantity only.*

## 1. Question
MG5 showed (Gate PASS) that a settled particle survives continual field-blind rewriting of the Delaunay graph under it
(single in 4/4 runs at ν = 0.002, 0.02, 0.2), and is mobile at ν = 0.2. MG6 adds the coupling: **the local G value
modulates the spatial rule**. Then:

1. **Deficit:** does a settled particle carry a measurable deficit (or excess) of space, here counted as nodes, around it,
   compared with a decoupled run and a randomized-coupling run?
2. **Shape:** how does the deficit fall off with distance from the particle? Does it follow the φ_G profile?
3. **Scaling:** does the deficit scale with the integrated G-deficit (§4, which depends on the definition of "mass", I2)?
4. **Feedback (non-trivial part):** does the coupling change the particle: its survival, depth, size or motion?

## 2. What MG6 builds on
- MG5 (rule v1.2 events S1/S2, tick-synchronous application, interpolating transfer, ROS2 with dt = 0.05, panel-calibrated
  h, MG4c start states at b = 28, MG5's detection, classes, balance and jump accounting). Everything not listed in §3 is
  unchanged from MG5.
- MG5 §3.3 already said the decoupled rule is the natural baseline, with birth rate νρ becoming νρ f(φ_G).

## 3. The coupling (rules amendment v1.3, MG6 only)
**What "volume" means here.** The operator L = (A − D)/h² uses a fixed edge length h per run. A region holding more nodes
therefore represents more physical area (area ∝ node count × h²). So a **node-count excess or deficit is the
volume-growth excess or deficit**. This is a modelling choice (the plan's "volume growth"), stated so the reader can
challenge it.

**Coupled births (thinning).** Births occur at rate ν ρ f(φ_G(x)) per unit area at position x; deaths stay at ν per node.
With uniform deaths the local stationary density is ρ f(φ_G) when φ_G varies slowly compared with the relaxation scale
1/ν.

Implementation that keeps events reproducible: candidate births are generated at rate ν ρ f_max from a `seed_for(6, …)`
stream, each with a uniform number u. At the tick boundary a candidate is accepted if u < f(φ_G(x)) / f_max, with φ_G(x)
the barycentric interpolation on the pre-edit mesh (the same interpolation §3.4 of MG5 uses for the transferred value).
Deaths are the MG5 stream. Accept decisions depend on the field, so the **committed replay specification covers the
candidate list; the accepted set is recorded per tick and hashed.** A candidate list is still a pure function of the seed.

**Coupling function [PI decision].**
| Option | f(φ_G) | Notes |
|---|---|---|
| **F1 (recommended)** | exp(κ φ_G) | Positive everywhere, no clipping. Density deficit where φ_G < 0. f_max = exp(κ·φ_max), with φ_max from the start states plus a registered guard (a tick with φ_G > φ_max is a flag). |
| F2 | max(0, 1 + κ φ_G) | Linear, so ∫(f − 1) = κ ∫φ_G = 0 exactly: no net volume change, a strictly local deficit with a compensating halo. Clipping is never active at the proposed κ (§5). |
| F3 | exp(κ min(φ_G, 0)) | Well only; the halo is ignored. f ≤ 1, so f_max = 1 and the accepted stream is a **subset of MG5's stream** with the same seed (see §6 on reusing MG5 as the decoupled arm). Net volume deficit is non-zero, which sidesteps the ∫φ_G = 0 issue below. |

**The ∫φ_G = 0 issue (assumption M6, I2).** With conservative diffusion the net ∫φ_G of a stationary particle is zero, so
a coupling linear in φ_G gives zero *net* volume change. The deficit is local: in the stored states φ_G runs from about
−0.28 (core) to about +0.05 (halo). Whether "integrated G-deficit" means the integral over the well only (φ_G < 0), the
well depth, or something else is exactly the open definition of mass (I2). **The scaling test in §4 cannot be written
until the PI picks this.** Options are in §4.

**Strength κ [PI decision].** Proposed panel κ ∈ {2, 4} for the main arms (§5 gives the expected signals). κ = 0 is the
decoupled arm.

## 4. Measurements and outcomes (all **[proposed]**, to be frozen after PI decisions)
**Deficit profile.** At each check (every 10 time units after settling), around the tracked particle centre c(t), count
nodes in annuli. Excess ΔN(R) = N(<R) − ρπR² (torus-wrapped, with the mean density measured on the whole torus at that
check). Average ΔN(R) over checks in the measurement window and over runs; report the mean and its standard error from
the spread across runs. Radii: R ∈ {1.0, 1.5, 2.0, 3.0, 5.0, 8.0}.

**Arms and null controls (per the roadmap).**
- **Coupled:** f as chosen, κ ∈ {2, 4}.
- **Decoupled:** κ = 0 (MG5 rule).
- **Randomized coupling:** the same f, but evaluated on φ_G taken from a *spatially unrelated place*. Two variants **[PI
  decision]**: (a) *shifted*: φ_G read at x + s with a random fixed offset s of at least 6 from the particle (keeps spatial
  structure, breaks the link to the particle's position); (b) *shuffled*: node-wise permutation of φ_G values (destroys
  structure, keeps the value distribution). I recommend (a) as the main null and (b) as a check.

**Registered outcomes [proposed].**
- **INVALID-NULL:** the decoupled or randomized arm shows |ΔN(1.5)| / se ≥ 3 at the particle. Then MG6 stops for the PI.
- **DEFICIT PRESENT:** in the coupled arm ΔN(1.5) < 0 with |ΔN(1.5)| / se ≥ 3 at both κ, and it differs from both nulls by
  ≥ 3 combined se.
- **NO DEFICIT:** otherwise.
- **Falloff (descriptive, not pass/fail):** compare ΔN(R) with the *direct prediction* ρ·∑_{nodes in R} (f(φ_G) − 1)·dA
  computed from the run's own measured φ_G. Report the ratio observed/predicted per radius. Also fit |ΔN(R)| or the local
  density excess against r with (i) an exponential and (ii) a power law over the stated range and report both with
  uncertainty. **I expect an exponential-like decay** (the φ_G profile of a stationary Turing-type soliton is short-range,
  with λ ≈ 3.1), **not a 1/r-type power law.** Gate MG6 should not presuppose a power law; if the PI wants one tested,
  that is a separate pre-registered hypothesis.
- **Scaling [definition awaiting PI, I2]:** the deficit amplitude in the κ = 4 arm over the κ = 2 arm is predicted to be 2
  for F1/F2 at small κφ_G. That is a check on the *coupling*, not on the *particle*. To test scaling with the particle's own
  integrated G-deficit, the arms must vary the particle (for instance b = 27.5, 28, 28.5 with each b's own stored state;
  MG4c has these). **[PI decision]**: include this b-scan (3 b × arms) or defer it.
- **Feedback (descriptive):** class (single / multiplied / dissolved), depth against the decoupled arm at matched times,
  and mobility with D_p, exactly as MG5 measured them.

**Built-in expectation, stated now (important for honesty).** With fixed h and a rule that sets the local density to
ρ f(φ_G), the first-order deficit is **close to built in** by the coupling. The positive result "a deficit is present"
is therefore a *validation that the implemented rule realises the intended coupling* (and that the particle's profile is
stable under it) rather than a physical discovery. The informative parts are the feedback on the particle, whether
the observed deficit matches the direct prediction (if the particle shifts or reshapes, it will not), the randomized-null
contrast, and the falloff form. The write-up must say this.

## 5. Signal against noise (design calculation, `spec/MG6_design_calc.py`)
Starting states, b = 28, four runs; ρ = 47.25 per unit area. Expected node-count excess within R if density were
ρ·exp(κφ_G), against the Poisson counting sd of one snapshot:

| κ | R = 1.5 (sd 18.3) | R = 3 (sd 36.6) | R = 5 (sd 60.9) | R = 8 (sd 97.5) |
|---|---|---|---|---|
| 1 | −20 to −24 | +1 to +2 | 0 to +0.5 | +2 |
| 2 | −37 to −45 | +6 to +8 | +4.5 to +5 | +8 to +9 |
| 4 | −63 to −77 | +29 to +31 | +23 to +26 | +29 to +35 |

- The deficit is concentrated within about 1.5 of the centre and **is essentially cancelled by R = 3**; at κ = 4 the halo
  overshoots to a net positive excess (the convexity of exp).
- One snapshot at κ = 2 is about 2σ at R = 1.5. With 4 runs × about 10 independent checks, the standard error falls to
  about sd/√40 ≈ 3 nodes, i.e. a 12σ signal if the density follows its target. The ≥ 3σ criterion has margin.
- The signal-to-noise at larger R is poor (excess ≪ sd), so the falloff beyond R ≈ 3 will be hard to measure at these κ.
  A larger κ buys signal at the cost of leaving the small-κφ_G regime; **κ = 4 gives φ_G in the exponent up to about
  −1.1** (f ≈ 0.33 at the core), which is a strong perturbation of the density at the particle, and feedback on the
  particle must be watched.
- This calculation assumes the local density has equilibrated. The density relaxation time is 1/ν (5 time units at
  ν = 0.2; 50 at ν = 0.02). A moving particle also drags its target density; at ν = 0.2 its RMS displacement over 100 units is
  about 0.6 (MSD = 4 D_p τ with D_p ~ 1e-3), small against its core size and slow against the relaxation time.

## 6. Runs, cost and budget [proposed]
- **Rate:** ν = 0.2 only for the main arms (relaxation 1/ν = 5, fastest equilibration; MG5's ν = 0.2 particle survived 4/4).
  ν = 0.02 as a secondary arm if budget allows **[PI decision]**.
- **Start states and horizon:** the four MG4c b = 28 states, settle then measure. Proposed t_end = 200, measurement window
  t ∈ [100, 200] (20 relaxation times at ν = 0.2).
- **Arms at b = 28, ν = 0.2:** coupled κ = 2, coupled κ = 4, randomized (shifted) at κ = 4, decoupled κ = 0. 4 arms × 4
  runs = 16 runs. Optional: shuffled null (+4), κ = 2 shifted (+4).
- **Reusing MG5 as the decoupled arm.** MG5's ν = 0.2 b = 28 runs share the starts and seeds. They are a valid decoupled arm
  *only if* the MG6 candidate list reproduces MG5's event stream when κ = 0, which requires f_max = 1 (option F3) or a
  design where κ = 0 reruns the MG5 stream exactly. MG5 ran to t = 500; MG6 would compare at matched times. **[PI
  decision]**: reuse (saves 4 runs ≈ 6 core-h at t = 200) or rerun for a clean paired control (recommended: a rerun is
  cheap enough and removes the cross-study dependence).
- **Cost model (measured in MG5, ledger):** ν = 0.2 at t = 500 cost 15.08 core-h for 4 runs (3.77 core-h each, about 27 s per
  time unit, the same as the Part 0 rate once the relaunch used 5 workers). At t = 200 and the MG5 ledger rate: about 1.5
  core-h per run. **16 runs ≈ 24 core-h; with 4 optional runs ≈ 30 core-h.** The coupling adds a per-tick interpolation
  of φ_G for accept decisions (small against the ROS2 cost). I will time it before budgeting (Part 0, §7).
- **Run cap and cutoff:** register run caps (the programme's convention): proposed cap 24 runs, 40 core-h launch cutoff on
  the ledger. Workers: 5, one BLAS thread each (MG017). **Order: decoupled and randomized arms first**, so the
  null controls always run before the cutoff can bind. (The MG016 ordering error cost 65 core-h and must not repeat.)

## 7. Part 0 (before any registered run, on an unregistered graph) [proposed]
1. **Accept-rate check:** with κφ_G fixed to a known synthetic field, the accepted-birth fraction equals mean f/f_max to within
   statistical error; the accepted candidate set is a deterministic function of the seed and field.
2. **κ = 0 equivalence:** κ = 0 reproduces MG5's event stream and class on an MG5 Part 0 check run (hash-identical lists,
   same counts). This anchors MG6's decoupled arm to the registered MG5 rule.
3. **Tick check:** Δ = 0.2 against 0.1 at κ = 4: counts and the centre/depth criteria as in MG5, plus ΔN(1.5) agreeing within
   its standard error.
4. **Density equilibration:** on a settled start with a synthetic frozen φ_G, local density reaches ρ f within 5 relaxation times.
5. **Flag guard:** no tick with φ_G > φ_max; jump accounting and interpolation audit as MG5.

## 8. Not in MG6
Token noise; non-geometric (combinatorial) rewriting; 3D; a weighted-operator arm; coupling of deaths; making h depend on
local density (the operator keeps its fixed h, §3); gravitational *motion* as a gate (descriptive only); a multiway/frame
test (MG7).

## 9. Decisions needed from the PI
1. **Coupling function:** F1 (recommended), F2 or F3.
2. **κ panel:** {2, 4} proposed.
3. **Randomized null:** shifted as the main null (recommended), shuffled optional.
4. **Definition of integrated G-deficit / mass (I2)** and whether to include the b-scan (3 b) now or defer scaling.
5. **Decoupled arm:** rerun (recommended) or reuse MG5.
6. **Rates:** ν = 0.2 only, or add ν = 0.02.
7. **Thresholds:** the 3σ rule, the radii, and the t_end = 200 / window 100–200.
8. **Budget:** run cap 24, 40 core-h cutoff.
