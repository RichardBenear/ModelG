# MG11-G5 conditional stochastic pattern-response pilot — readiness (MG100b)

G4 established that the new finite-resource stochastic mechanism behaves correctly on the registered small systems. This next pilot will ask whether the geometry-dependent field response found in G3 is still measurable when the full lattice carries token noise. Preparation now supports launching that bounded test: checkpointed paths reproduce uninterrupted paths exactly, accounting and analysis checks pass, and the local cost has been measured. No long-horizon pilot result exists yet. Even a positive pilot would establish a conditional classical response from four prepared contexts, not particle lifetime, quantum entanglement or a physical interpretation of A/B.

Status: **READY FOR CONCRETE PI REGISTRATION; production not run**.
Specification: [draft1](../../spec/MG11G5_pilot.md). Source and prospective seed/donor plans were committed and pushed before each numerical development phase.

## Proposed run

| Item | Fixed choice |
|---|---|
| Prepared contexts | G2 C1 fine parents at T10, LOCAL/CROSS origin x H/S; all seven stores retained |
| Token conversion | Omega1000, round-to-even populations and original reference denominators |
| Intervention | B LOCAL/CROSS only; A retains prepared network; active ring500 unchanged |
| Numerical settings | max tau1/24000 and1/48000; independent streams; simultaneous G4 step |
| Stochastic sample | four balanced blocks/setting, eight branches/block:64branches |
| Controls | eight contexts x continuum/finite-count drift:16deterministic controls |
| Horizon | post0..15 = absolute10..25;61quarter-unit frames |
| Primary | late20..25 X fixed-template projection, H-subtracted B interaction, equal origin weights |
| Inference | n4/setting, nominal97.5% t intervals; fine-minus-coarse90% Welch agreement inside +/-.02 |
| Directional outcome | both positive primary intervals AND step agreement AND complete/audited full run; otherwise unresolved |
| Location | PI local PC, max8workers, one thread/library |
| Forecast | **16.234798core-hours;121.76minutes ideal8worker wall**, including explicit safety/startup/IO allowances |
| Requested cap | **24core-hours aggregate;3600seconds/job**;75% forecast admission reserve, .25second polling |
| Storage/retention | full144main/control/excluded-warm raw files, prefixes/failures/seeds/analysis publicly retained; ~530MB uncompressed final arrays |

Allow roughly **2–2.5hours** local wall time, subject to load and state-dependent rates; the short-horizon extrapolation is not a guarantee. Forecast58445.273998worker-seconds is below64800reserve. Largest job forecast1023.241seconds is below3600cap. Free disk at qualification302GB. A forecast stop drains active jobs; cap/failure stops preserve files. No automatic retries or replacements. Scientific closure requires later PI sign-off.

## Completed validation

- 48small fixture/setting/path segmented-versus-monolithic comparisons, plus16held full-lattice comparisons, all physical states, edge flows and step statistics bit-exact. No arrivals/products act twice within a step; the inherited G4 kernel is unchanged.
- 60integer raw audits passed: nonnegative stores, exact reaction/hop/edge reconstruction and closed total. All checkpoint counters and the absolute postclock carry forward.
- 16short deterministic controls passed accounting;16separate Radau population references passed the predeclared2e-6scaled bound. Largest scaled discrepancy is 8.75568313537e-11.
- Same actual integer parents and fixed references verified across destinations/settings/blocks. No homogeneous reduction is used for noisy H backgrounds.
- Known signed interaction/background cancellation, independent dense projection, Student-t CDF checks, complete/prefix/missing/empty raw pipelines and separate interval recomputation passed. Missing branches invalidate balanced blocks and never become zeros.
- Six adversarial cases rejected: four invalid stream inputs, one corrupt ledger and one altered primary result.
- Original development stopped before any draws/integrations on an excessively strict inverse-t comparison (3.93543e-11absolute discrepancy). Original source/plan/failure preserved. Version2 verifies both quantiles by CDF residual<1e-12; confidence levels, physical model and scientific criteria unchanged.
- Version2 development65.564962seconds; supplement10.795452seconds; including10seconds reserved for original failure, accounted86.360414seconds<900development ceiling. No late-window development data or full pilot generated.

## Interpretation and remaining risks

The comparison estimates future stochastic response conditional on four deterministic prepared contexts; it does not sample noisy parent formation. Fixed A topology differs between origins and can influence fluctuations, so origin differences cannot isolate memory alone. Control closures expose rounding/finite-count effects but are not exact stochastic ensemble means. All noisy H backgrounds evolve across500nodes.

Four blocks give weak precision. Nominal t/Welch coverage relies on distributional approximation, not G4's finite-sample CDF guarantee. No joint95% coverage is claimed across primary and step-agreement clauses. Independent noisy trajectories are compared statistically, never by pointwise tolerance. Refinement failure or overlapping intervals are valid unresolved outcomes, with no adaptive replication or favorable filtering.

Projection can persist or increase after reference matching is lost; post-switch continuity is descriptive and cannot satisfy a20-unit persistence definition within15units. OFF effects, parent uncertainty and additional graph families are outside this pilot. Source/qualification/donor hashes and the full seed plan must remain unchanged before registration/launch. Same-author separate recomputation shares donors/template and accounting dependencies and is not external replication.

## Execution after approval

Use src/run_mg11g5_ready.py --register with the PI approval text, commit/push registration, then --production with MODELG_WORKERS=8. This entry point requires both qualification layers. The inherited budget controller launches the16controls first, then64stochastic branches, automatically evaluates/audits, and adds stock/uptake/post-reference-continuity diagnostics. The lead reports the result and verifies the published full archive before requesting closure. External review is optional at Richard's discretion.
