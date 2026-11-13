# MG162B — independent real-trajectory candidate-instrument pilot

Status: concrete proposal awaiting PI approval. MG162A replay/interior challenge is supported, with pronounced near-boundary sensitivity. This packet does not declare boundary classification qualified. Freeze MG161's actual fitting/classification/calibration unchanged; no retraining or new PSD estimates using fresh trajectories. No Model G trajectories have been drawn for this stage.

## Fixed sample and metrics

Eight new independent seed pairs:4prepared and4mass-matched homogeneous, each comparing0.0001/0.00002 through post-preparationT=0.1, observed every0.001. Same frozen500-site/Omega10 finite-source state, ring, stock references and retained products. All earlier trajectories excluded from the new denominator. Source/seed manifest and actual approval published before draws. Default historical global measurements may be retained for comparison but do not replace the new primary instrument.

Primary outputs are the complete MG161 point labels, fitted parameters and per-prepared-path all101frame compatibility/first noncompatible bracket; for each homogeneous marginal path, any compatible crossing. A crossing in a control is a possible false detection or new structure needing investigation, not automatically a physical birth or a known instrument error. Both marginal paths within a seed pair are correlated. Four pairs per condition cannot establish small population error probabilities or long-lived particle survival.

Frozen engineering boundary advisory: among core-compatible frames, flag relative amplitude<=0.6 or fitted scale within0.125 of either0.5 or2. This was already encoded in the draft entry point before A's results were read. It is a sensitivity warning based on the estimator/grid, not a statistical confidence interval. Report all original point labels unchanged alongside the advisory. No rescore of MG160-P1 or MG161.

Declare clear-case pilot separation supported only if every prepared marginal path is core-compatible at every frame with no boundary advisories, and no homogeneous marginal path has any core-compatible frame. Otherwise report not established, retaining every outcome. This is a strict finite-pilot operational criterion, not a population reliability qualification. A wider physical-survival interpretation remains unresolved even if the finite criterion passes. No sample increase, outcome-dependent stopping, exclusion or changed thresholds. NoT=1/T=15 extension is included.

## Execution and retention

Reuse the frozen MG160 chunk writer/paired engine through an isolated-process directory/source-binding adapter. Copy the unchanged prepared/control inputs into ignored local_raw/MG162_fresh_model, with the public registration mirrored there for the writer. Continue RNG state through100segments; no products/arrivals/catalysts reused within their own tick. Keep the existing one-million-event segment cap and exact count/histogram/conservation checks. Retain full observations and two numeric sentinel traces locally. Audit all checkpoint/RNG/boundary records and independently recompute instrument outputs before a final scientific report; candidate scores alone are not an audited verdict. Record all analysis attempts/timings.

Measured prior cost36912.7worker-seconds for the identical-size mixture;25% allowance yields46140.9worker-seconds (12.8169core-hours), about1.60hours on8workers. Proposed ceiling100000worker-seconds (27.7778core-hours),12000/job, max8localworkers/library-thread1, plus1800activeparent/audit seconds. Use forecast admission safety and no scientific retries. Reserve2GB raw local storage, publish only compact decisions/source/figures/seed and hash records. Keep an immutable local-only raw commit and prevent archive pushes. Retain any incomplete denominator as incomplete.

PI approval of these concrete criteria and budget is required under RESEARCH_PROTOCOL.md workflow3 before registration/draws. The entry point requires fresh_approval.json; none exists. External review is optional and not requested. Earlier stage closures are not inferred.
