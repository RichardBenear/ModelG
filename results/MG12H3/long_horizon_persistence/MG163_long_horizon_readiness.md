# MG163 — long-horizon test preparation

MG162A/B are signed off. The next proposed test asks whether the recognizable local Model G core persists for ten times longer, through T=1, while finite A/B supplies continue to evolve. Keep the successful short-pilot measurement fixed and compare four new prepared-pattern pairs against four matched homogeneous-control pairs at both tick sizes. A positive result would extend evidence of sampled short-term persistence; it would not prove permanent particles or physical lifetimes. No new simulations have started.

The complete prospective design is [MG163_long_horizon_persistence.md](../../../spec/MG163_long_horizon_persistence.md).

| Item | Proposal |
|---|---|
| Horizon | T=1, ten times MG162B |
| Independent pairs | 4 prepared + 4 homogeneous controls |
| Tick sizes | 0.0001 and 0.00002 |
| Observations | 1001/path, every0.001 |
| Primary requirement | Every prepared frame compatible/no boundary advisory; no compatible control frame |
| Estimated runtime with25%allowance | About19.6hours on8workers |
| Hard worker budget | 800000seconds =222.222core-hours |
| Per-job ceiling | 100000seconds |
| Active preparation/audit/report allowance | 3600seconds |
| Full local data projection/reserve | About1.63GB /4GBreserve |
| Public data | Compact exports only; full raw local |
| Full event traces | Two fixed sentinel jobs, predetermined10%ofsegments |

Sparse trace retention changes verification coverage and storage, not the dynamics or scoring. Every segment still retains checkpoint, count, histogram and RNG records. Current measured free disk space is about308GB, well above the proposed8GB launch preflight. The runtime estimate is an extrapolation and may change as sources deplete.

Status: DESIGN READY FOR PI DECISION; execution implementation/qualification remains a prelaunch gate. No production seeds, allocation approval, or trajectories exist. The historical writer hardcodes100segments and a11980-second deadline; it must be parameterized in a new stage-specific writer and audited before launch. This report does not claim those implementation checks have passed. No previous scientific scores changed.
