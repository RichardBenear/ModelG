# MG128 — H3 efficiency development

The first step toward a faster strict-tick engine passed a deterministic check: count-based event channels reproduce both the changes to populations and the participants that must be reserved in every one of the 896 saved small-system states. This supports replacing enormous lists of individual token combinations with population-based calculations, while keeping the requirement that a token can participate only once per tick. It does not yet establish runtime improvement or a stochastic particle.

## Completed check

All **896 states** and **4442 count channels** passed comparisons to explicit identity enumeration. Joint output increments and participant reservations matched; rate comparisons used relative and absolute tolerance 1e-12. In particular, R4 reserves both X catalysts even though their species do not change. All tested updates conserve token count, subtract reservations from the available snapshot and apply products only to the separate output. No random paths or integrations were run.

The [development specification](../../../spec/MG12H3_efficiency_development.md) gives the exponential-clock pruning argument: permanently blocked matches may be removed; surviving independent exponential clocks retain memoryless residuals. This supports an accepted-event race without materializing every rejected match. It preserves the accepted-event distribution, not the old seed-by-seed trajectory or rejected-event trace. Current linear transport channels can be combined for count observables; cell provenance needs additional treatment.

## Remaining development

Implement the accepted-event executor and test boundary handling, competition and snapshot reservation. Then prepare measured timing and stochastic-law qualification together where independent, with a concrete budget for approval. Particle-scale design depends on that qualification and performance evidence. The count arithmetic shares the existing independently assembled reference formulas, so subsequent validation must disclose this and retain the explicit identity engine as a comparison. No stochastic production is registered or started.

The preceding mixed-law qualification is PI CLOSED as MG127. Full numerical precision is retained in raw records; future plain-language reports use at most six significant digits unless more are necessary to determine an outcome.
