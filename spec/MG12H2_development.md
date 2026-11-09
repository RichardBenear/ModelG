# MG12-H2 fixed incidence transport — MG110 development

2026-11-09, following PI closure MG109. Build/test fixed higher-arity transport and known-dimensional diagnostics. Closed H1 engine and qualification are unchanged. No particle or topology-evolution test is authorized by this development step.

## Local law

Use typed unordered ActiveCell/ACell/BCell relations with distinct incident sites, arity k>=2 and positive weight w. For a participating species S token at site i, one match per other incident site j has rate D_S*w/(k-1), or kappa_S*w/(k-1) for A/B. Thus each incident cell supplies total exit rate D_S*w, independent of arity. A node incident on several cells has the sum of their exit rates; this is not degree-normalized. Equal site volumes/capacities only. The inverse-square spacing factor belongs explicitly in w for regular-grid fixtures.

Matches retain original cell provenance even when two cells yield identical source/destination channels. Tick conflicts reserve the token, not shared read-only cells/sites. Chemical matches and catalyst reservations retain the frozen H1 law. Overlapping cells remain separate channels; no matching multiplicity is lost through dictionary keys. Pure pairwise H1 links may coexist with cells.

This linear single-token transport is **exactly reducible to a weighted pairwise generator**, with off-diagonal rates summed over common cells. That is an expected control, not a failure and not evidence of irreducible higher-order dynamics. Higher-arity incidence is explicit, but cooperative multi-token transport would be a separately registered new kinetic model. Do not add it solely to force a distinction.

## Bounded development

One worker/library-thread1,300elapsed-second cap including sampling/integration/IO/checks; no retry after failure. Freeze source/config/seed manifest before RNG or matrix integration. Six fixtures:4site binary ring;3site ternary cell;4sites overlapping ternary/quaternary ActiveCells;3x3periodic2D grid;3x3x3periodic3Dgrid;4site A/B tokens on distinct overlapping resource cells. Active fixtures carry one G with chemistry off; resource fixture carries one A and one B. Fixed times0,.05,.1; exact and strict ticks.005,256paths per fixture/method=3072development paths. Disjoint fixed seeds with deterministic collision resolution against prior plans. Retain counts/traces/seeds/config and checkpoints; no particle outcomes read.

Independent incidence-matrix construction derives pairwise rates by vector algebra rather than enumerating engine matches. Compare complete generator/event increments for each fixture's finite count states to relative/absolute1e-12. Check constant field invariance, row sums, symmetry for this equal-volume undirected domain, binary rollback to H1, overlapping-cell exit rates, relabeling/corresponding marks, conservation, strict one-use and source snapshot. Compare reference matrix exponentials to analytic regular-grid Fourier-mode decay within1e-11. No hidden renormalization of degree or time.

Endpoint exact/strict distances are development diagnostics, not formal qualification. Exact executor CDF smoke comparison to the independent finite law uses alpha.01 across6fixtures and256paths; strict distances descriptive. Conditional on registered PRNG/reference assumptions only. No favorable step/sample selection, no monotonic Monte Carlo requirement.

Dimension calibration is deterministic: Manhattan balls on known1D/2D/3D integer lattices, measured at prespecified radii4and8 without periodic wrapping. Independently count by direct coordinate enumeration and closed combinatorial formulas. Report finite-scale log growth slopes; they need not equal the asymptotic dimension at small radii. Arity, site count and embedding dimension are not interchangeable. No dimension is inferred for the tiny irregular cells from this calibration.

H2development completion supports a subsequent formal transport qualification design, not immediate particle runs. Any extension of strict-tick accuracy beyond H1's fixtures requires its own measured criteria/replication budget. Explicit high-count match expansion still needs later optimization and equivalence checks. All results begin with a plain-language project contribution/limits summary; PI closure and further runs remain separate.
