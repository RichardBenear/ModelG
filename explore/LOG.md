# Exploratory log

Exploratory-track entries under RESEARCH_PROTOCOL.md version 2. Write each entry **before** running: date, question,
method, compute cap (2 core-hours default; more only with PI approval), who runs it, and the PI's go-ahead. Append method
changes with reasons as they happen, then the result and a link to the note and records. Exploratory results are
labelled EXPLORATORY and cannot support a programme claim on their own.

| Date | ID | Question | Runner | Cap | PI go-ahead | Status / result |
|---|---|---|---|---|---|---|
| 2026-11-16 | EX001 | Does MG192's gravity drift check (linear-response mobility vs measured ETDRK4 drift) hold at MG193's operating point a = 14, b = 29.20526, including off-centre placements? Method: MG192 code with b = 29.20526; adjoint mobility mu_a; dynamic runs eps 0.001/0.002 at x0 = 0, +L/8, -L/8 plus nulls (eps 0 at x0 = 0, L/8); speed fit t in [100, 400]; compare with mu_a da/dx (MG192's 10% rule used as a reference, not a gate). Does not rescore MG192. | Claude | 2 core-h | Richard, 2026-11-16 09:56 America/Denver: "Go ahead and do the recheck of MG192 gravity dynamic with the new (a,b) values." | Running |
