# ChatGPT review of MG002 — verbatim
Received 2026-10-27 via PI. Reviewed @ 7c4b194.

**I independently reproduced the 3D amplitudes and wavelength.** The requested quantities pass their registered tolerances. The overall **MG1 “FAIL as registered” verdict remains correct**, because I also reproduced the 1D stationarity failure.

I reviewed [the evaluation](https://github.com/RichardBenear/ModelG/blob/main/results/MG002/MG1_evaluation.md), [solver](https://github.com/RichardBenear/ModelG/blob/main/src/mg_continuum.py), and [gate result](https://github.com/RichardBenear/ModelG/blob/main/GATES.md). My implementation independently expands the perturbation equations, uses centered radial finite differences rather than Claude’s finite-volume operator, and uses Radau rather than BDF. Domain, seed, parameters, and final time match the specification.

At \(h=0.0125,\ t=100\):

| Quantity | Claude | Independent | Paper |
|---|---:|---:|---:|
| \(\phi_G(0)\) | −0.41147648 | −0.41148105 | −0.411 |
| \(\phi_X(0)\) | −14.63361600 | −14.63365708 | −14.6 |
| \(\phi_Y(0)\) | 1.69934082 | 1.69934939 | 1.70 |
| Pooled wavelength | 3.07930633 | 3.07930594 | 3.08 |
| \(\max|\partial_t\phi|\) | \(1.3090\times10^{-8}\) | \(1.3090\times10^{-8}\) | — |

The amplitudes agree between implementations within **0.0012%**. Refining my grid from \(0.025\) to \(0.0125\) changes the amplitudes, core radii, and wavelength by less than **0.01%**.

My additional 1D run gives \(\max|\partial_t\phi|=3.49318\times10^{-4}\), confirming failure against \(10^{-4}\).

Three review notes:

1. **Correct the seeded integral balance.**  
   In `analyze()`, the summed balance needs \(+\int_V\chi\,dV\) on its right-hand side. The documentation’s “every instant” claim is incorrect during seeding. This does **not** affect the reported \(t=100,300\) diagnostics, where the seed is negligible.

2. **Soften the inference about the paper.**  
   “The paper’s 1D column is itself a pre-stationary snapshot” is stronger than the evidence. Say: “Our reproduction of the published protocol remains nonstationary at \(t=100\), while matching the tabulated values.” Likewise, tail quadrature is a plausible explanation for the full-integral discrepancy, not an established cause.

3. **Report wavelength-window sensitivity.**  
   I reproduced Claude’s estimator exactly. Restricting instead to \(10<r<25\) gives **3.08540**; both estimates pass. The difference shows that the quoted four decimal places describe a particular window, rather than a uniquely determined wavelength.

**Disposition:** the 3D reproduction is independently confirmed. I support a separately documented PI exception or revised gate for longer 1D settling time, while preserving the original failure record. I have not independently replicated every integral or null-control result.

[Independent replication script](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review_work/independent_mg1.py) · [3D results](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review_work/independent_3d_h0.0125.json)

No repository files were changed.
