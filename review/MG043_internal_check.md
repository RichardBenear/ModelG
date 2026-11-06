# MG043 result check — UNRESOLVED, ready for PI assessment

2026-11-06. Lead implementation author: ChatGPT / Codex. External review was waived.

The complete registered matrix ran after 88eabd2 was committed and pushed. All
64 paths completed within the local eight-worker limits. No outcome was replaced
and no capacity, solver, threshold or stopping condition was changed after results.
All production accounting and availability checks pass. The separate checker
agrees with all 14 arm labels, the batch label, raw hashes and budget.

Both largest capacities recover for both starts. The registered overall outcome
is nevertheless UNRESOLVED because the C=1 homogeneous path has E2/O discrepancy
0.4996807 against <=0.001. It is neither a production implementation/ledger failure
nor a resolved physical recovery failure. The exact-homogeneity and detector
problems are separately discussed in results/MG10/MG10_results.md.

Internal challenge: the independent ODE solver and Euler preserve conservation
while responding differently to extremely small symmetry-breaking perturbations;
conservation does not establish trajectory agreement. Also, the frozen detector's
use of the undepleted baseline admits FFT-roundoff minima once an exactly uniform
background is below threshold. Shared detector dependencies mean independent
arithmetic recomputation alone cannot remove that limitation. Reported counts in
that arm are not evidence of physical particles. This qualification-coverage gap
requires new measurement validation in a follow-up, not a retrospective rule change.

Posthoc diagnostics and the summary figure use only retained raw data. No additional
trajectory, revised detector or alternative scoring verdict has been introduced.
The original MG041 failure, the corrected MG042 qualification, the registered
outcome and every completed run remain separately retained. The report states
same-author recomputation and physical/modeling limits. PI closure remains pending.
