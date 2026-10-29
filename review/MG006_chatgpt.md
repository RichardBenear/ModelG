# ChatGPT review of the MG3 result — summary (verbatim, relayed by the PI 2026-10-29)

Note (Claude): the full review, answering all eight questions, is on the PI's machine (`review/MG006_replication/README.md`
in ChatGPT's workspace) and is not yet committed here. This file holds the summary the PI relayed.

---

I recommend preserving MG3's FAIL and proceeding with a separately registered follow-up. The main numerical results reproduce, and independent checks support spontaneous formation of persistent structures.
The full review, answering all eight questions, is in [MG006 replication review](C:/Users/rbene/.codex/.chatgpt-projects/g-p-6ae095a6c94c819198efcd6d84f02c03/review/MG006_replication/README.md).
I independently checked 168 completed run records and found no exposure, censoring, confirmation-time or early-stop bookkeeping mismatch. C1–C4 and V5 calculations reproduce. All 36 saved Stage-B first-event fields pass my separately implemented deterministic solver. A replay also reproduced two excursion records exactly; both independently decay to homogeneous.
The important qualifications are:

* C2(a) correctly fails: confirmation fractions are 5/7, 15/21 and 16/26. These describe outcomes through the first confirmed event, as specified before running—not detector quality over full trajectories.
* C3 supports an activated trend over Ω = 100–150. Three cells and five events in the weakest cell do not establish the mechanism. The reported "known variances" are estimated approximations.
* Stationarity remains weakly checked. Ω = 100 and b = 29.3 have no late-window exposure; their p = 1 values are uninformative.
* C4's registered PASS reproduces, but its highest rate rests on three surviving runs with waits of 5, 4 and 2 units after burn-in.
* V5 passed its registered comparison; timestep independence remains unresolved. Its broad interval and the reproduced startup asymmetry warrant follow-up.
* The interval formulas reproduce, but their claimed statistical exactness needs justification for stopped first-passage data. I have not demonstrated undercoverage here.

The Stage-A selection followed the rules. Its 4/4 trigger was poorly powered: at most a 2.2% probability per cell under the stated exponential model. The budget overrun should remain explicit.
For the follow-up, prioritize burn-in/timestep checks, more events per fitted cell, and uncertainty calculations calibrated for the actual sampling design. Keep any revision of the 80% threshold separate from this gate's result.
