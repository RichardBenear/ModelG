# ChatGPT review of Amendment 4 (verbatim, relayed by the PI 2026-10-28)

I recommend accepting the stated residual risk and proceeding to MG3, rather than requiring a stronger pooled audit. Section 2 now honestly describes what the replacement can and cannot establish. This is a PI risk decision, not an unresolved implementation error. [Amendment 4](https://github.com/RichardBenear/ModelG/blob/main/spec/phase2_amendment_4.md)
I checked the replacement implementation and calibration script. The leave-one-out statistic includes the correct prediction-variance factor, and the Bonferroni correction counts every supplied run × field × direction test. The reported calibration results match the saved output. [Implementation](https://github.com/RichardBenear/ModelG/blob/main/src/stats_mg2.py), [calibration](https://github.com/RichardBenear/ModelG/blob/main/src/calibrate_audit.py)
The important distinction is:

* The old audit was defective: its statistic could never reach the rejection threshold.
* The replacement works but adds no detection power: anything it flags would already fail the less stringent per-cell run-level check.
* Rare modes remain an assumption: passing these diagnostics does not establish their absence or validate confidence bounds against them.

A stronger pooled test could improve sensitivity by borrowing scale information, but that information must be justified. Within-run block variance is not automatically the variance of a run mean: temporal correlation and variation between runs matter. Adding that approach now would introduce another validation problem, with no demonstrated improvement for this experiment.
I suggest two small wording qualifications before closure:

1. The reported power curve applies to one displaced run in one scalar projection under the simulated Gaussian model. It does not measure sensitivity to every kind of rare mode.
2. "A more powerful pooled audit would need an independent scale estimate" is too categorical. Prefer: "The proposed route to a stronger pooled audit would require a separately justified scale estimate." Other statistical designs are possible.

Neither qualification requires another calibration campaign. My recommendation is to approve Amendment 4 with those clarifications, preserve MG2's registered FAIL, and carry the explicit A1 limitation into MG3's eventual interpretation.
