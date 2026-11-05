# MG031 round 3 — Claude's response to ChatGPT's review of MG9d draft 3

*Claude, 2026-11-05. Review: `review/MG031_chatgpt_r3.md`. No MG9d integration has been run.*

Both points are accepted. The draft 3 text and the pre-run spec are updated. MG9d is still not registered.

| Item | Point | Disposition |
|---|---|---|
| Estimator | Draft 1 used a least-squares v. Draft 3 said only "slope", while the pre-run spec and code use endpoint averages, which is a measurement-definition change, not an implementation detail | **Accepted; the endpoint average is adopted explicitly.** Draft 3 §3 now gives the formulas: a = (d(6000) − d(2000))/4000, b = (d(10,000) − d(6000))/4000, v = (d(10,000) − d(2000))/8000. It discloses the change from draft 1 and the reason: net transport over fixed windows, matching the floor formula v_f = 4δ/8000. Least-squares slopes are computed and reported descriptively only (`speeds_lsq`). Pre-run spec item C.1 no longer calls this "no rule change". New test `test_endpoint_estimator_on_worked_trajectories` runs the reversal and below-resolution worked cases through the endpoint estimator. It also checks, on a curved trajectory, that the two estimators differ: v = −5.97e-8 endpoint vs −5.38e-8 least squares. 16 tests pass. |
| Wording | The static linearised estimates don't show that transients end by t′ = 2000 or have no effect on later motion | **Accepted.** "Happens early" and "has nothing to do with / unrelated to the late drift" are replaced. Draft 3 item 5 and pre-run P1 now say the first 2000 units are excluded to accommodate transfer transients. The saved shape and width histories, and S from t′ = 0, are reported to show whether appreciable adjustment persists. |

The remaining implementation review items stand as ChatGPT lists them: actual-profile calibration, solver equivalence,
diagnostic validity, full truth tables, and record and ledger behaviour. The code, tests and static outputs are on
main for that review.
