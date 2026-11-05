# MG9d — pre-run implementation specification (stage MG032)

*Claude, 2026-11-05. Status: written **before any MG9d integration**. MG9d is not yet registered. The PI approved
draft 3 for implementation (decision card, 01:42 UTC). This document covers ChatGPT's pre-run implementation review
and two proposed changes to validity items that need PI decisions before registration. Everything here comes from
static computation only: maps, linear algebra and recovery of known shifts. No time-stepping of an MG9d state was
done.*

Code: `src/run_mg9d.py` (driver, diagnostics, rules), `src/evaluate_mg9d.py` (re-evaluation from the records), and
`tests/test_mg9d.py` (16 tests, all pass, about 5 s). Static outputs: `results/MG9d/transfer.json` and
`results/MG9d/floor.json`. `python src/run_mg9d.py runs` refuses to start unless `results/MG9d/REGISTERED` exists. That
file is written only after registration.

## A. Static results (before freezing)

**State-transfer checks (§2): all pass** (`transfer.json`).

| Check | Value | Limit |
|---|---:|---:|
| Source SHA-256 against the MG9c manifest | match (`f22ca280…`) | — |
| y₀ − source (Nyquist filter) | 1.56e-5 | 2e-5 |
| Shift round-trips h/4, h/2, 3h/4 | 8.0e-15, 8.9e-15, 8.0e-15 | 1e-12 |
| Reflection applied twice | 0 (exact) | 1e-12 |
| R4 sampled back at the D02 nodes | 5.3e-15 | 1e-12 |

**Measurement floor (§3): δ = 1.11e-9, so v_f = max(5e-9, 4δ/8000) = 5e-9** and the resolved threshold is 1.5e-8. The
known drift (about 4e-7) is about 27 times that threshold. Each measure's largest static recovery error over the 24
cases (both grids, ±0.001 to ±3h/4) is:

| Measure | d (alignment) | C_edge | C_amp | C_phi2 | m (minima) |
|---|---:|---:|---:|---:|---:|
| Largest error | 1.1e-9 | **2.3e-4** | 1.0e-9 | 6.9e-9 | 7.6e-5 |

All 16 interior minima matched in every case, and no case was ALIGNMENT-AMBIGUOUS.

**Initial right-hand side, max |RHS| (recorded as §2 asks).**

| State | Value |
|---|---:|
| Source (unfiltered) | 1.09e-5 |
| y₀ (R0, R1, R3) | 1.65e-3 |
| R2a, R2b, R2c | 0.26, 0.37, 0.26 |
| R4 (D01) | 4.85 |

A linearised static estimate (pseudo-inverse correction −J⁺·RHS, dropping the near-zero translation mode) gives the
size of the shape adjustment each start state needs, as ‖Δφ_X‖ / ‖φ_X‖:

| State | Relative adjustment |
|---|---:|
| R0 | 3.8e-6 |
| R2a | 8.5e-4 |
| R2b | 1.0e-3 |
| R4 | 1.8e-2 |

## B. Findings and changes (PI decisions)

**PI decision (2026-11-05 01:50 UTC, decision card "Adopt P1 and P2"): P1 and P2 adopted; P0 written into draft 3 §0.**
Draft 3 §3 items 5–6 now read as below, and the code already implements them.

**P0. A disclosure correction, no rule change.** Draft 3 §0 says the filter changes the initial residual "slightly". In
fact it raises max |RHS| from 1.09e-5 to 1.65e-3, which is 150 times higher. The cause is that the removed Nyquist
component is the most strongly damped lattice mode. The shape adjustment it implies is small (3.8e-6 relative), and the
2000-unit exclusion before the speed windows covers it. The registered text should give the numbers in place of
"slightly".

**P1. Measure validity items 5 and 6 from t′ = 2000, not from t′ = 0 (proposed).** The sub-grid-shifted and
finer-mesh starts are not lattice equilibria. The shape adjustment estimated above (8.5e-4 to 1.0e-3 for R2, 1.8e-2 for
R4) is expected to occur mainly as an early transfer transient. Item 5 as written, S(T) ≤ 1e-3 measured from t′ = 0, would
probably mark R4 INVALID and could mark R2b INVALID by construction. Item 6, consistency at T from t′ = 0, has the same
problem: the early adjustment moves C_edge and the centroids differently from d. The speed windows already start at
2000.

Proposal:
- Item 5 becomes S_ref(T) = ‖φ_X(T) − T_{d(T)−d(2000)} φ_X(2000)‖ / ‖φ_X(2000)‖ ≤ 1e-3.
- Item 6 compares the changes in each measure over [2000, T] with d(T) − d(2000).
- Identity (item 3) stays measured from t′ = 0.
- S from t′ = 0 is still recorded and reported descriptively.
- The alternative is to keep draft 3 as written and accept that R4 (and perhaps R2b) may come out INVALID because of
  the start-up adjustment rather than the late motion. An invalid R4 limits only C5, but an invalid R2 arm makes C2 INCONCLUSIVE.

The code implements the proposal (`VALID_REF_T = 2000`). Draft 3's version is `validity(rec, t_ref=0)`.

**P2. Raise item 6's absolute floor from 2e-4 to 5e-4 (proposed).** C_edge's static recovery error reaches 2.3e-4, which
is above draft 3's absolute floor of 2e-4. The half-level crossing is interpolated linearly on the grid, so its error
depends on sub-grid position. With the expected late displacement of about 0.0016 over [2000, 6000] and 0.0032 over
[2000, 10000], the band 10%|Δd| + 2e-4 is 5.2e-4. That is only about twice C_edge's static error, so a valid run could
fail item 6 from edge-measurement error alone. 5e-4 is about twice the largest static error of any measure.

## C. Implementation choices (for review; no rule changes)

1. **Slopes are endpoint averages** of d: a = (d(6000) − d(2000))/4000, and likewise b and v. This is the estimator
   the floor formula v_f = 4δ/8000 assumes: two endpoint errors of at most δ over 4000 units. **This is a
   measurement-definition change from draft 1's least-squares v (ChatGPT MG031 r3), now written into draft 3 §3 and
   disclosed there.** Least-squares slopes are reported descriptively only. The static linearised estimates in §A do not
   show that every transient ends by t′ = 2000; the saved shape and width histories will show whether appreciable
   adjustment persists.
2. **Tracker** (item 2) is checked at every 5-unit sample, not only every 50.
3. **Minima m(t′):** the reference is the interior minima at t′ = 0 (all but the outermost minimum on each side). They
   are matched by nearest position within λ/4. m uses the one subset matched at every sample, and the identities are
   stored. For R3, all quantities are in R3's own (reflected) coordinates, and the sign is flipped once, in
   `evaluate` (`reflect_adjust`).
4. **Alignment:** bounded search ±λ/4 around the previous sample's d, with xatol 1e-10. A sample is
   ALIGNMENT-AMBIGUOUS within 1e-3 of a bound. T_s is the Fourier shift on the run's own grid.
5. **Solver:** `mg9d_ode` reproduces `run_mg9.mg9_ode` bit for bit at max_step 0.5 (tested on MG003's b = 29
   reference over 5 units). Segments are 1000 units, and t_eval is every 5 units in every run.
6. **Ledger:** an integration is counted when it starts, so an interrupted run still counts toward the cap of 10. A
   run that completes with complete = False is INVALID (item 1) and is not re-run. The three spare integrations are
   only for infrastructure interruptions, each disclosed.
7. **Snapshots:** the full three-field state every 1000 units including t′ = 0 (11 per run), saved as `.npy` during
   the run (gitignored), then packed into `results/MG9d/snapshots.npz` with a SHA-256 manifest. About 1 MB.
8. **E_pat** = 6.60086115019687, read at run time from the MG9c record's top-level field `E_pat`. Draft 3 §3 says
   `info.E_pat`; that is a wording slip, now corrected in the draft text.
9. **Platform:** cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1. It is recorded in every record.

## D. Tests (`tests/test_mg9d.py`, 16; the 16th, added after MG031 r3, runs the worked cases through the endpoint estimator)

- Wrapper bit-for-bit match.
- Exact E_pat and the source hash.
- The Nyquist filter.
- The shift, reflection (x → −x) and D01 maps.
- The transfer checks.
- Alignment recovery and the ambiguity flag.
- The floor rule.
- The AGREE/DIFFER bands, including a tiny sign change near zero.
- The truth tables for C1, C2, C3 (sign applied once, NOT EVALUATED), C4 (every branch, the ρ gap and the R1-must-match
  rule) and C5 (every branch, including draft 1's overlap case).
- The overall order, with ChatGPT's three worked cases. This covers mesh evidence not being hidden by a reversal,
  NO RESOLVED DRIFT requiring R1, and the qualification when C5 is unavailable.
- An end-to-end run with a fake ODE (exact translation at −4e-7). It checks the registration guard, the record
  structure (201 diagnostics, 2001 track samples, 11 snapshots), speeds recovered to 1e-11, validity, and the cap
  counted at start.
