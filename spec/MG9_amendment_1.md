# MG9 amendment 1 — envelope by complex demodulation (replaces the draft 3 §6.2 kernel)

*Stage MG027. Claude, 2026-11-04. Status: **PROPOSED, awaiting PI sign-off.** It amends `spec/MG9_proposal.md` draft 3
(registered, frozen at `7810b8b`). No Model G MG9 run has been made, before or after this amendment.*

## Why

Draft 3 §6.2 requires a synthetic validation of the envelope before any Model G run. It also says: "Each failure blocks
the study until the kernel is fixed and re-registered." While implementing it (`tests/test_mg9.py`), the registered
kernel **failed** that validation:

| Envelope | Unit cosine, 10 phases (limit [0.98, 1.02]) | Worst tracked-edge error on known patches (limit λ/10) |
|---|---|---|
| Draft 3: E = √(2 K_λ[(φ_X − K_λφ_X)²]) | 0.9996–1.0011 (D02); 0.9999–1.0000 (D01): pass | 0.27 λ on both meshes: **fail** |
| Amendment 1 (below) | 1.0000–1.0000 on both: pass | 0.036 λ on both meshes: pass |

The known patches are cos(kx + phase) · ½[tanh((x + a)/ℓ) − tanh((x − a)/ℓ)], with:
- 10 phases;
- a = 25 and 17.3;
- ℓ = λ/2 (the production splice), λ/8 and 2λ.

The true edge is the half-amplitude point x = ±a. The failure has two parts:
- **A phase-dependent error** of up to 0.15 λ. Near an edge, a single λ-wide average of φ_X² contains a partial
  oscillation of varying amplitude.
- **A bias** of up to λ/4 for sharp edges. A constant bias would cancel in displacements, but the phase-dependent part
  would not, and it is comparable to the λ/4 arrest margin.

These numbers are synthetic only (no Model G field). ChatGPT's round 2 named the alternative: "a specified
demodulation and low-pass procedure".

## Change

§6.2's envelope definition is replaced by:

- K_λ is the physical-width boxcar of draft 3 (width exactly λ, fractional end weights, periodic), and K² = K_λ ∘ K_λ
  (a triangular low-pass of width 2λ).
- m = K²[φ_X]; z(x) = (φ_X − m) · e^{−ikx} with k = 2π/λ and x = (j − N/2) h; **E = 2 |K²[z]|** (K² applied to the real
  and imaginary parts).

For a(x) cos(kx + phase) with slowly varying a, E is the low-passed a, independent of phase. K² removes the e^{−2ikx}
term exactly for constant a (K_λ has zeros at all multiples of 2π/λ). e^{−ikx} is periodic on every domain, because
kL = 2π·34 (D02, D01) or 2π·68 (D02L).

Everything that uses E is unchanged: E_pat, the ½ E_pat patterned set, components, interfaces, the cell-count
prominence, and the classification. The synthetic validation and its limits are unchanged. Its patch set is now fixed
as listed above.

## Effect

The observable changes before any data exist. No gate criterion, threshold or outcome rule changes. The draft 3 kernel is
kept in the code as `envelope_rms`, with a test documenting its failure (`test_registered_rms_envelope_fails_its_validation`).
