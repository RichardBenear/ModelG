# MG193 — Continuum 1D existence and stability window of the Model G soliton in (a, b)

Status: **REGISTERED 2026-11-16 (UTC); frozen before computation.** Author: Claude, running it at the PI's request.
PI approval: Richard Benear, 2026-11-16 06:14 America/Denver: "Go ahead with the existence-window mapping." This is a
descriptive mapping study. It has numerical qualification checks but no hypothesis gates.

## 1. Purpose

MG192 found that the 1D continuum soliton at the published point (a = 14, b = 29) sits about 0.011% below a fold in a.
That was already partly known: the Phase 2 protocol (MG003, §2 of `spec/phase2_protocol_v1.3.4.md`) gives an operational
continuum lower edge b_s = 28.993 in 1D, 0.024% below b = 29. It also notes that b = 29 is close to the edge only in 1D
(2D 27.343, 3D 26.972). The upper edge in the continuum, and the shape of the window in (a, b), are not known.

MG193 maps where a single stationary 1D soliton exists and is linearly stable, as both a and b vary. The output is a
well-inside operating point for later 1D studies (gravity dynamics, two-particle force), with the margins to every edge
stated.

## 2. Method

- **Numerics:** the same equations and Fourier spectral discretization as MG192 (`src/mg192_mobility.py`, reused
  unchanged): periodic L = 100, N = 1024; u = s = w = 0, p = q = 1, g = 0.1, d_x = 1, d_y = 12; A and B constant.
- **Continuation:** natural-parameter continuation from the MG192 stationary soliton, with Newton in the even
  subspace.
  - Converged means max residual ≤ 1e-10, positivity holds, and the core depth φ_X(0) ≤ −1, so the solution has not
    slid back to the homogeneous state.
  - The step halves on failure. An edge is bracketed when the step falls below 1e-5 (absolute in b; relative in a).
- **Stability at each converged point:** full eigenvalues of the even- and odd-subspace Jacobians (dense), after
  removing the translation zero mode, which is the odd eigenvalue nearest 0. Linearly stable means every remaining
  eigenvalue has Re λ < 0. Where stability changes sign between steps, the crossing is located by bisection to 1e-5.
- **Lines mapped:**
  1. b varied at a = 12, 13, 14, 15, 16, from b = 29 (or the nearest stable point found) both up and down. The upward
     search stops at the homogeneous Turing threshold b_c(a).
  2. a varied at b = 29, both up and down.
  - If the soliton does not exist at b = 29 for some a, that line starts from the a = 14 solution, continued in a at
    the b value that is the a = 14 window centre.
- **Edge labels** at each end of the window:
  - **FOLD:** Newton fails as the even eigenvalue approaches 0.
  - **STABILITY LOSS:** an eigenvalue crosses 0 while the solution still exists; even or odd is reported.
  - **TURING:** the line reaches b_c.
  - **DECAY:** the core depth rises above −1.

## 3. Outputs

- For each a: b_low, b_high, the type of each edge, and the stable window width.
- The window in a at b = 29.
- A recommended 1D operating point: on the a = 14 line, the window centre in b. Where it lies in a, and its margins to
  all four edges, are also reported.
- Comparison with the Phase 2 operational estimates (continuum b_s = 28.993; lattice h = 0.2 [28.929, 29.397];
  h = 0.1 [28.978, 29.412]). The definitions differ, so this is context only.

## 4. Numerical qualification (reported; failure means the affected edge is labelled UNQUALIFIED)

- Both a = 14 edges in b are recomputed at N = 2048 and at L = 200 (N = 2048). Each must agree with N = 1024 / L = 100
  within 1e-3 in b.
- At the recommended operating point, the soliton is checked by a time-dependent ETDRK4 run from the stationary profile
  plus a 1% random even-and-odd perturbation (seed 20261116), to T = 200. The core depth must stay within 1% of the
  stationary value at T = 200.

## 5. Limits

- 1D continuum only. 2D and 3D windows are known only from the Phase 2 operational lower edges.
- Single-soliton existence and linear stability only. No basin sizes, no noise, no multi-particle states.
- The mapped edges depend on the periodic domain L = 100 only weakly; the L = 200 check measures this.

## 6. Execution and records

- **Where:** Claude's cloud container, at most 4 processes, cap 4 core-hours.
- **Code:** `src/mg193_window.py`.
- **Records:** `results/MG193_existence_window/` (continuation tables as JSON, a figure and the report).
- **Separate check:** each a = 14 edge in b, and both b = 29 edges in a, are reconfirmed by time-dependent ETDRK4
  integration from the nearest converged profile, run once just inside and once just outside the edge. The offset is
  0.01 in b, or 0.05% relative in a. Inside, the core depth must stay within 1% of its start at T = 500. Outside, the
  particle must dissolve (depth above −1) or change shape (another deep minimum, or depth changed by more than 10%) by
  T = 500. An edge whose check is inconclusive is reported as such.
- No result is rescored after the data are seen.
