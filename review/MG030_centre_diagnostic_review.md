# MG030 supplement: exploratory centre diagnostic review

**I accept the main interpretation: late motion is predominantly coherent translation of the sampled X-field profile. Amplitude redistribution or threshold-crossing motion alone does not explain it.** This supports the earlier three-field instantaneous diagnostic, but does not identify the cause or change MG030's registered outcomes.

Reviewed `results/MG9c/diagnostics/` and `src/diag_mg9c_centre.py` at revision `2417833460d0209c6a6f256d3487ef877226524f`, 2026-11-04. No integrations were run.

## Independent reproduction

I independently calculated both weighted centers, totals, variances and periodic-cut weight fractions at all 24 saved field times, using FFT convolution and no production imports. The largest moment discrepancy is 3.7e-12. I independently fitted the late translation with a coarse search followed by local refinement, and recalculated minima shifts. The production synthetic validation was also rerun and reproduces its results.

Changes over [15000,30000]:

| b | Boundary midpoint | E-squared center | X-potential-squared center | Profile alignment | Unexplained change norm |
|---|---:|---:|---:|---:|---:|
| 29.25 | -0.008105 | -0.008539 | -0.008134 | -0.008345 | 1.46% |
| 29.30 | -0.008826 | -0.009332 | -0.008706 | -0.009050 | 1.32% |
| 29.35 | -0.008693 | -0.009169 | -0.008536 | -0.008893 | 1.69% |
| 29.22, descriptive | -0.006545 | -0.006772 | -0.006587 | -0.006671 | 2.00% |

All center/mean-minimum changes are within approximately 4.1% of the alignment estimate. Crucially, interior minima move too. This is stronger evidence than agreement between centroids alone. The synthetic asymmetric-amplitude example illustrates why a centroid shift alone is insufficient.

The unexplained fraction is the Euclidean norm after optimal translation divided by the unaligned difference norm. It is not an energy fraction or conventional explained variance. Some deformation remains: use “predominantly translation,” not a claim that no redistribution occurs.

## Small corrections

1. **Use the same matched minima for mean and spread.** The late mean uses minima matched to the t = 1800 interior subset, but the late standard deviation/range discard the outermost minima at the late times. These subsets differ in the 18-cell cases. At b = 29.25 the reported anchored mean is -0.00836409, whereas the subset used for the spread has mean -0.00834537. Both support the conclusion. Select one matched subset and compute all its statistics together. Store its identities explicitly; `interior_min_x` currently contains all detected minima.
2. **Standard deviation is not a maximum bound.** “Each” minimum shifts by the same amount “to within about 1.5e-4” is too strong. Standard deviations are approximately 0.9–1.35e-4, but maximum deviations from the corresponding subset means reach 3.45e-4. Say “standard deviation about 1e-4” and retain ranges.
3. **Scope synthetic precision.** Better-than-1e-6 translation recovery applies to the chosen synthetic profiles. It is not a proven error bound on actual peak locations. An independent Fourier-interpolated minimum-location check changes individual late peak shifts by up to 1.51e-4 relative to the three-point parabolic estimates. Mean shifts remain close (within about 8e-6 for the late-order subsets), preserving the conclusion. Do not interpret all tiny peak-to-peak dispersion as physical deformation.

These corrections do not undermine the observed shift of roughly 0.008–0.009.

## Interpretation

The user's boundary-artifact explanation was a useful alternative to test. The late measurements argue against it as the dominant explanation: the weighted distributions, interior minima and aligned profiles move together. They do not distinguish preparation-induced phase relaxation, lattice alignment, persistent travel or numerical effects.

The plot's minima curves retain an early offset from t = 1800. A second panel rebased at t = 15000 would show the late agreement more clearly; the cumulative early discrepancy does not contradict the late table.

Six sparse snapshot times cannot show every intervening deformation or establish asymptotic speed. Keep this exploratory, preserve MG030's outcome, and do not assign physical energy/information meaning to the weights. A causal mesh/phase test requires separately specified shift, reflection, solver and mesh checks.

Replication: `review/MG030_replication/centre_audit.py` and `centre_audit.json`. Run `python review/MG030_replication/centre_audit.py CHECKOUT` with NumPy and SciPy. The checkout must include the committed MG9c snapshots, tracks and diagnostic JSON. The script imports no production modules and performs no ODE solves; Fourier minimum localization is an interpolation-sensitivity check, not ground truth.
