# Projection diagnostics

Saved fits are replayed without minimization. Full-fit normalization is preserved in every window.
Signal windows contain the union of the central 99% intervals of PR and NP separately. Core windows contain the central 99% of PR; they need not contain 99% of NP.
Prefit lifetime windows use the central 99% of their own PDF. Mass uses the full fit window.

Display bin count is max(400, 2*sqrt(N_eff)), rounded up to a multiple of 50; full displays also retain at least the native bin count. The rule depends on statistics, never on chi-square. N_eff=(sum weights)^2/sum weights^2.
Dense signal views use twice the display bins. Native/full 60/120/240 and signal 60/120/240/480 scans remain in `summary.csv`, so fine binning cannot hide a broad residual. The table reports the dense display and the fixed 120-bin comparison together.
Expected bin contents are exact PDF integrals. max_center_bias_pull measures the largest bin-center approximation error in units of sqrt(expected); that approximation is not used for residuals.
Roo chi-square excludes empty observed bins and uses asymmetric data errors. Pearson/deviance include empty bins. Grouping merges adjacent expected counts to at least 5, with the final remainder merged backwards.
NDF = bins minus projection-relevant floating parameters is a descriptive approximation for constrained fits and restricted projections. No calibrated p-value or automatic fit acceptance is inferred.

| Stage | y | pT | Window | Display bins | Roo chi2/ndf | Roo / 120 bins | Grouped Pearson/ndf | Center bias (sigma) | Sparse bins | PDF |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|
| fit2d_ctau | 0.00–1.60 | 7.00–8.00 | -0.137–2.182 | 400 | 0.667 | 0.758 | 0.993 | 0.0345 | 243/400 | [plot](fit2d_ctau_y0.00_1.60_pt7.00_8.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 8.00–9.00 | -0.125–2.096 | 400 | 0.643 | 0.775 | 0.961 | 0.0367 | 259/400 | [plot](fit2d_ctau_y0.00_1.60_pt8.00_9.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 9.00–10.00 | -0.112–1.861 | 400 | 0.679 | 0.852 | 1.111 | 0.0298 | 269/400 | [plot](fit2d_ctau_y0.00_1.60_pt9.00_10.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 10.00–12.00 | -0.091–2.160 | 400 | 0.631 | 0.818 | 0.951 | 0.0484 | 259/400 | [plot](fit2d_ctau_y0.00_1.60_pt10.00_12.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 12.00–14.00 | -0.081–2.073 | 400 | 0.570 | 0.670 | 0.995 | 0.0438 | 292/400 | [plot](fit2d_ctau_y0.00_1.60_pt12.00_14.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 14.00–16.00 | -0.072–1.913 | 400 | 0.612 | 0.816 | 1.002 | 0.0253 | 341/400 | [plot](fit2d_ctau_y0.00_1.60_pt14.00_16.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 16.00–20.00 | -0.065–2.054 | 400 | 0.615 | 0.690 | 1.117 | 0.0434 | 338/400 | [plot](fit2d_ctau_y0.00_1.60_pt16.00_20.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 20.00–30.00 | -0.055–2.078 | 400 | 0.608 | 0.520 | 1.025 | 0.039 | 367/400 | [plot](fit2d_ctau_y0.00_1.60_pt20.00_30.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 1.00–2.00 | -0.696–3.119 | 1300 | 0.906 | 1.100 | 1.033 | 0.00124 | 0/1300 | [plot](fit2d_ctau_y1.60_2.40_pt1.00_2.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 2.00–3.00 | -0.592–2.278 | 1150 | 0.912 | 2.353 | 1.118 | 0.00198 | 169/1150 | [plot](fit2d_ctau_y1.60_2.40_pt2.00_3.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 3.00–4.00 | -0.289–2.333 | 950 | 0.852 | 1.517 | 1.174 | 0.00487 | 305/950 | [plot](fit2d_ctau_y1.60_2.40_pt3.00_4.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 4.00–5.00 | -0.212–2.224 | 700 | 0.712 | 1.033 | 0.992 | 0.0119 | 291/700 | [plot](fit2d_ctau_y1.60_2.40_pt4.00_5.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 5.00–6.00 | -0.182–1.937 | 450 | 0.974 | 1.505 | 1.428 | 0.0285 | 200/450 | [plot](fit2d_ctau_y1.60_2.40_pt5.00_6.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 6.00–7.00 | -0.169–1.797 | 400 | 0.801 | 1.097 | 1.173 | 0.0304 | 209/400 | [plot](fit2d_ctau_y1.60_2.40_pt6.00_7.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 7.00–8.00 | -0.144–1.850 | 400 | 0.762 | 1.260 | 1.219 | 0.0327 | 247/400 | [plot](fit2d_ctau_y1.60_2.40_pt7.00_8.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 8.00–9.00 | -0.148–2.023 | 400 | 0.722 | 0.834 | 1.157 | 0.0246 | 282/400 | [plot](fit2d_ctau_y1.60_2.40_pt8.00_9.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 9.00–10.00 | -0.125–2.054 | 400 | 0.637 | 0.601 | 1.121 | 0.025 | 313/400 | [plot](fit2d_ctau_y1.60_2.40_pt9.00_10.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 10.00–12.00 | -0.117–1.750 | 400 | 0.789 | 0.891 | 1.316 | 0.0333 | 291/400 | [plot](fit2d_ctau_y1.60_2.40_pt10.00_12.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 12.00–14.00 | -0.100–1.854 | 400 | 0.684 | 0.770 | 1.061 | 0.0321 | 344/400 | [plot](fit2d_ctau_y1.60_2.40_pt12.00_14.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 14.00–20.00 | -0.093–1.840 | 400 | 0.711 | 1.230 | 1.400 | 0.0411 | 349/400 | [plot](fit2d_ctau_y1.60_2.40_pt14.00_20.00.pdf) |
| fit2d_mass | 0.00–1.60 | 7.00–8.00 | 2.600–3.500 | 400 | 0.725 | 0.635 | 0.884 | 0.00396 | 0/400 | [plot](fit2d_mass_y0.00_1.60_pt7.00_8.00.pdf) |
| fit2d_mass | 0.00–1.60 | 8.00–9.00 | 2.600–3.500 | 400 | 0.703 | 0.821 | 0.943 | 0.00386 | 0/400 | [plot](fit2d_mass_y0.00_1.60_pt8.00_9.00.pdf) |
| fit2d_mass | 0.00–1.60 | 9.00–10.00 | 2.600–3.500 | 400 | 0.671 | 0.891 | 1.049 | 0.00403 | 247/400 | [plot](fit2d_mass_y0.00_1.60_pt9.00_10.00.pdf) |
| fit2d_mass | 0.00–1.60 | 10.00–12.00 | 2.600–3.500 | 400 | 0.612 | 0.715 | 0.927 | 0.00498 | 134/400 | [plot](fit2d_mass_y0.00_1.60_pt10.00_12.00.pdf) |
| fit2d_mass | 0.00–1.60 | 12.00–14.00 | 2.600–3.500 | 400 | 0.545 | 0.686 | 1.046 | 0.00379 | 324/400 | [plot](fit2d_mass_y0.00_1.60_pt12.00_14.00.pdf) |
| fit2d_mass | 0.00–1.60 | 14.00–16.00 | 2.600–3.500 | 400 | 0.538 | 0.736 | 1.202 | 0.00248 | 347/400 | [plot](fit2d_mass_y0.00_1.60_pt14.00_16.00.pdf) |
| fit2d_mass | 0.00–1.60 | 16.00–20.00 | 2.600–3.500 | 400 | 0.501 | 0.527 | 0.833 | 0.00194 | 347/400 | [plot](fit2d_mass_y0.00_1.60_pt16.00_20.00.pdf) |
| fit2d_mass | 0.00–1.60 | 20.00–30.00 | 2.600–3.500 | 400 | 0.446 | 0.504 | 0.718 | 0.00155 | 359/400 | [plot](fit2d_mass_y0.00_1.60_pt20.00_30.00.pdf) |
| fit2d_mass | 1.60–2.40 | 1.00–2.00 | 2.600–3.500 | 1350 | 1.102 | 1.089 | 1.098 | 4.46e-05 | 0/1350 | [plot](fit2d_mass_y1.60_2.40_pt1.00_2.00.pdf) |
| fit2d_mass | 1.60–2.40 | 2.00–3.00 | 2.600–3.500 | 1150 | 1.019 | 0.955 | 1.014 | 0.000103 | 0/1150 | [plot](fit2d_mass_y1.60_2.40_pt2.00_3.00.pdf) |
| fit2d_mass | 1.60–2.40 | 3.00–4.00 | 2.600–3.500 | 950 | 0.977 | 0.994 | 0.986 | 0.000182 | 0/950 | [plot](fit2d_mass_y1.60_2.40_pt3.00_4.00.pdf) |
| fit2d_mass | 1.60–2.40 | 4.00–5.00 | 2.600–3.500 | 700 | 0.959 | 1.028 | 0.985 | 0.000383 | 0/700 | [plot](fit2d_mass_y1.60_2.40_pt4.00_5.00.pdf) |
| fit2d_mass | 1.60–2.40 | 5.00–6.00 | 2.600–3.500 | 450 | 0.955 | 1.182 | 1.021 | 0.000888 | 0/450 | [plot](fit2d_mass_y1.60_2.40_pt5.00_6.00.pdf) |
| fit2d_mass | 1.60–2.40 | 6.00–7.00 | 2.600–3.500 | 400 | 0.849 | 0.932 | 0.963 | 0.000994 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt6.00_7.00.pdf) |
| fit2d_mass | 1.60–2.40 | 7.00–8.00 | 2.600–3.500 | 400 | 0.780 | 1.009 | 0.954 | 0.000813 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt7.00_8.00.pdf) |
| fit2d_mass | 1.60–2.40 | 8.00–9.00 | 2.600–3.500 | 400 | 0.749 | 0.877 | 0.957 | 0.00077 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt8.00_9.00.pdf) |
| fit2d_mass | 1.60–2.40 | 9.00–10.00 | 2.600–3.500 | 400 | 0.620 | 0.837 | 0.938 | 0.000627 | 122/400 | [plot](fit2d_mass_y1.60_2.40_pt9.00_10.00.pdf) |
| fit2d_mass | 1.60–2.40 | 10.00–12.00 | 2.600–3.500 | 400 | 0.629 | 0.831 | 1.006 | 0.000674 | 123/400 | [plot](fit2d_mass_y1.60_2.40_pt10.00_12.00.pdf) |
| fit2d_mass | 1.60–2.40 | 12.00–14.00 | 2.600–3.500 | 400 | 0.429 | 0.634 | 0.938 | 0.000428 | 322/400 | [plot](fit2d_mass_y1.60_2.40_pt12.00_14.00.pdf) |
| fit2d_mass | 1.60–2.40 | 14.00–20.00 | 2.600–3.500 | 400 | 0.486 | 0.659 | 0.934 | 0.000445 | 327/400 | [plot](fit2d_mass_y1.60_2.40_pt14.00_20.00.pdf) |
