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
| fit2d_ctau | 0.00–1.60 | 7.00–8.00 | -0.129–2.172 | 400 | 0.681 | 0.755 | 1.124 | 0.0319 | 243/400 | [plot](fit2d_ctau_y0.00_1.60_pt7.00_8.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 8.00–9.00 | -0.126–2.136 | 400 | 0.733 | 0.682 | 1.176 | 0.0355 | 261/400 | [plot](fit2d_ctau_y0.00_1.60_pt8.00_9.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 9.00–10.00 | -0.114–1.930 | 400 | 0.720 | 0.976 | 1.145 | 0.0312 | 272/400 | [plot](fit2d_ctau_y0.00_1.60_pt9.00_10.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 10.00–12.00 | -0.094–2.155 | 400 | 0.705 | 0.802 | 1.120 | 0.0462 | 259/400 | [plot](fit2d_ctau_y0.00_1.60_pt10.00_12.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 12.00–14.00 | -0.078–2.064 | 400 | 0.598 | 0.779 | 1.050 | 0.0426 | 291/400 | [plot](fit2d_ctau_y0.00_1.60_pt12.00_14.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 14.00–16.00 | -0.070–1.902 | 400 | 0.608 | 0.773 | 1.051 | 0.024 | 340/400 | [plot](fit2d_ctau_y0.00_1.60_pt14.00_16.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 16.00–20.00 | -0.064–2.058 | 400 | 0.606 | 0.754 | 1.048 | 0.0422 | 338/400 | [plot](fit2d_ctau_y0.00_1.60_pt16.00_20.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 20.00–30.00 | -0.054–2.085 | 400 | 0.610 | 0.640 | 0.904 | 0.0369 | 368/400 | [plot](fit2d_ctau_y0.00_1.60_pt20.00_30.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 1.00–2.00 | -0.696–3.119 | 1300 | 0.906 | 1.100 | 1.033 | 0.00124 | 0/1300 | [plot](fit2d_ctau_y1.60_2.40_pt1.00_2.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 2.00–3.00 | -0.592–2.278 | 1150 | 0.912 | 2.353 | 1.118 | 0.00198 | 169/1150 | [plot](fit2d_ctau_y1.60_2.40_pt2.00_3.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 3.00–4.00 | -0.283–2.332 | 950 | 0.844 | 1.545 | 1.154 | 0.00494 | 306/950 | [plot](fit2d_ctau_y1.60_2.40_pt3.00_4.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 4.00–5.00 | -0.214–2.225 | 700 | 0.718 | 0.958 | 0.996 | 0.0114 | 291/700 | [plot](fit2d_ctau_y1.60_2.40_pt4.00_5.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 5.00–6.00 | -0.182–1.937 | 450 | 0.974 | 1.505 | 1.428 | 0.0285 | 200/450 | [plot](fit2d_ctau_y1.60_2.40_pt5.00_6.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 6.00–7.00 | -0.169–1.797 | 400 | 0.801 | 1.097 | 1.173 | 0.0304 | 209/400 | [plot](fit2d_ctau_y1.60_2.40_pt6.00_7.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 7.00–8.00 | -0.144–1.850 | 400 | 0.762 | 1.260 | 1.219 | 0.0327 | 247/400 | [plot](fit2d_ctau_y1.60_2.40_pt7.00_8.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 8.00–9.00 | -0.147–2.002 | 400 | 0.683 | 0.732 | 1.109 | 0.0259 | 282/400 | [plot](fit2d_ctau_y1.60_2.40_pt8.00_9.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 9.00–10.00 | -0.125–2.054 | 400 | 0.637 | 0.601 | 1.121 | 0.025 | 313/400 | [plot](fit2d_ctau_y1.60_2.40_pt9.00_10.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 10.00–12.00 | -0.080–1.728 | 400 | 0.734 | 0.887 | 1.219 | 0.0309 | 301/400 | [plot](fit2d_ctau_y1.60_2.40_pt10.00_12.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 12.00–14.00 | -0.069–1.798 | 400 | 0.646 | 0.896 | 1.002 | 0.0287 | 345/400 | [plot](fit2d_ctau_y1.60_2.40_pt12.00_14.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 14.00–20.00 | -0.074–1.836 | 400 | 0.722 | 1.055 | 1.299 | 0.0408 | 349/400 | [plot](fit2d_ctau_y1.60_2.40_pt14.00_20.00.pdf) |
| fit2d_mass | 0.00–1.60 | 7.00–8.00 | 2.600–3.500 | 400 | 0.725 | 0.635 | 0.884 | 0.00396 | 0/400 | [plot](fit2d_mass_y0.00_1.60_pt7.00_8.00.pdf) |
| fit2d_mass | 0.00–1.60 | 8.00–9.00 | 2.600–3.500 | 400 | 0.703 | 0.821 | 0.943 | 0.00386 | 0/400 | [plot](fit2d_mass_y0.00_1.60_pt8.00_9.00.pdf) |
| fit2d_mass | 0.00–1.60 | 9.00–10.00 | 2.600–3.500 | 400 | 0.671 | 0.891 | 1.053 | 0.00403 | 248/400 | [plot](fit2d_mass_y0.00_1.60_pt9.00_10.00.pdf) |
| fit2d_mass | 0.00–1.60 | 10.00–12.00 | 2.600–3.500 | 400 | 0.612 | 0.715 | 0.927 | 0.00498 | 134/400 | [plot](fit2d_mass_y0.00_1.60_pt10.00_12.00.pdf) |
| fit2d_mass | 0.00–1.60 | 12.00–14.00 | 2.600–3.500 | 400 | 0.545 | 0.686 | 1.046 | 0.00379 | 324/400 | [plot](fit2d_mass_y0.00_1.60_pt12.00_14.00.pdf) |
| fit2d_mass | 0.00–1.60 | 14.00–16.00 | 2.600–3.500 | 400 | 0.538 | 0.736 | 1.181 | 0.00248 | 347/400 | [plot](fit2d_mass_y0.00_1.60_pt14.00_16.00.pdf) |
| fit2d_mass | 0.00–1.60 | 16.00–20.00 | 2.600–3.500 | 400 | 0.501 | 0.527 | 0.833 | 0.00194 | 347/400 | [plot](fit2d_mass_y0.00_1.60_pt16.00_20.00.pdf) |
| fit2d_mass | 0.00–1.60 | 20.00–30.00 | 2.600–3.500 | 400 | 0.446 | 0.505 | 0.718 | 0.00155 | 359/400 | [plot](fit2d_mass_y0.00_1.60_pt20.00_30.00.pdf) |
| fit2d_mass | 1.60–2.40 | 1.00–2.00 | 2.600–3.500 | 1350 | 1.102 | 1.089 | 1.098 | 4.46e-05 | 0/1350 | [plot](fit2d_mass_y1.60_2.40_pt1.00_2.00.pdf) |
| fit2d_mass | 1.60–2.40 | 2.00–3.00 | 2.600–3.500 | 1150 | 1.019 | 0.955 | 1.014 | 0.000103 | 0/1150 | [plot](fit2d_mass_y1.60_2.40_pt2.00_3.00.pdf) |
| fit2d_mass | 1.60–2.40 | 3.00–4.00 | 2.600–3.500 | 950 | 0.978 | 1.003 | 0.987 | 0.000181 | 0/950 | [plot](fit2d_mass_y1.60_2.40_pt3.00_4.00.pdf) |
| fit2d_mass | 1.60–2.40 | 4.00–5.00 | 2.600–3.500 | 700 | 0.955 | 1.005 | 0.982 | 0.000388 | 0/700 | [plot](fit2d_mass_y1.60_2.40_pt4.00_5.00.pdf) |
| fit2d_mass | 1.60–2.40 | 5.00–6.00 | 2.600–3.500 | 450 | 0.955 | 1.182 | 1.021 | 0.000888 | 0/450 | [plot](fit2d_mass_y1.60_2.40_pt5.00_6.00.pdf) |
| fit2d_mass | 1.60–2.40 | 6.00–7.00 | 2.600–3.500 | 400 | 0.849 | 0.932 | 0.963 | 0.000994 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt6.00_7.00.pdf) |
| fit2d_mass | 1.60–2.40 | 7.00–8.00 | 2.600–3.500 | 400 | 0.780 | 1.009 | 0.954 | 0.000813 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt7.00_8.00.pdf) |
| fit2d_mass | 1.60–2.40 | 8.00–9.00 | 2.600–3.500 | 400 | 0.749 | 0.878 | 0.956 | 0.000769 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt8.00_9.00.pdf) |
| fit2d_mass | 1.60–2.40 | 9.00–10.00 | 2.600–3.500 | 400 | 0.620 | 0.837 | 0.938 | 0.000627 | 122/400 | [plot](fit2d_mass_y1.60_2.40_pt9.00_10.00.pdf) |
| fit2d_mass | 1.60–2.40 | 10.00–12.00 | 2.600–3.500 | 400 | 0.629 | 0.837 | 1.005 | 0.000672 | 123/400 | [plot](fit2d_mass_y1.60_2.40_pt10.00_12.00.pdf) |
| fit2d_mass | 1.60–2.40 | 12.00–14.00 | 2.600–3.500 | 400 | 0.427 | 0.634 | 0.916 | 0.000427 | 322/400 | [plot](fit2d_mass_y1.60_2.40_pt12.00_14.00.pdf) |
| fit2d_mass | 1.60–2.40 | 14.00–20.00 | 2.600–3.500 | 400 | 0.485 | 0.659 | 0.934 | 0.000445 | 327/400 | [plot](fit2d_mass_y1.60_2.40_pt14.00_20.00.pdf) |
