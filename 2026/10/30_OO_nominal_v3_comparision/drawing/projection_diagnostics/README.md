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
| fit2d_ctau | 0.00–1.60 | 7.00–8.00 | -0.107–2.104 | 400 | 0.810 | 1.117 | 1.302 | 0.0575 | 243/400 | [plot](fit2d_ctau_y0.00_1.60_pt7.00_8.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 8.00–9.00 | -0.116–2.114 | 400 | 0.641 | 0.935 | 0.981 | 0.0505 | 261/400 | [plot](fit2d_ctau_y0.00_1.60_pt8.00_9.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 9.00–10.00 | -0.102–1.886 | 400 | 0.679 | 0.772 | 1.126 | 0.0393 | 272/400 | [plot](fit2d_ctau_y0.00_1.60_pt9.00_10.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 10.00–12.00 | -0.078–2.117 | 400 | 0.752 | 0.922 | 1.194 | 0.0852 | 260/400 | [plot](fit2d_ctau_y0.00_1.60_pt10.00_12.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 12.00–14.00 | -0.075–2.052 | 400 | 0.680 | 0.898 | 1.215 | 0.0626 | 292/400 | [plot](fit2d_ctau_y0.00_1.60_pt12.00_14.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 14.00–16.00 | -0.066–1.882 | 400 | 0.666 | 0.874 | 1.093 | 0.0284 | 340/400 | [plot](fit2d_ctau_y0.00_1.60_pt14.00_16.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 16.00–20.00 | -0.056–2.030 | 400 | 0.591 | 0.724 | 1.037 | 0.049 | 338/400 | [plot](fit2d_ctau_y0.00_1.60_pt16.00_20.00.pdf) |
| fit2d_ctau | 0.00–1.60 | 20.00–30.00 | -0.048–2.089 | 400 | 0.582 | 0.557 | 0.884 | 0.0827 | 367/400 | [plot](fit2d_ctau_y0.00_1.60_pt20.00_30.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 1.00–2.00 | -1.160–3.258 | 1350 | 0.857 | 0.990 | 0.990 | 0.00162 | 0/1350 | [plot](fit2d_ctau_y1.60_2.40_pt1.00_2.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 2.00–3.00 | -0.603–2.281 | 1150 | 0.959 | 2.494 | 1.183 | 0.00203 | 168/1150 | [plot](fit2d_ctau_y1.60_2.40_pt2.00_3.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 3.00–4.00 | -0.287–2.331 | 950 | 0.851 | 1.559 | 1.175 | 0.00515 | 305/950 | [plot](fit2d_ctau_y1.60_2.40_pt3.00_4.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 4.00–5.00 | -0.241–2.225 | 700 | 0.641 | 0.980 | 0.899 | 0.0111 | 286/700 | [plot](fit2d_ctau_y1.60_2.40_pt4.00_5.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 5.00–6.00 | -0.181–1.925 | 450 | 0.986 | 1.749 | 1.463 | 0.028 | 199/450 | [plot](fit2d_ctau_y1.60_2.40_pt5.00_6.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 6.00–7.00 | -0.169–1.782 | 400 | 0.820 | 1.071 | 1.269 | 0.0298 | 209/400 | [plot](fit2d_ctau_y1.60_2.40_pt6.00_7.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 7.00–8.00 | -0.154–1.899 | 400 | 0.714 | 0.872 | 1.095 | 0.0299 | 248/400 | [plot](fit2d_ctau_y1.60_2.40_pt7.00_8.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 8.00–9.00 | -0.149–2.029 | 400 | 0.750 | 0.804 | 1.163 | 0.0251 | 283/400 | [plot](fit2d_ctau_y1.60_2.40_pt8.00_9.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 9.00–10.00 | -0.108–2.070 | 400 | 0.610 | 0.565 | 1.000 | 0.0237 | 316/400 | [plot](fit2d_ctau_y1.60_2.40_pt9.00_10.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 10.00–12.00 | -0.116–1.764 | 400 | 0.782 | 0.971 | 1.347 | 0.032 | 291/400 | [plot](fit2d_ctau_y1.60_2.40_pt10.00_12.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 12.00–14.00 | -0.100–1.885 | 400 | 0.686 | 0.832 | 1.016 | 0.0329 | 344/400 | [plot](fit2d_ctau_y1.60_2.40_pt12.00_14.00.pdf) |
| fit2d_ctau | 1.60–2.40 | 14.00–20.00 | -0.086–1.843 | 400 | 0.762 | 1.035 | 1.587 | 0.0394 | 349/400 | [plot](fit2d_ctau_y1.60_2.40_pt14.00_20.00.pdf) |
| fit2d_mass | 0.00–1.60 | 7.00–8.00 | 2.600–3.500 | 400 | 0.727 | 0.642 | 0.882 | 0.00394 | 0/400 | [plot](fit2d_mass_y0.00_1.60_pt7.00_8.00.pdf) |
| fit2d_mass | 0.00–1.60 | 8.00–9.00 | 2.600–3.500 | 400 | 0.703 | 0.822 | 0.943 | 0.00386 | 0/400 | [plot](fit2d_mass_y0.00_1.60_pt8.00_9.00.pdf) |
| fit2d_mass | 0.00–1.60 | 9.00–10.00 | 2.600–3.500 | 400 | 0.670 | 0.892 | 1.045 | 0.00403 | 245/400 | [plot](fit2d_mass_y0.00_1.60_pt9.00_10.00.pdf) |
| fit2d_mass | 0.00–1.60 | 10.00–12.00 | 2.600–3.500 | 400 | 0.612 | 0.718 | 0.924 | 0.00497 | 134/400 | [plot](fit2d_mass_y0.00_1.60_pt10.00_12.00.pdf) |
| fit2d_mass | 0.00–1.60 | 12.00–14.00 | 2.600–3.500 | 400 | 0.544 | 0.686 | 1.036 | 0.00379 | 324/400 | [plot](fit2d_mass_y0.00_1.60_pt12.00_14.00.pdf) |
| fit2d_mass | 0.00–1.60 | 14.00–16.00 | 2.600–3.500 | 400 | 0.536 | 0.735 | 1.201 | 0.00248 | 347/400 | [plot](fit2d_mass_y0.00_1.60_pt14.00_16.00.pdf) |
| fit2d_mass | 0.00–1.60 | 16.00–20.00 | 2.600–3.500 | 400 | 0.499 | 0.526 | 0.832 | 0.00193 | 347/400 | [plot](fit2d_mass_y0.00_1.60_pt16.00_20.00.pdf) |
| fit2d_mass | 0.00–1.60 | 20.00–30.00 | 2.600–3.500 | 400 | 0.446 | 0.504 | 0.717 | 0.00155 | 359/400 | [plot](fit2d_mass_y0.00_1.60_pt20.00_30.00.pdf) |
| fit2d_mass | 1.60–2.40 | 1.00–2.00 | 2.600–3.500 | 1350 | 1.102 | 1.077 | 1.099 | 4.45e-05 | 0/1350 | [plot](fit2d_mass_y1.60_2.40_pt1.00_2.00.pdf) |
| fit2d_mass | 1.60–2.40 | 2.00–3.00 | 2.600–3.500 | 1150 | 1.019 | 0.955 | 1.014 | 0.000103 | 0/1150 | [plot](fit2d_mass_y1.60_2.40_pt2.00_3.00.pdf) |
| fit2d_mass | 1.60–2.40 | 3.00–4.00 | 2.600–3.500 | 950 | 0.982 | 1.034 | 0.990 | 0.000179 | 0/950 | [plot](fit2d_mass_y1.60_2.40_pt3.00_4.00.pdf) |
| fit2d_mass | 1.60–2.40 | 4.00–5.00 | 2.600–3.500 | 700 | 0.954 | 1.001 | 0.981 | 0.00039 | 0/700 | [plot](fit2d_mass_y1.60_2.40_pt4.00_5.00.pdf) |
| fit2d_mass | 1.60–2.40 | 5.00–6.00 | 2.600–3.500 | 450 | 0.955 | 1.182 | 1.021 | 0.000888 | 0/450 | [plot](fit2d_mass_y1.60_2.40_pt5.00_6.00.pdf) |
| fit2d_mass | 1.60–2.40 | 6.00–7.00 | 2.600–3.500 | 400 | 0.849 | 0.932 | 0.963 | 0.000994 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt6.00_7.00.pdf) |
| fit2d_mass | 1.60–2.40 | 7.00–8.00 | 2.600–3.500 | 400 | 0.776 | 0.996 | 0.954 | 0.000819 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt7.00_8.00.pdf) |
| fit2d_mass | 1.60–2.40 | 8.00–9.00 | 2.600–3.500 | 400 | 0.749 | 0.877 | 0.957 | 0.00077 | 0/400 | [plot](fit2d_mass_y1.60_2.40_pt8.00_9.00.pdf) |
| fit2d_mass | 1.60–2.40 | 9.00–10.00 | 2.600–3.500 | 400 | 0.620 | 0.837 | 0.939 | 0.000628 | 122/400 | [plot](fit2d_mass_y1.60_2.40_pt9.00_10.00.pdf) |
| fit2d_mass | 1.60–2.40 | 10.00–12.00 | 2.600–3.500 | 400 | 0.628 | 0.828 | 1.007 | 0.000676 | 123/400 | [plot](fit2d_mass_y1.60_2.40_pt10.00_12.00.pdf) |
| fit2d_mass | 1.60–2.40 | 12.00–14.00 | 2.600–3.500 | 400 | 0.429 | 0.634 | 0.938 | 0.000429 | 322/400 | [plot](fit2d_mass_y1.60_2.40_pt12.00_14.00.pdf) |
| fit2d_mass | 1.60–2.40 | 14.00–20.00 | 2.600–3.500 | 400 | 0.484 | 0.659 | 0.976 | 0.000445 | 327/400 | [plot](fit2d_mass_y1.60_2.40_pt14.00_20.00.pdf) |
