# Separation of prompt and nonprompt J/ψ mesons

아래 영문은 현재 구현과 저장된 결과에 맞춘 논문용 방법 서술 초안이다. 수치 설정은 [구간별 정책표](FIT_POLICY.md), 검증 수치는 [검증 기록](VALIDATION.md), 그림은 [빈별 검토 문서](FIT_PLOTS_BY_BIN.md)에 정리했다.

## Signal extraction

The prompt and nonprompt J/ψ contributions are separated with an extended unbinned maximum-likelihood fit to the dimuon invariant mass, m, and the reconstructed pseudo-proper decay length, ℓ. The fits are performed independently in twelve transverse-momentum intervals in the forward rapidity region, 1.6 < |y| < 2.4, and eight intervals in the midrapidity region, |y| < 1.6. The invariant-mass range is 2.6–3.5 GeV/c². Both muons are required to satisfy |η| < 2.4, with pT > 3.5 GeV/c for |η| < 1.6 and pT > 1.5 GeV/c otherwise. The same muon acceptance requirements are used for data and the simulated samples entering the calibration fits.

The fitted event intensity in each analysis interval is written as

\[
\mathcal I(m,\ell)=N_{\mathrm{sig}}S(m)
\left[(1-f_{\mathrm{NP}})R(\ell)+f_{\mathrm{NP}}F_{\mathrm{NP}}(\ell)\right]
+N_{\mathrm{bkg}}B_m(m)B_\ell(\ell).
\]

Here, S and B_m describe the signal and background mass distributions, R describes the prompt decay-length resolution, and F_NP describes the nonprompt decay-length distribution. All probability density functions are normalized over their respective fit domains. The signal and background yields and the nonprompt fraction f_NP are fitted, with 0 ≤ f_NP ≤ 1. The reported fraction refers to the signal yield in the applied fit domain; no post-fit extrapolation of the fraction is applied. The mass and decay-length dependence is factorized within the signal and background components in the fitted model.

The signal and background mass-shape parameters are obtained from the corresponding calibration fits and held fixed in the two-dimensional fit. The prompt decay-length distribution is described by a Gaussian mixture with a common mean and ordered widths. The independent width coordinates are σ₁ and the positive differences Δᵢ = σᵢ − σᵢ₋₁. The nonprompt model is a mixture of positive-decay-length exponential functions convolved with the same resolution model. The component numbers are selected from the fits to the corresponding simulated samples and need not be identical across analysis intervals. The finite fit domain is taken into account in the normalization and mixture coefficients.

## Constraints and fixed parameters

The objective function includes Gaussian constraint terms for selected calibration parameters:

\[
-\ln\mathcal L_{\mathrm{tot}}
=N_{\mathrm{sig}}+N_{\mathrm{bkg}}
-\sum_{i}\ln\mathcal I(m_i,\ell_i)
+\frac12(\boldsymbol\theta-\boldsymbol\theta_0)^T
V^{-1}(\boldsymbol\theta-\boldsymbol\theta_0)
+Q_{\mathrm{ratio}},
\]

up to terms independent of the fit parameters. The constrained coordinates, central values, covariance matrices, and width-ratio terms are specified for each kinematic region in the configuration table. A constrained parameter remains a fit parameter; a fixed parameter is held at its calibration value throughout the minimization. Where a covariance matrix is rescaled, the scale applies to the covariance matrix, rather than directly to the parameter standard deviations. For the joint prompt constraints, the correlation matrix is taken from the prompt simulation fit, while the applied standard deviations are specified explicitly in the configuration table.

The prompt mixture fractions are fixed to their simulated values in intervals where they are not released by the regional configuration. Selected core-width coordinates are also fixed in the forward region above 10 GeV/c. No nonprompt lifetime parameter is fixed in the selected configuration. For forward 2 < pT < 3 GeV/c, the nonprompt constraint is expressed in terms of the mean lifetime and total normalization obtained from simulation. In the midrapidity interval 8 < pT < 9 GeV/c, the nonprompt lifetimes and mixture fraction are fitted without an additional simulation constraint.

The background decay-length model is initialized using invariant-mass sideband data. For forward 1 < pT < 2 GeV/c and 3 < pT < 5 GeV/c, the background lifetime and component-normalization parameters are fitted with constraints derived from the sideband fit. In the other intervals, the background resolution core is held fixed and the remaining background decay-length parameters are profiled in the two-dimensional fit. The complete list of fixed and floating coordinates is retained with the fit configuration.

## Fit domains and numerical validation

The decay-length fit ranges are selected by inspecting the prompt and nonprompt distributions obtained from the two-dimensional fits to data. The adopted ranges are [−4, 6.5] mm for forward 1 < pT < 3 GeV/c, [−3.5, 4.5] mm for forward 3 < pT < 5 GeV/c, [−1, 4] mm for forward 5 < pT < 10 GeV/c, and [−0.5, 3] mm for forward 10 < pT < 20 GeV/c. The range [−1, 4] mm is used throughout the midrapidity region. The fitted prompt and nonprompt curves decrease below the displayed yield threshold before the boundaries of these ranges.

The calibration fits are assessed using their minimization status, covariance quality, estimated distance to the minimum, and binned goodness-of-fit diagnostics. The final fits are examined using the invariant-mass and decay-length projections. The corresponding χ² values, numbers of degrees of freedom, and diagnostic p-values are reported for each analysis interval. These projection diagnostics are not interpreted as a single goodness-of-fit p-value for the full two-dimensional distribution. Numerical convergence and the agreement with reference fit values are recorded separately.

The selected configuration was rerun over all twenty analysis intervals using the recorded initial values, initial parameter errors, and parameter bounds. The fitted prompt and nonprompt probability densities were reproduced from the saved fit results and compared with an independent numerical evaluation. The retained result has a valid minimization status and covariance matrix in every interval.

## Fit uncertainties and component interpretation

The quoted fit errors are obtained from the covariance matrix conditional on the selected model, constraints, and fixed calibration parameters. Contributions associated with fixed calibration inputs, model choices, and constraint variations require separate uncertainty evaluation.

In the forward interval 2 < pT < 3 GeV/c, the fitted nonprompt lifetimes converge to approximately 0.42418 mm. A nonprompt component has a negligible fitted fraction in forward 1 < pT < 2 GeV/c and midrapidity 8 < pT < 9 GeV/c. The individual exponential components are therefore not interpreted as separately measured physical lifetimes. Reproduction of the selected forward 2–3 and midrapidity 8–9 results requires the recorded initialization; convergence to a unique solution from arbitrary starting points is not assumed.

## Suggested figure captions

**Invariant-mass and decay-length projections.** Projections of the simultaneous invariant-mass and decay-length fit in the specified transverse-momentum and rapidity interval. The curves show the total model and its components, with the prompt and nonprompt contributions displayed separately in the decay-length projection. The lower panels show the residuals normalized by their uncertainties. The quoted χ²/ndf values refer to the individual projections.

**Calibration fits.** Distributions used to determine the invariant-mass and decay-length calibration parameters. The prompt and nonprompt signal calibrations are obtained from simulation, while the background decay-length calibration is obtained from invariant-mass sideband data. The displayed goodness-of-fit quantities refer to the stated binning and calibration sample.

**Fit-domain inspection.** Prompt and nonprompt decay-length distributions evaluated with the parameters obtained from the two-dimensional fit to data. The vertical black lines indicate the adopted fit boundaries. The data points include all invariant masses in the fit interval; the displayed signal components alone are not the total model for these points. The green lines are diagnostic range suggestions and were not used as the fit boundaries.

**Reference comparison.** Relative differences between the fitted nonprompt fractions and the reference values in each kinematic interval. The shaded band denotes the ±2.1% reference-comparison interval used in the configuration study. This comparison is distinct from the goodness-of-fit evaluation of the data distributions.
