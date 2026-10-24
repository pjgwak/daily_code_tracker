# Mid, 7 < pT < 8 GeV/c, |y| < 1.6

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md) · [적용한 정책](../FIT_POLICY.md)

B = **0.236488 ± 0.010234**. 기준값 대비 상대 차이 **-0.82269%**. 오차는 고정 입력과 선택된 모형을 조건으로 한 fit covariance에서 얻었다.

수렴 status=0, covQual=3, EDM=1.39e-08. ctau 적합 범위 [-1.0, 4.0] mm.

## 적용 모형과 매개변수

PR Gaussian 성분 수: 4. NP 지수 성분 수: 2.

PR 구성비를 MC 값에 고정하고 평균·폭을 추가 Gaussian 제약 없이 적합한다. NP에는 MC 공분산을 적용하되 최장 수명 좌표는 제약에서 제외한다.

[최종 ROOT](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y0.00_1.60/fit2d/fit2d_result_y0.00_1.60_pt7.00_8.00.root) · [초기값·범위](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y0.00_1.60/fit2d/diagnostics_y0.00_1.60_pt7.00_8.00/ACTUAL_INITIAL_START_0.tsv) · [최종 변수](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y0.00_1.60/fit2d/diagnostics_y0.00_1.60_pt7.00_8.00/FINAL_OBJECTIVE_PARAMETERS.tsv) · [실제 제약 항](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y0.00_1.60/fit2d/diagnostics_y0.00_1.60_pt7.00_8.00/EXTERNAL_CONSTRAINTS.txt)

## Prefit 적합도

| Prefit | χ²/ndf | status/covQual | EDM | Equal-population χ²/ndf |
|---|---:|:---:|---:|---:|
| MC_MASS | 1.0540 | 0/3 | 0.00016 | — |
| PR | 0.9573 | 0/3 | 0.000576 | — |
| NP | 0.9232 | 0/3 | 0.000448 | 0.8202 (ndf 96) |
| MASS | 1.0249 | 0/3 | 0.000624 | — |
| BKG | 1.4087 | 0/3 | 1.46e-05 | 4.1594 (ndf 21) |

원래 저장된 binning의 χ²/ndf와 equal-population 진단을 구별한다. 후자는 추가 binning 진단이며 최소값을 골라 통과 여부를 바꾸지 않았다.

## MC mass prefit

![MC mass prefit](../figures/prefits/mid_7_8_mc_mass.png)

[PDF](../figures/prefits/mid_7_8_mc_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y0.00_1.60/mc_mass/mc_mass_model_y0.00_1.60_pt7.00_8.00.root)

## Prompt MC decay-length prefit

![Prompt MC decay-length prefit](../figures/prefits/mid_7_8_pr.png)

[PDF](../figures/prefits/mid_7_8_pr.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y0.00_1.60/ctau_pr/ctau_resolution_y0.00_1.60_pt7.00_8.00.root)

## Nonprompt MC truth prefit

![Nonprompt MC truth prefit](../figures/prefits/mid_7_8_np.png)

[PDF](../figures/prefits/mid_7_8_np.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_np/roots/y0.00_1.60/ctau_np/ctau_np_model_y0.00_1.60_pt7.00_8.00.root)

## Data mass prefit

![Data mass prefit](../figures/prefits/mid_7_8_mass.png)

[PDF](../figures/prefits/mid_7_8_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y0.00_1.60/mass/mass_model_y0.00_1.60_pt7.00_8.00.root)

## Data sideband decay-length prefit

![Data sideband decay-length prefit](../figures/prefits/mid_7_8_bkg.png)

[PDF](../figures/prefits/mid_7_8_bkg.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y0.00_1.60/ctau_bkg/ctau_bkg_fitresult_y0.00_1.60_pt7.00_8.00.root)

## Invariant-mass projection

χ²/ndf = 83.111/98 = 0.8481; diagnostic p-value = 0.85869.

![Invariant-mass projection](../figures/projections/mid_7_8_mass.png)

[PDF](../figures/projections/mid_7_8_mass.pdf)

## Decay-length projection

χ²/ndf = 127.632/185 = 0.6899; diagnostic p-value = 0.99956.

![Decay-length projection](../figures/projections/mid_7_8_lifetime.png)

[PDF](../figures/projections/mid_7_8_lifetime.pdf)

## Data-fit PDF boundary review

회색 점은 전체 mass 범위의 데이터이며 빨강·파랑 곡선은 각각 피팅된 PR·NP 신호 성분이다. 배경이 포함된 데이터와 신호 성분만의 높이를 직접적인 전체 적합도 검사로 해석하지 않는다. 검은 점선은 실제 적합 경계, 초록 점선은 진단 코드가 제안한 경계이며 적용하지 않았다.

![Data-fit PR and NP boundaries](../figures/boundaries/mid_7_8.png)

[PDF](../figures/boundaries/mid_7_8.pdf) · [경계 수치](../figures/boundaries/windows.tsv)

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md)
