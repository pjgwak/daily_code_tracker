# Forward, 1 < pT < 2 GeV/c, 1.6 < |y| < 2.4

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md) · [적용한 정책](../FIT_POLICY.md)

B = **0.095960 ± 0.020464**. 기준값 대비 상대 차이 **+0.13923%**. 오차는 고정 입력과 선택된 모형을 조건으로 한 fit covariance에서 얻었다.

수렴 status=0, covQual=3, EDM=1.13e-07. ctau 적합 범위 [-4.0, 6.5] mm.

## 적용 모형과 매개변수

PR Gaussian 성분 수: 4. NP 지수 성분 수: 3.

PR 평균·폭·구성비를 적합한다. PR 폭 좌표의 Gaussian 제약 표준편차는 (0.0075, 0.015, 0.030, 0.060) mm이며, 구성비에는 MC 공분산을 적용한다. NP3의 공분산 배율은 175이다. 배경 수명·yield는 sideband prefit 제약 아래 적합한다.

## 해석 시 확인할 사항

- NP 혼합모형의 한 성분 비중이 0 경계에 접근한다. 해당 성분의 수명을 독립적으로 결정한 결과로 해석하지 않는다.

[최종 ROOT](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/fit2d_result_y1.60_2.40_pt1.00_2.00.root) · [초기값·범위](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt1.00_2.00/ACTUAL_INITIAL_START_0.csv) · [최종 변수](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt1.00_2.00/FINAL_OBJECTIVE_PARAMETERS.csv) · [실제 제약 항](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt1.00_2.00/EXTERNAL_CONSTRAINTS.txt)

## Prefit 적합도

| Prefit | χ²/ndf | status/covQual | EDM | Equal-population χ²/ndf |
|---|---:|:---:|---:|---:|
| MC_MASS | 0.6818 | 0/3 | 0.000952 | — |
| PR | 0.7366 | 0/3 | 3.56e-05 | — |
| NP | 0.9546 | 0/3 | 9.13e-05 | 0.8688 (ndf 94) |
| MASS | 1.0456 | 0/3 | 5.04e-05 | — |
| BKG | 0.6856 | 0/3 | 0.000348 | 1.0447 (ndf 86) |

원래 저장된 binning의 χ²/ndf와 equal-population 진단을 구별한다. 후자는 추가 binning 진단이며 최소값을 골라 통과 여부를 바꾸지 않았다.

## MC mass prefit

![MC mass prefit](../figures/prefits/fwd_1_2_mc_mass.png)

[PDF](../figures/prefits/fwd_1_2_mc_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y1.60_2.40/mc_mass/mc_mass_model_y1.60_2.40_pt1.00_2.00.root)

## Prompt MC decay-length prefit

![Prompt MC decay-length prefit](../figures/prefits/fwd_1_2_pr.png)

[PDF](../figures/prefits/fwd_1_2_pr.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y1.60_2.40/ctau_pr/ctau_resolution_y1.60_2.40_pt1.00_2.00.root)

## Nonprompt MC truth prefit

![Nonprompt MC truth prefit](../figures/prefits/fwd_1_2_np.png)

[PDF](../figures/prefits/fwd_1_2_np.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_np/roots/y1.60_2.40/ctau_np/ctau_np_model_y1.60_2.40_pt1.00_2.00.root)

## Data mass prefit

![Data mass prefit](../figures/prefits/fwd_1_2_mass.png)

[PDF](../figures/prefits/fwd_1_2_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y1.60_2.40/mass/mass_model_y1.60_2.40_pt1.00_2.00.root)

## Data sideband decay-length prefit

![Data sideband decay-length prefit](../figures/prefits/fwd_1_2_bkg.png)

[PDF](../figures/prefits/fwd_1_2_bkg.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y1.60_2.40/ctau_bkg/ctau_bkg_fitresult_y1.60_2.40_pt1.00_2.00.root)

## Invariant-mass projection

χ²/ndf = 85.848/98 = 0.8760; diagnostic p-value = 0.80471.

![Invariant-mass projection](../figures/projections/fwd_1_2_mass.png)

[PDF](../figures/projections/fwd_1_2_mass.pdf)

## Decay-length projection

χ²/ndf = 294.902/425 = 0.6939; diagnostic p-value = 1.

![Decay-length projection](../figures/projections/fwd_1_2_lifetime.png)

[PDF](../figures/projections/fwd_1_2_lifetime.pdf)

## Data-fit PDF boundary review

회색 점은 전체 mass 범위의 데이터이며 빨강·파랑 곡선은 각각 피팅된 PR·NP 신호 성분이다. 배경이 포함된 데이터와 신호 성분만의 높이를 직접적인 전체 적합도 검사로 해석하지 않는다. 검은 점선은 실제 적합 경계, 초록 점선은 진단 코드가 제안한 경계이며 적용하지 않았다.

![Data-fit PR and NP boundaries](../figures/boundaries/fwd_1_2.png)

[PDF](../figures/boundaries/fwd_1_2.pdf) · [경계 수치](../figures/boundaries/windows.csv)

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md)
